import streamlit as st
import zipfile
import os
import tempfile
from io import BytesIO
import pandas as pd
from pipeline.carga_ab1 import cargar_ab1_zip
from pipeline.cortar_ab1 import cortar_ab1
from pipeline.QC import qc_plots
from pipeline.consenso import generar_consensos, fasta_consenso
from pipeline.blast import taxonomia

#Se establecen las características de la app.
st.set_page_config(
    page_title="Analisis de ab1",
    layout="wide"
)

st.title("Preprocesado de ficheros ab1")
st.write("Aplicación para el análisis de secuencias Sanger (.ab1)")

#Se cargan los ficheros a la app.
carga_zip = st.file_uploader(
    "Carga el archivo ZIP que contiene los ficheros .ab1",
    type=["zip"])

if carga_zip:
    st.success(f"✅ ZIP {carga_zip.name} cargado correctamente")
    zip_name = os.path.splitext(carga_zip.name)[0]

    #Se almacenan temporalmente los ficheros.
    with tempfile.TemporaryDirectory() as tmpdir:
        ruta_zip = os.path.join(tmpdir, carga_zip.name)

        with open(ruta_zip, "wb") as f:
            f.write(carga_zip.getbuffer())

        try:
            with zipfile.ZipFile(ruta_zip, "r") as zip_ref:
                zip_ref.extractall(tmpdir)

        except zipfile.BadZipFile:
            st.error("El fichero no es archivo ZIP válido")
            st.stop

        #Se identifican los ficheros encontrados.
        ab1_files = []
        for root, dirs, files in os.walk(tmpdir):
            for f in files:
                if f.lower().endswith(".ab1"):
                    ab1_files.append(os.path.join(root, f))        
        
        st.write(f"Se han encontrado {len(ab1_files)} ficheros ab1")

        #Se procesan los datos cargados.
        if not ab1_files:
            st.warning("⚠️No se han encontrado ficheros ab1 en el ZIP")

        else:
            #Se permite establecer el umbral de calidad deseado.
            umbral_usuario = st.slider("Selecciona el umbral de calidad Phred",
                min_value=0, max_value=40, value=20, step=1)
            st.session_state["umbral"] = umbral_usuario

            try:
                #Se extrae la información de los ficheros ab1.
                df = cargar_ab1_zip(ab1_files, umbral = umbral_usuario)
                st.session_state["df"] = df
                st.success("✅ Los datos se han cargado correctamente")
                if "df" in st.session_state:
                    mostrar_df = st.checkbox("📋 Mostrar dataframe AB1")
                    if mostrar_df:
                        st.dataframe(st.session_state["df"])
                        mostrar_plot = st.checkbox("📊 Mostrar gráficas AB1")

                        if mostrar_plot:
                            figs = qc_plots(
                                st.session_state["df"],
                                umbral=st.session_state["umbral"],
                                trimmed=False)

                            for fig in figs:
                                st.pyplot(fig)

                        #Se permite la descarga de los datos en forma de excel.
                        output = BytesIO()
                        with pd.ExcelWriter(output, engine='openpyxl') as writer:
                            df.to_excel(writer, index=False, sheet_name="AB1 Results")
                        processed_data = output.getvalue()

                        st.download_button(label="📥 Descargar resultados Excel",
                        data=processed_data, file_name=f"{zip_name}_ab1.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

            except ValueError as e:
                st.error(str(e))

#Se incluye el trimming de secuencias.
if "df" in st.session_state:
    st.subheader("Limpieza de secuencias")
    
    try:
        trimmed_df = cortar_ab1(df, umbral = umbral_usuario)
        st.session_state["trimmed_df"] = trimmed_df
        st.success("✅ Las secuencias han sido cortadas correctamente")

        if "trimmed_df" in st.session_state:
            mostrar_trimming = st.checkbox("📋 Mostrar secuencias recortadas")
            if mostrar_trimming:
                st.dataframe(st.session_state["trimmed_df"])
                mostrar_plot_trimming = st.checkbox("📊 Mostrar gráficas recortadas")
                
                if mostrar_plot_trimming:
                    figs_trimmed = qc_plots(
                        st.session_state["trimmed_df"],
                        umbral=st.session_state["umbral"],
                        trimmed=True)

                    for fig_trim in figs_trimmed:
                        st.pyplot(fig_trim)
                

                #Se incluye el botón de descarga.
                output_trim = BytesIO()
                with pd.ExcelWriter(output_trim, engine='openpyxl') as writer:
                    trimmed_df.to_excel(writer, index=False, sheet_name="Trimmed Results")
                processed_trim = output_trim.getvalue()

                st.download_button(label="📥 Descargar resultados del trimming",
                    data=processed_trim, file_name=f"{zip_name}_ab1_trimmed.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    except ValueError as e:
        st.error(str(e)) 

#Se incluye la construcción de secuencias consenso.
if "trimmed_df" in st.session_state:
    st.subheader("Secuencias consenso")

    try:
        cons_df = generar_consensos(trimmed_df, umbral=umbral_usuario)
        st.session_state["cons_df"] = cons_df
        
        if not st.session_state["cons_df"].empty:
            st.session_state["cons_fasta"] = fasta_consenso(cons_df)
            st.success("✅ Las secuencias consenso se han generado correctamente")

        if "cons_df" in st.session_state:
            mostrar_consenso = st.checkbox("📋 Mostrar secuencias consenso")
            if mostrar_consenso:
                st.write(f"Se han generado {len(cons_df)} secuencias consenso")
                st.dataframe(st.session_state["cons_df"])

                output_cons = BytesIO()
                with pd.ExcelWriter(output_cons, engine='openpyxl') as writer:
                    cons_df.to_excel(writer, index=False, sheet_name="Trimmed Results")
                processed_cons = output_cons.getvalue()

                st.download_button(label="📥 Descargar resultados del consenso",
                    data=processed_cons, file_name=f"{zip_name}_consensos.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    
                st.download_button(label="📥 Descargar FASTA multiconsenso",
                                   data=st.session_state["cons_fasta"].encode(), file_name=f"{zip_name}_multiconsenso.fasta", 
                                   mime="text/plain")
               
        else:
            st.warning("⚠️ No se han generado consensos válidos")
                
    except ValueError as e:
        st.error(str(e)) 

st.markdown("---")
st.title("Alineamiento de secuencias")
usar_cons = st.checkbox("Usar secuencias preprocesadas", value = True)

fasta_path = None

#Se permite el uso de las secuencias consenso definidas previamente:
if usar_cons:
    if "cons_df" not in st.session_state or st.session_state["cons_df"].empty:
        st.warning("⚠️ No hay secuencias consenso disponibles")
    else:
        fasta_str = st.session_state["cons_fasta"]
        
        with tempfile.NamedTemporaryFile(mode = "w", suffix = ".fasta", delete = False) as tmp_fasta:
            tmp_fasta.write(fasta_str)
            fasta_path = tmp_fasta.name     

#En caso de no activar la checkbox el usuario debe cargar su propio fasta.
else:
    upload_fasta = st.file_uploader("Carga el fichero fasta con secuencias consenso", type = ["fasta", "fa", "fna"])

    if upload_fasta:
        with tempfile.NamedTemporaryFile(mode = "w", suffix = ".fasta", delete = False) as tmp_fasta:
            tmp_fasta.write(upload_fasta.getbuffer())
            fasta_path = tmp_fasta.name

#Se pide el email para usar NCBI.
email = st.text_input("Introduce el email para NCBI Entrez", value = "")    

if fasta_path and st.button("🔬 Alineamiento contra NCBI"):
    if email:
        with st.spinner("Ejecutando BLAST remoto contra NCBI..."):
            try:
                blast_df = taxonomia(fasta_path, email)
                st.session_state["blast_df"] = blast_df
                st.success("✅ Se ha realizado el BLAST correctamente")                    
            
            except Exception as e:
                st.error(f"❌ Error al ejecutar BLAST: {e}")        

if  "blast_df" in st.session_state and not st.session_state["blast_df"].empty:
    mostrar_blast = st.checkbox("📋 Mostrar resultado de BLAST contra NCBI")
    
    if mostrar_blast:
        st.dataframe(st.session_state["blast_df"])
        
        output_blast = BytesIO()
        with pd.ExcelWriter(output_blast, engine='openpyxl') as writer:
            st.session_state["blast_df"].to_excel(writer, index=False, sheet_name="Trimmed Results")
        processed_blast = output_blast.getvalue()
            
        st.download_button(label="📥 Descargar resultados de BLAST",
                            data=processed_blast, file_name=f"Resultado_blast_{zip_name}.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
