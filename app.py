import streamlit as st
import zipfile
import os
import tempfile
from io import BytesIO
import pandas as pd
from pipeline.carga_ab1 import cargar_ab1_zip
from pipeline.cortar_ab1 import cortar_ab1
from pipeline.QC import qc_plots

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
            
            try:
                #Se extrae la información de los ficheros ab1.
                df = cargar_ab1_zip(ab1_files, umbral = umbral_usuario)
                st.session_state["df"] = df
                st.session_state["umbral"] = umbral_usuario
                st.success("✅ Los datos se han cargado correctamente")
                if "df" in st.session_state:
                    mostrar_df = st.checkbox("📋 Mostrar dataframe AB1")
                    if mostrar_df:
                        st.dataframe(st.session_state["df"])
                        mostrar_plot = st.checkbox("📊 Mostrar gráficas AB1")

                        if mostrar_plot:
                            figs = plot_qc_metrics(
                                st.session_state["df"],
                                threshold=st.session_state["umbral"],
                                trimmed=False)

                            for fig in figs:
                                st.pyplot(fig)

                        #Se permite la descarga de los datos en forma de excel.
                        output = BytesIO()
                        with pd.ExcelWriter(output, engine='openpyxl') as writer:
                            df.to_excel(writer, index=False, sheet_name="AB1 Results")
                        processed_data = output.getvalue()

                        st.download_button(
                        label="📥 Descargar resultados Excel",
                        data=processed_data,
                        file_name="AB1_results.xlsx",
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
                
                if mostrar_plot:
                    figs_trimmed = plot_qc_metrics(
                        st.session_state["trimmed_df"],
                        threshold=st.session_state["umbral"],
                        trimmed=True)

                    for fig_trim in figs_trimmed:
                        st.pyplot(fig_trim)
                

                #Se incluye el botón de descarga.
                output_trim = BytesIO()
                with pd.ExcelWriter(output_trim, engine='openpyxl') as writer:
                    trimmed_df.to_excel(writer, index=False, sheet_name="Trimmed Results")
                processed_trim = output_trim.getvalue()

                st.download_button(
                    label="📥 Descargar resultados del trimming Excel",
                    data=processed_trim,
                    file_name="AB1_trimmed_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    except ValueError as e:
        st.error(str(e)) 
