import streamlit as st

st.set_page_config(
    page_title="Analisis de ab1",
    layout="wide"
)

st.title("Preprocesado de ficheros ab1")
st.write("Aplicación para el análisis de secuencias Sanger (.ab1)")

carga_zip = st.file_uploader(
    "Carga el archivo ZIP que contiene los ficheros .ab1",
    type=["zip"]
)

if carga_zip is not None:
    st.success(f"Archivo cargado: {uploaded_zip.name}")
