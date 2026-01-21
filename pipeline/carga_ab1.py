from Bio import SeqIO
import pandas as pd
import statistics
import os

def cargar_ab1_zip(filepaths, umbral=20):

  datos = []

  if not filepaths:
    raise ValueError("No se han proporcionado ficheros ab1")

  for filepath in filepaths:
    try:
      record = SeqIO.read(filepath, "abi")
    except Exception as e:
      print(f"⚠️No se ha podido leer{filepath}: {e}")

    secuencia = str(record.seq)
    sec_length = len(secuencia)
    calidad = record.letter_annotations.get("phred_quality", [])

#Media y mediana de Phred.
    calidad_media = statistics.mean(calidad) if calidad else 0
    calidad_mediana = statistics.median(calidad) if calidad else 0

#% bases con Phred >= umbral
    calidad_alta = sum(1 for q in calidad if q >= umbral)
    pct_calidad_alta = (calidad_alta / sec_length * 100) if sec_length > 0 else 0

#nº de bases ambiguas (N)
    num_amb = secuencia.upper().count("N")

#Tramo continuo con bases de alta calidad más largo.
    larga_calidad = 0
    actual = 0
    for q in calidad:
      if q >= 20:
        actual += 1
        larga_calidad = max(larga_calidad, actual)
      else:
        actual = 0

    datos.append({
      "ID": os.path.basename(filepath),
      "Secuencia Fasta": secuencia,
      "Longitud secuencia": sec_length,
      "Puntuación calidad": calidad,
      "Calidad media": calidad_media,
      "Mediana calidad": calidad_mediana,
      "% Alta calidad": pct_calidad_alta,
      "Bases ambiguas(N)": num_amb,
      "Tramo de alta calidad": larga_calidad})
    
  if not datos:
    raise ValueError("No se han encontrado ficheros ab1 válidos")
      
  return pd.DataFrame(datos)
