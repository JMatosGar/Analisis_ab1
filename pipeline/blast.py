import os
import pandas as pd
import subprocess
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML

#Almacenamiento de los datos registrados.
cache_tax = {}

#Se crean funciones para obtener la taxonomía de las muestras introducidas.
def taxid_from_acc(accesion):
  try:
    handle = Entrez.esummary(db = "nucleotide", id = accession, retmode = "xml")
    record = Entrez.read(handle)
    handle.close()

    return record[0].get("TaxID")

  except Exception:
    return None

def taxon_from_taxid(taxid):
  try:
    handle = Entrez.esummary(db = "taxonomy", id = taxid, retmode = "xml")
    record = Entrez.read(handle)
    handle.close()  

    lineage = {rango: None for rango in["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]}
    if not record:
      return lineage

    tax_rec = record[0]
    for nodo in tax_rec.get("LineageEx", []):
      rango = nodo["Rank"].capitalize()
      if rango in lineage:
          lineage[rango] = nodo["ScientificName"]

    if tax_rec["Rank"] == "species":
      lineage["Species"] = tax_rec["ScientificName"]

    return lineage

#Se construye una función que almacena la taxonomía de las muestras y la almacena en forma de dataframe.
def taxonomia(fasta, email, output = None):
  #Se introduce el correo.
  Entrez.email = email

  #Se limpian los datos almacenados previamente.
  cache_tax.clear()

  if output is None:
    output = os.path.join(os.path.dirname(fasta), "resultado BLAST")
  os.makedirs(output, exist_ok = True)
  
  filas = []
  for record in SeqIO.parse(fasta, "fasta"):
    query_name = record.id
    query_fasta = os.path.join(output, f"{query_name}.fasta")
    SeqIO.write(record, query_fasta, "fasta")

    xml_path = os.path.join(output, f"{query_name}.blast.xml")

    cmd = ["blastn", "-query", query_fasta, "-db", "nt", 
           "-remote", "-outfmt", "5", "-max_target_seqs", 
           "1", "-task", "megablast", "-out", xml_path]

    #Se ejecuta el blast remoto.
    try:
      subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
      st.error(f"❌ Fallo al ejecutar BLAST remoto: {e}")
      return None
    except Exception as e:
      st.error(f"❌ Error inesperado: {e}")
      return None

    with open(xml_path) as handle:
      blast_record = NCBIXML.read(handle)
      for aln in blast_record.alignments:
        hsp = aln.hsps[0]

        hit_acc = aln.accession
        
        if hit_acc not in cache_tax:
          taxid = taxid_from_acc(hit_acc)
          
          if taxid:
            cache_tax[hit_acc] = taxon_from_taxid(taxid)
          else:
            cache_tax[hit_acc] = {rango: None for rango in["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]}

        taxonomia = cache_tax[hit_acc]
        cobertura = hsp.align_length / len(record.seq) * 100

        filas.append({
          "ID": query_name,
          "Dominio": taxonomia.get("Domain"),
          "Filo": taxonomia.get("Phylum"),
          "Clase": taxonomia.get("Class"),
          "Orden": taxonomia.get("Order"),
          "Familia": taxonomia.get("Family"),
          "Género": taxonomia.get("Genus"),
          "Especie": taxonomia.get("Species"),
          "ID hit": aln.hit_def.split()[0],
          "Def ID": aln.hit_def,
          "Longitud alineamiento": hsp.align_length,
          "Bases idénticas": hsp.identities,
          "% Identidad": hsp.identities / hsp.align_length * 100,
          "Gaps": hsp.gaps,
          "e Valor": hsp.expect,
          "bit Score": hsp.bits,
          "Cobertura": cobertura})

  df = pd.DataFrame(filas)
  
  return df
