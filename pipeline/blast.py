import os
import pandas as pd
import subprocess
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML, NCBIWWW

#Almacenamiento de los datos registrados.
cache_tax = {}

#Se crean funciones para obtener la taxonomía de las muestras introducidas.
def taxid_from_acc(accession):
  try:
    handle = Entrez.esummary(db = "nucleotide", id = accession, retmode = "xml")
    record = Entrez.read(handle)
    handle.close()

    return record[0].get("TaxID")

  except Exception:
    return None

def taxon_from_taxid(taxid):
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

#Se construye una función que almacena la omía de las muestras y la almacena en forma de dataframe.
def taxonomia(fasta, email):
    Entrez.email = email
    
    cache_tax.clear()
    filas = []

    for record in SeqIO.parse(fasta, "fasta"):
        try:
            # Ejecutar BLAST remoto
            result_handle = NCBIWWW.qblast("blastn", "nt", record.seq, hitlist_size=1)
            blast_record = NCBIXML.read(result_handle)

            # Procesar resultados
            for aln in blast_record.alignments:
                hsp = aln.hsps[0]
                hit_acc = aln.accession

                if hit_acc not in cache_tax:
                    taxid = taxid_from_acc(hit_acc, email)
                    if taxid:
                        cache_tax[hit_acc] = taxon_from_taxid(taxid, email)
                    else:
                        cache_tax[hit_acc] = {rango: None for rango in ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]}

                taxon = cache_tax[hit_acc]
                cobertura = hsp.align_length / len(record.seq) * 100

                filas.append({
                    "ID": record.id,
                    "Dominio": taxon.get("Domain"),
                    "Filo": taxon.get("Phylum"),
                    "Clase": taxon.get("Class"),
                    "Orden": taxon.get("Order"),
                    "Familia": taxon.get("Family"),
                    "Género": taxon.get("Genus"),
                    "Especie": taxon.get("Species"),
                    "ID hit": aln.hit_def.split()[0],
                    "Def ID": aln.hit_def,
                    "Longitud alineamiento": hsp.align_length,
                    "Bases idénticas": hsp.identities,
                    "% Identidad": hsp.identities / hsp.align_length * 100,
                    "Gaps": hsp.gaps,
                    "e Valor": hsp.expect,
                    "bit Score": hsp.bits,
                    "Cobertura": cobertura
                })

        except Exception as e:
            print(f"Error BLAST para {record.id}: {e}")
            filas.append({
                "ID": record.id,
                "Dominio": None, "Filo": None, "Clase": None, "Orden": None,
                "Familia": None, "Género": None, "Especie": None,
                "ID hit": None, "Def ID": None, "Longitud alineamiento": None,
                "Bases idénticas": None, "% Identidad": None, "Gaps": None,
                "e Valor": None, "bit Score": None, "Cobertura": None
            })

    df = pd.DataFrame(filas)
    return df
