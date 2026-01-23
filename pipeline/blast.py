import os
import pandas as pd
import subprocess
import time
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML, NCBIWWW

#Almacenamiento de los datos registrados.
cache_tax = {}

#Se crean funciones para obtener la taxonomía de las muestras introducidas.
def taxid_from_acc(accession):
    try:
        handle = Entrez.elink(
            dbfrom="nucleotide",
            db="taxonomy",
            id=accession
        )
        record = Entrez.read(handle)
        handle.close()

        time.sleep(0.34)

        links = record[0].get("LinkSetDb", [])
        if not links:
            return None

        return links[0]["Link"][0]["Id"]

    except Exception:
        return None

def taxon_from_taxid(taxid):
    lineage = {
        "Domain": None,
        "Phylum": None,
        "Class": None,
        "Order": None,
        "Family": None,
        "Genus": None,
        "Species": None}

    try:
        handle = Entrez.efetch(
            db="taxonomy",
            id=taxid,
            retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()

        time.sleep(0.4)
    except Exception:
        return lineage

    tax_rec = record[0]

    for node in tax_rec.get("LineageEx", []):
        rank = node["Rank"].capitalize()
        if rank in lineage:
            lineage[rank] = node["ScientificName"]

    if tax_rec.get("Rank") == "species":
        lineage["Species"] = tax_rec["ScientificName"]

    return lineage

#Se construye una función que almacena la omía de las muestras y la almacena en forma de dataframe.
def taxonomia(fasta, email):
    Entrez.email = email
    cache_tax.clear()
    filas = []

    for record in SeqIO.parse(fasta, "fasta"):
        try:
            #BLAST remoto
            result_handle = NCBIWWW.qblast(program="blastn", database="nt",
                sequence=record.seq, hitlist_size=1, expect=1e-10)

            blast_record = NCBIXML.read(result_handle)

            #Sin hits BLAST
            if not blast_record.alignments:
                filas.append({
                    "ID": record.id,
                    "Dominio": None,
                    "Filo": None,
                    "Clase": None,
                    "Orden": None,
                    "Familia": None,
                    "Género": None,
                    "Especie": None,
                    "ID hit": None,
                    "Def ID": "Sin hits BLAST",
                    "Longitud alineamiento": 0,
                    "Bases idénticas": 0,
                    "% Identidad": 0,
                    "Gaps": None,
                    "e Valor": None,
                    "bit Score": None,
                    "Cobertura": 0
                })
                continue

            #Mejor hit
            aln = blast_record.alignments[0]
            hsp = aln.hsps[0]
            hit_acc = aln.accession.split(".")[0]

            #Taxonomía
            if hit_acc not in cache_tax:
                taxid = taxid_from_acc(hit_acc)
                if taxid:
                    cache_tax[hit_acc] = taxon_from_taxid(taxid)
                else:
                    cache_tax[hit_acc] = {r: None for r in
                        ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]}

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
                "ID hit": aln.accession,
                "Def ID": aln.hit_def,
                "Longitud alineamiento": hsp.align_length,
                "Bases idénticas": hsp.identities,
                "% Identidad": round(hsp.identities / hsp.align_length * 100, 2),
                "Gaps": hsp.gaps,
                "e Valor": hsp.expect,
                "bit Score": hsp.bits,
                "Cobertura": round(cobertura, 2)
            })

        except Exception as e:
            filas.append({
                "ID": record.id,
                "Dominio": None,
                "Filo": None,
                "Clase": None,
                "Orden": None,
                "Familia": None,
                "Género": None,
                "Especie": None,
                "ID hit": None,
                "Def ID": f"Error BLAST: {e}",
                "Longitud alineamiento": None,
                "Bases idénticas": None,
                "% Identidad": None,
                "Gaps": None,
                "e Valor": None,
                "bit Score": None,
                "Cobertura": None
            })

    return pd.DataFrame(filas)

