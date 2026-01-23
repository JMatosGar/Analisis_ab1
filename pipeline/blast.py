import os
import pandas as pd
import time
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML, NCBIWWW

# Cachés para acelerar la búsqueda de taxonomía
taxonomy_cache = {}  # por accession
seq_cache = {}       # por secuencia completa

def taxid_from_acc(accession):
    """Obtener TaxID a partir de un accession NCBI"""
    try:
        handle = Entrez.elink(dbfrom="nucleotide", db="taxonomy", id=accession)
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.34)  # evitar sobrecarga NCBI

        links = record[0].get("LinkSetDb", [])
        if not links:
            return None
        return links[0]["Link"][0]["Id"]
    except Exception:
        return None

def taxon_from_taxid(taxid):
    """Obtener taxonomía completa a partir de TaxID"""
    lineage = {
        "Domain": None, "Phylum": None, "Class": None, "Order": None,
        "Family": None, "Genus": None, "Species": None
    }
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.4)  # evitar sobrecarga NCBI
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

def taxonomia(fasta_path, email):
    """
    Ejecuta BLAST remoto y devuelve un DataFrame con taxonomía.
    Incluye columna FASTA y caché de hits repetidos.
    """
    Entrez.email = email
    taxonomy_cache.clear()
    seq_cache.clear()
    filas = []

    records = list(SeqIO.parse(fasta_path, "fasta"))

    for record in records:
        query_seq = str(record.seq)
        
        # Si la secuencia ya fue procesada
        if query_seq in seq_cache:
            filas.append(seq_cache[query_seq])
            continue

        try:
            # BLAST remoto
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=query_seq,
                hitlist_size=1,
                expect=1e-10
            )
            blast_record = NCBIXML.read(result_handle)

            # Sin hits BLAST
            if not blast_record.alignments:
                row = {
                    "ID": record.id,
                    "FASTA": query_seq,
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
                }
                cached_row = seq_cache[query_seq].copy()
                cached_row["FASTA"] = query_seq
                filas.append(cached_row)
                continue

            # Mejor hit
            aln = blast_record.alignments[0]
            hsp = aln.hsps[0]
            hit_acc = aln.accession.split(".")[0]

            # Taxonomía: usar caché si ya se obtuvo antes
            if hit_acc in taxonomy_cache:
                taxon = taxonomy_cache[hit_acc]
            else:
                taxid = taxid_from_acc(hit_acc)
                taxon = taxon_from_taxid(taxid) if taxid else {
                    r: None for r in ["Domain","Phylum","Class","Order","Family","Genus","Species"]
                }
                taxonomy_cache[hit_acc] = taxon

            cobertura = hsp.align_length / len(record.seq) * 100

            row = {
                "ID": record.id,
                "FASTA": f"{query_seq}\n{record.id}",
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
            }

            filas.append(row)
            row["FASTA"] = query_seq
            seq_cache[query_seq] = row
            
        except Exception as e:
            row = {
                "ID": record.id,
                "FASTA": query_seq,
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
            }
            filas.append(row)
            row["FASTA"] = query_seq
            seq_cache[query_seq] = row


    return pd.DataFrame(filas)
