import os
import pandas as pd
import time
import tempfile
import subprocess
from collections import OrderedDict
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
        fasta_str = f">{record.id}\n{query_seq}"
        
        #Si la secuencia ya fue procesada
        if query_seq in seq_cache:
            filas.append(seq_cache[query_seq])
            continue

        try:
            #BLAST remoto
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
                    "FASTA": fasta_str,
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
                "FASTA": fasta_str,
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
                "FASTA": fasta_str,
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

def taxonomia_local(fasta_muestra, fasta_db, max_hits = 25):
    filas = []
    sec_fasta = OrderedDict()

    #Se establece una base de datos temporal.
    db_temp = tempfile.NamedTemporaryFile(delete = False).name
    subprocess.run(["makeblastdb",
                    "-in", fasta_db,
                    "-dbtype", "nucl",
                    "-out", db_temp], check = True)
    
    #Se ejecuta el BLAST local.
    blast_out = tempfile.NamedTemporaryFile(delete = False, suffix = ".xml").name
    subprocess.run(["blastn",
                    "-query", fasta_muestra,
                    "-db", db_temp,
                    "-out", blast_out,
                    "-outfmt", "5",
                    "-max_target_seqs", str(max_hits),
                    "-evalue", "1e-10"], check = True)
    
    #Se cargan los queries y se adicionan al fasta final.
    query_rec = {rec.id: rec for rec in SeqIO.parse(fasta_muestra, "fasta")}

    for qid, rec in query_rec.items():
        sec_fasta[f"QUERY_{qid}"] = str(rec.seq)

    #Se parsean los resultados del BLAST.
    with open(blast_out) as handle:
        for blast_rec in NCBIXML.parse(handle):
            id_rec = blast_rec.query

            if not blast_rec.alignments:
                continue

            for aln in blast_rec.alignments[:max_hits]:
                hsp = aln.hsps[0]
                id_hit = aln.hit_id.split()[0]
                def_hit = aln.hit_def

                identidad = round(hsp.identities/hsp.align_length*100, 4)
                cobertura = round(hsp.align_length/blast_rec.query_length*100,4)

                filas.append({
                    "ID": id_rec,
                    "ID hit": id_hit,
                    "Def hit": def_hit,
                    "% Identidad": identidad,
                    "% Cobertura": cobertura,
                    "eValue": hsp.expect,
                    "Bit Score": hsp.bits})

                #Se almacenan las secuencias hit.
                hit_key = f"HIT_{id_hit}"
                if hit_key not in sec_fasta:
                    sec_fasta[hit_key] = hsp.sbjct.replace("-","")

    df = pd.DataFrame(filas)

    #Se genera el fichero FASTA final.
    fasta_out = ""

    for id_sec, sec in sec_fasta.items():
        fasta_out += f">{id_sec}\n{sec}\n"

    return df, fasta_out
                    

    
