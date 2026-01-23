import os
import subprocess
import time
import pandas as pd
from io import StringIO
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML

# Cache de taxonomía
taxonomy_cache = {}

def get_taxid_from_accession(accession):
    try:
        handle = Entrez.esummary(
            db="nucleotide",
            id=accession,
            retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.34)
        return record[0].get("TaxId")
    except Exception:
        return None


def get_taxonomy_from_taxid(taxid):
    if taxid in taxonomy_cache:
        return taxonomy_cache[taxid]

    lineage = {
        "Domain": None,
        "Phylum": None,
        "Class": None,
        "Order": None,
        "Family": None,
        "Genus": None,
        "Species": None
    }

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
        taxonomy_cache[taxid] = lineage
        return lineage

    tax_record = record[0]

    for node in tax_record.get("LineageEx", []):
        rank = node["Rank"].capitalize()
        if rank in lineage:
            lineage[rank] = node["ScientificName"]

    if tax_record.get("Rank") == "species":
        lineage["Species"] = tax_record["ScientificName"]

    taxonomy_cache[taxid] = lineage
    return lineage


def taxonomia(fasta_path, email, top_hits=1, db="nt", task="megablast"):
    Entrez.email = email
    taxonomy_cache.clear()

    rows = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        query_id = record.id

        # Crear FASTA temporal por secuencia
        with open("query_tmp.fasta", "w") as f:
            SeqIO.write(record, f, "fasta")

        # Ejecutar BLAST remoto (XML)
        cmd = [
            "blastn",
            "-query", "query_tmp.fasta",
            "-db", db,
            "-remote",
            "-task", task,
            "-max_target_seqs", str(top_hits),
            "-outfmt", "5"
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        if not result.stdout.strip():
            continue

        blast_record = NCBIXML.read(StringIO(result.stdout))

        if not blast_record.alignments:
            rows.append({
                "ID": query_id,
                "Domain": None,
                "Phylum": None,
                "Class": None,
                "Order": None,
                "Family": None,
                "Genus": None,
                "Species": None,
                "Hit": None,
                "Percent_identity": 0,
                "Coverage": 0
            })
            continue

        aln = blast_record.alignments[0]
        hsp = aln.hsps[0]

        accession = aln.accession
        taxid = get_taxid_from_accession(accession)

        taxonomy = (
            get_taxonomy_from_taxid(taxid)
            if taxid else
            {k: None for k in
             ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]}
        )

        rows.append({
            "ID": query_id,
            **taxonomy,
            "Hit": aln.hit_def,
            "Percent_identity": round(hsp.identities / hsp.align_length * 100, 2),
            "Coverage": round(hsp.align_length / len(record.seq) * 100, 2),
            "Evalue": hsp.expect,
            "Bit_score": hsp.bits
        })

        time.sleep(1)  # pausa BLAST

    return pd.DataFrame(rows)
