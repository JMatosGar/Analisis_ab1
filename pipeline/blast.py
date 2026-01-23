import time
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

# Cache global
taxonomy_cache = {}

def taxid_from_accession(accession):
    try:
        handle = Entrez.esummary(db="nucleotide", id=accession, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.34)
        return record[0].get("TaxId")
    except Exception:
        return None

def taxonomy_from_taxid(taxid):
    lineage = {
        "Domain": None,
        "Phylum": None,
        "Class": None,
        "Order": None,
        "Family": None,
        "Genus": None,
        "Species": None,
    }

    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.4)
    except Exception:
        return lineage

    if not record:
        return lineage

    tax_rec = record[0]

    for node in tax_rec.get("LineageEx", []):
        rank = node["Rank"].capitalize()
        if rank in lineage:
            lineage[rank] = node["ScientificName"]

    if tax_rec.get("Rank") == "species":
        lineage["Species"] = tax_rec["ScientificName"]

    return lineage

def taxonomia(fasta, email):
    Entrez.email = email
    taxonomy_cache.clear()
    rows = []

    for record in SeqIO.parse(fasta, "fasta"):
        try:
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=record.seq,
                hitlist_size=1,
                expect=1e-10
            )

            blast_record = NCBIXML.read(result_handle)

            if not blast_record.alignments:
                rows.append({"ID": record.id})
                continue

            aln = blast_record.alignments[0]
            hsp = aln.hsps[0]
            accession = aln.accession.split(".")[0]

            if accession not in taxonomy_cache:
                taxid = taxid_from_accession(accession)
                taxonomy_cache[accession] = (
                    taxonomy_from_taxid(taxid) if taxid else
                    {k: None for k in ["Domain","Phylum","Class","Order","Family","Genus","Species"]}
                )

            tax = taxonomy_cache[accession]

            rows.append({
                "ID": record.id,
                **tax,
                "Hit accession": aln.accession,
                "Hit description": aln.hit_def,
                "% Identidad": round(hsp.identities / hsp.align_length * 100, 2),
                "e-value": hsp.expect,
                "Bit score": hsp.bits,
                "Cobertura (%)": round(hsp.align_length / len(record.seq) * 100, 2),
            })

        except Exception as e:
            rows.append({
                "ID": record.id,
                "Error": str(e)
            })

    return pd.DataFrame(rows)

