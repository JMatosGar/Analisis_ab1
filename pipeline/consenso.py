import statistics
import os
from Bio import pairwise2
from Bio.Seq import Seq
import pandas as pd

#Se definen los primers conocidos y sus direcciones.

def primers(filename):
    fname = os.path.basename(filename).replace(".ab1", "")
    parts = fname.split("_")

    primer = parts[1].upper() if len(parts) > 1 else ""

    if fname.lower().endswith("(reversed)"):
        return "reverse", primer

    for_list = ["ITS1", "EUKA", "AB1", "TETRAS PM", "RHULVA", "NRSSU178",
                    "CFD", "LSUB", "U1", "FPU", "DAITOMC", "ITSCIANO1"]
    rev_list = ["ITS4", "EUKB", "AB2", "TETRAS TH", "B", "DPRBCL7",
                    "U3", "DAITOMD", "ITSCIANO2", "ITSCIANO3", "DP"]

    if primer.endswith("F") or primer in for_list:
        return "forward", primer
    if primer.endswith("R") or primer in rev_list:
        return "reverse", primer

    return None, None

#Se crea una función que determina la secuencia consenso.
def consenso(forward, reverse, forward_cal, reverse_cal,
    umbral=20, min_overlap=50, delta=5):
    
    reverse_rc = str(Seq(reverse).reverse_complement())
    reverse_qual = reverse_qual[::-1]

    alignments = pairwise2.align.localms(
        forward, reverse_rc,
        2, -1, -5, -2)

    if not alignments:
        return "", 0, True, 0

    aln_f, aln_r = alignments[0][0], alignments[0][1]

    consensus = []
    consensus_cal = []

    f_i = r_i = 0

    for a_f, a_r in zip(aln_f, aln_r):

        f_q = forward_cal[f_i] if a_f != "-" else None
        r_q = reverse_cal[r_i] if a_r != "-" else None

        if a_f != "-" and a_r != "-" and a_f == a_r:
            if f_q + r_q >= 2 * umbral:
                consensus.append(a_f)
                consensus_cal.append(max(f_q, r_q))
            else:
                consensus.append("N")
                consensus_cal.append(0)

        elif a_f != "-" and a_r != "-" and a_f != a_r:
            if f_q >= r_q + delta and f_q >= umbral:
                consensus.append(a_f)
                consensus_cal.append(f_q)
            elif r_q >= f_q + delta and r_q >= umbral:
                consensus.append(a_r)
                consensus_cal.append(r_q)
            else:
                consensus.append("N")
                consensus_cal.append(0)

        elif a_f != "-" and a_r == "-":
            consensus.append(a_f if f_q >= umbral else "N")
            consensus_cal.append(f_q if f_q >= umbral else 0)

        elif a_r != "-" and a_f == "-":
            consensus.append(a_r if r_q >= umbral else "N")
            consensus_cal.append(r_q if r_q >= umbral else 0)

        if a_f != "-":
            f_i += 1
        if a_r != "-":
            r_i += 1

    consensus_seq_raw = "".join(consensus)

    leading_n = len(consensus_seq_raw) - len(consensus_seq_raw.lstrip("N"))
    trailing_n = len(consensus_seq_raw) - len(consensus_seq_raw.rstrip("N"))

    consensus_seq = consensus_seq_raw.strip("N")
    consensus_cal = consensus_cal[leading_n: len(consensus_cal) - trailing_n]

    bad_overlap = len(consensus_seq) < min_overlap

    valid_cals = [q for q in consensus_cal if q > 0]
    mean_q = statistics.mean(valid_cals) if valid_cals else 0

    return consensus_seq, mean_q, bad_overlap, leading_n + trailing_n


#Se crea una función que integra estas funciones y limita los resultados generados.
def generar_consensos(trimmed_df, umbral=20):
    df = trimmed_df.copy()

    df["Orientation"] = None
    df["Primer"] = None

    for i, row in df.iterrows():
        orient, primer = primers(row["ID"])
        df.at[i, "Orientation"] = orient
        df.at[i, "Primer"] = primer

    df = df[df["Orientation"].notna()]

    df["Sample"] = df["ID"].str.split("_").str[0]

    results = []

    for sample, group in df.groupby("Sample"):
        fwd = group[group["Orientation"] == "forward"]
        rev = group[group["Orientation"] == "reverse"]

        if fwd.empty or rev.empty:
            continue

        for _, f in fwd.iterrows():
            for _, r in rev.iterrows():

                cons_seq, mean_q, bad_ov, n_trim = construir_consenso(
                    f["Secuencia cortada"],
                    r["Secuencia cortada"],
                    f["Puntuación cortada"],
                    r["Puntuación cortada"],
                    umbral=umbral,
                    min_overlap=min_overlap
                )

                pct_n = cons_seq.count("N") / max(len(cons_seq), 1) * 100

                results.append({
                    "Sample": sample,
                    "Forward primer": f["Primer"],
                    "Reverse primer": r["Primer"],
                    "Consensus sequence": cons_seq,
                    "Consensus length": len(cons_seq),
                    "Mean quality": mean_q,
                    "Percent N": pct_n,
                    "Bad overlap": bad_ov,
                    "Trimmed Ns": n_trim})

    cons_df = pd.DataFrame(results)
    cons_df = cons_df[(cons_df["Consensus length"] >= 150) & (cons_df["Percent N"] <= 5) & (cons_df["Mean quality"] >= umbral)]

    return cons_df
