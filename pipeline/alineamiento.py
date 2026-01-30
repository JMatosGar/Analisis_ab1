import subprocess
import tempfile
from collections import defaultdict
from io import StringIO
from Bio import SeqIO
import shutil

def agrupar(fasta_str):
    records = list(SeqIO.parse(StringIO(fasta_str), "fasta"))
    grupos = defaultdict(list)

    for rec in records:
        muestra = rec.id.split("_Hit_")[0]
        grupos[muestra].append(rec)

    return grupos

def mafft(fasta, modo="Auto"):
    if not shutil.which("mafft"):
        raise RuntimeError("❌ MAFFT no está disponible en el PATH del sistema")

    grupos = agrupar(fasta)
    final = ""

    for muestra, recs in grupos.items():
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tmp_in:
            SeqIO.write(recs, tmp_in, "fasta")
            tmp_in.flush()

            cmd = ["mafft"]

            if modo == "FFT-NS-2":
                cmd.append("--fftns")
            elif modo == "L-INS-i":
                cmd.extend(["--localpair", "--maxiterate", "1000"])

            cmd.append(tmp_in.name)

            result = subprocess.run(cmd, capture_output=True, text=True, check=True, shell=True)

            final += result.stdout + "\n"

    return final


def muscle(fasta):
    if not shutil.which("muscle"):
        raise RuntimeError("❌ MUSCLE no está disponible en el PATH del sistema")

    grupos = agrupar(fasta)
    final = ""

    for muestra, recs in grupos.items():
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tmp_in:
            SeqIO.write(recs, tmp_in, "fasta")
            tmp_in.flush()

            cmd = ["muscle", "-align", tmp_in.name]

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            final += result.stdout + "\n"

    return final
