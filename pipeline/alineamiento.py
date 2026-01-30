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
    grupos = agrupar(fasta)
    final = ""

    for muestra, recs in grupos.items():
        if not recs:
            continue

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tmp_in:
            SeqIO.write(recs, tmp_in, "fasta")
            tmp_in.flush()

            #Ruta completa a MAFFT en Windows
            cmd_mafft = r"C:\MAFFT\mafft.bat"

            #Construir comando según modo
            cmd = [cmd_mafft]
            if modo == "FFT-NS-2":
                cmd.append("--fftns")
            elif modo == "L-INS-i":
                cmd.extend(["--localpair", "--maxiterate", "1000"])

            cmd.append(tmp_in.name)

            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True, shell=True)
            
            except subprocess.CalledProcessError as e:
                print("STDOUT:", e.stdout)
                print("STDERR:", e.stderr)
                raise RuntimeError(f"❌ Error ejecutando MAFFT para {muestra}")

            final += result.stdout + "\n"

    return final



def muscle(fasta):
    grupos = agrupar(fasta)
    final = ""

    cmd_muscle = r"C:\MUSCLE\muscle.exe" 

    for muestra, recs in grupos.items():
        if not recs:
            continue

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tmp_in:
            SeqIO.write(recs, tmp_in, "fasta")
            tmp_in.flush()

            cmd = [cmd_muscle, "-align", tmp_in.name]

            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True, shell=True)
            
            except subprocess.CalledProcessError as e:
                print("STDOUT:", e.stdout)
                print("STDERR:", e.stderr)
                raise RuntimeError(f"❌ Error ejecutando MUSCLE para {muestra}")

            final += result.stdout + "\n"

    return final

