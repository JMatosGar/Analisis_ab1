import tempfile
from collections import defaultdict
from io import StringIO
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline, MuscleCommandline

#Se crea una función que realiza el agrupado de las secuencias por ID de muestra.
def _agrupar_por_muestra(fasta_str):
    io_fasta = StringIO(fasta_str)
    records = list(SeqIO.parse(io_fasta, "fasta"))

    grupos = defaultdict(list)
    for rec in records:
        muestra = rec.id.split("_Hit_")[0]
        grupos[muestra].append(rec)

    return grupos


def mafft(fasta, modo="Auto"):
    #Se aplica el agrupado de las muestras.  
    grupos = _agrupar_por_muestra(fasta)
    final = ""

    #Se itera por cada muestra y sus respectivos hits.
    for muestra, recs in grupos.items():
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tmp_in:
            SeqIO.write(recs, tmp_in, "fasta")
            tmp_in.flush()

            #Se configura MAFFT en base a su modo.
            if modo == "Auto":
                cmd = MafftCommandline(input=tmp_in.name)
            elif modo == "FFT-NS-2":
                cmd = MafftCommandline(input=tmp_in.name, fftns=True)
            elif modo == "L-INS-i":
                cmd = MafftCommandline(input=tmp_in.name, localpair=True, maxiterate=1000)

            stdout, stderr = cmd()
            final += stdout + "\n"

    return final

#Se repite el proceso con MUSCLE.
def muscle(fasta):
    grupos = _agrupar_por_muestra(fasta)
    final = ""

    for muestra, recs in grupos.items():
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tmp_in:
            SeqIO.write(recs, tmp_in, "fasta")
            tmp_in.flush()

            cmd = MuscleCommandline(input=tmp_in.name)
            stdout, stderr = cmd()
            final += stdout + "\n"

    return final
