import subprocess
import tempfile
import os
import shutil
from io import BytesIO, StringIO
from Bio import Phylo
import matplotlib.pyplot as plt

#Función para obtener las características del arbol generado
def resumen_arbol(tree):
  resumen = {}

  for linea in tree.splitlines():
    linea = linea.strip()

    if linea.startswith("Best-fit model according to"):
      resumen["Modelo evolutivo"] = linea.split(":")[-1].strip()

    elif linea.startswith("Log-likelihood of the tree"):
      resumen["Log-similitud"] = float(linea.split(":")[-1])

    elif linea.startswith("Total tree length"):
      resumen["Longitud del árbol"] = float(linea.split(":")[-1])

    elif linea.startswith("Number of sequences"):
      resumen["Número de secuencias"] = int(linea.split(":")[-1])

    elif linea.startswith("Alignment length"):
      resumen["Longitud del alineamiento"] = int(linea.split(":")[-1])

    elif "Bootstrap" in linea and "replicates" in linea:
      resumen["Bootstrap"] = linea.split(":")[-1].strip()

  return resumen

#Función de IQTree desde el PC.
def iqtree(fasta_alin, bootstrap = 1000, modelo = "MFP"):
  #Se comprueba que está IQTree en el PATH del sistema.
  if not shutil.which("iqtree3"):
    raise RuntimeError("❌ IQ-Tree no se encuentra en el PATH")

  with tempfile.TemporaryDirectory() as tmpdir:
    fasta_path = os.path.join(tmpdir, "alineamiento.fasta")
    with open(fasta_path, "w") as f:
      f.write(fasta_alin)

    #Se establecen las condiciones de IQ-TREE y se ejecuta.
    cmd = ["iqtree3",
           "-s", fasta_path,
           "-m", modelo,
           "-B", str(bootstrap),
           "-nt", "AUTO",
           "--quiet"]

    resultado = subprocess.run(cmd, cwd=tmpdir, capture_output=True, text=True)

    if resultado.returncode !=0:
      raise RuntimeError(f"❌ Error al ejecutar IQ-TREE:\n{resultado.stderr}")

    #Se almacenan los resultados.
    treefile = fasta_path + ".treefile"
    iqtree_file = fasta_path + ".iqtree"

    if not os.path.exists(treefile) or not os.path.exists(iqtree_file):
      raise RuntimeError("❌ IQ-TREE no ha podido generar los ficheros esperados")

    with open(treefile) as f:
      arbol_newick = f.read().strip()

    with open(iqtree_file) as f:
      iqtree = f.read()

    resumen = resumen_arbol(iqtree)

  return arbol_newick, resumen, iqtree

#Se incluye una función para representar el árbol generado.
def mostrar_arbol(newick, distancia = False):
  #Se lee el árbol y se toman sus características.
  tree = Phylo.read(StringIO(newick), "newick")

  fig = plt.figure(figsize = (10, 10))
  ax = fig.add_subplot(1, 1, 1)

  if distancia:
    Phylo.draw(tree, axes = ax, do_show = False, show_confidence = True,
                branch_labels=lambda c: f"{c.branch_length:.3f}" if c.branch_length else "")
  else:
    Phylo.draw(tree, axes = ax, do_show = False, show_confidence = True)

  buffer = BytesIO()
  plt.savefig(buf, format="png", bbox_inches="tight", dpi=300)
  plt.close(fig)
  buf.seek(0)

  return buffer
                                     
