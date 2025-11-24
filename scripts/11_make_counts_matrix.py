#!/usr/bin/env python3
# -------------------------------------------------------------------------
# Script: 11_make_counts_matrix.py
# Descripción:
#   Lee todos los archivos de conteos generados por coverageBed:
#
#       results/<BIOPROJECT_ID>/counts/*.counts.txt
#
#   y construye una matriz de conteos con formato tabular:
#
#       gene_id <TAB> sample1 <TAB> sample2 <TAB> ...
#       geneA   <TAB>   c_1   <TAB>   c_2   <TAB> ...
#
#   Donde:
#     - Cada archivo <sample>.counts.txt corresponde a una muestra.
#     - Cada fila en los .counts.txt proviene de coverageBed con un GFF.
#     - Se usa la columna 9 (atributos GFF) para extraer un identificador
#       de gen estable (locus_tag / Name / ID).
#     - Se usa la columna 10 (primer número de coverageBed) como conteo.
#
#   El archivo de salida se escribe en:
#
#       results/<BIOPROJECT_ID>/counts/All_counts_rbohb.txt
#
# Uso:
#   python3 11_make_counts_matrix.py
#   python3 11_make_counts_matrix.py --bioproject 482464
#
# Flags:
#   --bioproject/-b  ID de subcarpeta bajo results/ (por defecto "482464")
#
# Dependencias:
#   - Haber corrido antes:
#       * el alineamiento (para generar BAMs),
#       * el script de counts (coverageBed) para producir *.counts.txt.
# -------------------------------------------------------------------------
import sys
import pathlib

# Detectar raíz del proyecto (script en scripts/Alineamiento)
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

# Parámetro por defecto, se puede cambiar por línea de comandos
BIOPROJECT_ID = "482464"

# Paso de parametros por línea de comandos
args = sys.argv[1:]
if "--help" in args or "-h" in args:
    print("Uso: python3 11_make_counts_matrix.py [--bioproject ID]")
    print()
    print("  --bioproject/-b  ID de subcarpeta bajo results/ (default: 482464)")
    print()
    sys.exit(0)

i = 0
while i < len(args):
    arg = args[i]
    if arg in ("--bioproject", "-b"):
        i += 1
        BIOPROJECT_ID = args[i]
    else:
        print(f"ERROR: Argumento desconocido: {arg}", file=sys.stderr)
        sys.exit(1)
    i += 1

# Rutas de entrada/salida
COUNTS_DIR = PROJECT_ROOT / "results" / BIOPROJECT_ID / "counts"
OUT_COUNTS = COUNTS_DIR / "All_counts_rbohb.txt"
SAMPLES_TSV = PROJECT_ROOT / "samples.tsv"

# Información inicial para el usuario
print("[INFO] make_counts_matrix.py")
print(f"  - PROJECT_ROOT: {PROJECT_ROOT}")
print(f"  - COUNTS_DIR:   {COUNTS_DIR}")
print(f"  - OUT_COUNTS:   {OUT_COUNTS}")
print(f"  - SAMPLES_TSV:  {SAMPLES_TSV}")

# Validar directorio de counts
if not COUNTS_DIR.is_dir():
    print(f"ERROR: No existe el directorio de counts: {COUNTS_DIR}", file=sys.stderr)
    print("       Corre primero el script de counts (coverageBed).", file=sys.stderr)
    sys.exit(1)

# Buscar archivos *.counts.txt
count_files = sorted(COUNTS_DIR.glob("*.counts.txt"))
if not count_files:
    print(f"ERROR: No se encontraron *.counts.txt en {COUNTS_DIR}", file=sys.stderr)
    print("       Corre primero el script de counts (coverageBed).", file=sys.stderr)
    sys.exit(1)
# Mostrar archivos encontrados
print("[INFO] Archivos de entrada:")
for p in count_files:
    print("  -", p.name)

# Función para extraer gene_id de atributos GFF
def get_gene_id(attr_str: str):
    """
    Atributos GFF (columna 9) -> un ID tipo:
      prioridad: locus_tag > Name > ID
    Devuelve solo el valor (sin 'ID=' al inicio).
    """
    kv = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        kv[k.strip()] = v.strip()

    if "locus_tag" in kv:
        return kv["locus_tag"]
    if "Name" in kv:
        return kv["Name"]
    if "ID" in kv:
        return kv["ID"]
    return None


# Estructura para almacenar conteos:
#  counts[gene_id][sample] = count
counts: dict[str, dict[str, int]] = {}
# Lista de nombres de muestras
sample_names: list[str] = []

for path in count_files:
    # # Nombre de muestra = nombre del archivo sin el sufijo .counts.txt
    sample = path.name.replace(".counts.txt", "")
    sample_names.append(sample)
    # Conteos específicos de esta muestra
    sample_counts: dict[str, int] = {}

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 10:
                continue

            feature = cols[2]
            if feature != "gene":
                # sólo usamos líneas donde la columna 3 es "gene"
                continue

            attr = cols[8]
            gene_id = get_gene_id(attr)
            if gene_id is None:
                continue

            # columna 10 = primer número de coverageBed (nº lecturas)
            c_str = cols[9]
            try:
                c = int(float(c_str))
            except ValueError:
                continue
            # Acumulamos conteos por gene_id dentro de la muestra
            sample_counts[gene_id] = sample_counts.get(gene_id, 0) + c

    # mezclar en la estructura global
    for gene_id, c in sample_counts.items():
        if gene_id not in counts:
            counts[gene_id] = {}
        counts[gene_id][sample] = counts[gene_id].get(sample, 0) + c

# Ordenar genes y muestras para que quede estable
genes = sorted(counts.keys())
sample_names = sorted(set(sample_names))

print(f"[INFO] Nº de genes en la matriz: {len(genes)}")
print(f"[INFO] Nº de muestras en counts: {len(sample_names)}")

# Escribir matriz
with OUT_COUNTS.open("w") as out:
    # cabecera
    out.write("gene_id\t" + "\t".join(sample_names) + "\n")

    # filas : un renglón por gen, columnas = conteos por muestra
    for gene in genes:
        row = [gene]
        for sample in sample_names:
            # Si algún gene no está presente en una muestra → 0
            row.append(str(counts[gene].get(sample, 0)))
        out.write("\t".join(row) + "\n")

print("[OK] Matriz de conteos generada:")
print("     ", OUT_COUNTS)