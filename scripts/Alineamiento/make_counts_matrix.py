#!/usr/bin/env python3
import sys
import pathlib

# Detectar raíz del proyecto (script en scripts/Alineamiento)
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

COUNTS_DIR = PROJECT_ROOT / "results" / "counts"
OUT_COUNTS = COUNTS_DIR / "All_counts_rbohb.txt"
SAMPLES_TSV = PROJECT_ROOT / "samples.tsv"

print("[INFO] make_counts_matrix.py")
print(f"  - PROJECT_ROOT: {PROJECT_ROOT}")
print(f"  - COUNTS_DIR:   {COUNTS_DIR}")
print(f"  - OUT_COUNTS:   {OUT_COUNTS}")
print(f"  - SAMPLES_TSV:  {SAMPLES_TSV}")

#  Validar directorio de counts
if not COUNTS_DIR.is_dir():
    print(f"ERROR: No existe el directorio de counts: {COUNTS_DIR}", file=sys.stderr)
    print("       Corre primero el script de counts (coverageBed).", file=sys.stderr)
    sys.exit(1)

#  Buscar archivos *.counts.txt
count_files = sorted(COUNTS_DIR.glob("*.counts.txt"))
if not count_files:
    print(f"ERROR: No se encontraron *.counts.txt en {COUNTS_DIR}", file=sys.stderr)
    print("       Corre primero el script de counts (coverageBed).", file=sys.stderr)
    sys.exit(1)

print("[INFO] Archivos de entrada:")
for p in count_files:
    print("  -", p.name)


def get_gene_id(attr_str: str):
    """
    Atributos GFF (columna 9) -> un ID tipo:
      ID=PHAVU_001G000400   (prioridad locus_tag > Name > ID)
    """
    kv = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        kv[k.strip()] = v.strip()

    if "locus_tag" in kv:
        return "ID=" + kv["locus_tag"]
    if "Name" in kv:
        return "ID=" + kv["Name"]
    if "ID" in kv:
        return kv["ID"]
    return None


# counts[gene_id][sample] = count
counts: dict[str, dict[str, int]] = {}
sample_names: list[str] = []

for path in count_files:
    sample = path.name.replace(".counts.txt", "")
    sample_names.append(sample)

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

    # filas
    for gene in genes:
        row = [gene]
        for sample in sample_names:
            row.append(str(counts[gene].get(sample, 0)))
        out.write("\t".join(row) + "\n")

print("[OK] Matriz de conteos generada:")
print("     ", OUT_COUNTS)

