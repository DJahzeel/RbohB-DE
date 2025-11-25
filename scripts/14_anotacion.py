"""
kegg_annotation.py

Objetivo:
---------
- Construir la relación entre genes (locus_tag) y KO IDs de KEGG.
- Anotar las tablas de resultados de DESeq2 con KO IDs.
- Generar listas de KO IDs de genes diferencialmente expresados (UP/DOWN)
  por contraste, listas que servirán para análisis posteriores (enriquecimiento KEGG).

Supone que:
- Ya corriste el pipeline de expresión diferencial y tienes archivos
  *_deseq2_results.tsv en el directorio de DE.
- Tienes un GFF de referencia para PHAVU (con protein_id, locus_tag, etc.).
- Tienes un archivo PHAVU_to_KO.tsv (query -> KO). En el caso de PHAVU, se realizo una anotacion previa con eggnog-mapper.
    eggNOG-mapper se ocupo para hacer una correspondencia mediante ortologia de genes a ids KO de KEGG.
    El archivo PHAVU_to_KO.tsv, se proporciona en el repositorio para terminos practicos.

Uso típico (como script):
-------------------------
python kegg_annotation.py \
    --bioproject PRJNA656060 \
    --gff data/PHAVU/GCF_000499845.1_PhaVul1.0_genomic.gff \
    --kegg-orthologs data/PHAVU/PHAVU_to_KO.tsv \
    --padj-thr 0.01 \
    --lfc-thr 2.0
"""

import argparse
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

from paths.pathsval import diffexp_directory, project_root


# ============================================================
# 1. Funciones para leer GFF y construir mapas ID -> locus_tag
# ============================================================

def parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    """
    Parsea la columna de atributos de un GFF3 y regresa un diccionario
    clave=valor (ID, Parent, protein_id, locus_tag, etc.).
    """
    attrs: Dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def build_protein_to_locus_map(gff_path: Path) -> Dict[str, str]:
    """
    A partir de un archivo GFF3, construye un diccionario:

        protein_id -> locus_tag

    La lógica sigue el esquema del notebook:
    - Primer pase: mapear ID (gene/mRNA) -> locus_tag.
    - Segundo pase: para cada CDS, mapear protein_id -> locus_tag
      usando locus_tag directo o el de su Parent.
    """
    if not gff_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo GFF: {gff_path}")

    # 1er pase: ID (gene/mRNA) -> locus_tag
    id2locus: Dict[str, str] = {}

    with gff_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            ftype = cols[2]
            attrs = parse_gff_attributes(cols[8])

            if ftype in ("gene", "mRNA"):
                _id = attrs.get("ID")
                locus = attrs.get("locus_tag")
                if _id and locus:
                    id2locus[_id] = locus

    # 2º pase: protein_id -> locus_tag
    prot2locus: Dict[str, str] = {}

    with gff_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            ftype = cols[2]
            attrs = parse_gff_attributes(cols[8])

            if ftype == "CDS":
                prot = attrs.get("protein_id")
                parent = attrs.get("Parent")
                locus = attrs.get("locus_tag") or (id2locus.get(parent) if parent else None)

                if prot and locus:
                    prot2locus[prot] = locus

    print(f"[INFO] proteínas mapeadas protein_id -> locus_tag: {len(prot2locus)}")
    return prot2locus


# ============================================================
# 2. De PHAVU_to_KO.tsv a tabla gen (locus_tag) -> KO IDs
# ============================================================

def build_gene_to_ko_table(
    kegg_orthologs_path: Path,
    prot2locus: Dict[str, str],
) -> pd.DataFrame:
    """
    Lee el archivo de ortólogos de KEGG (p.ej. PHAVU_to_KO.tsv) y construye
    una tabla por gen (locus_tag) con sus KO IDs.

    El archivo de ortólogos se espera con formato:
        query<TAB>ko:KXXXXX
    sin encabezado (header=None) y puede tener líneas comentadas con '#'.
    """
    if not kegg_orthologs_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo de ortólogos KEGG: {kegg_orthologs_path}")

    ko_df = pd.read_csv(
        kegg_orthologs_path,
        sep="\t",
        comment="#",
        names=["query", "ko_id"],
        header=None,
    )

    # Mapear query (protein_id) -> locus_tag usando el diccionario del GFF
    ko_df["locus_tag"] = ko_df["query"].map(prot2locus)
    ko_df = ko_df.dropna(subset=["locus_tag"])

    # Limpiar prefijo "ko:" si está presente
    ko_df["ko_id"] = ko_df["ko_id"].str.replace("ko:", "", regex=False)

    # Si hay varios KO por proteína separados por coma, separarlos
    ko_df = ko_df.assign(ko_id=ko_df["ko_id"].str.split(",")).explode("ko_id")

    # Colapsar a 1 fila por gen con KOs únicos (unidos por coma)
    ko_gene = (
        ko_df.groupby("locus_tag")["ko_id"]
        .apply(lambda x: ",".join(sorted(set([k for k in x if k and k != "-"]))))
        .reset_index()
    )

    print(f"[INFO] genes con al menos un KO anotado: {ko_gene.shape[0]}")
    return ko_gene


# ============================================================
# 3. Anotar resultados DESeq2 con KO IDs
# ============================================================

def annotate_deseq_results_with_ko(
    diffexp_dir: Path,
    ko_gene: pd.DataFrame,
    pattern: str = "*_deseq2_results.tsv",
    gene_id_col: str = "gene_id",
) -> Dict[Path, Path]:
    """
    Para cada archivo *_deseq2_results.tsv en diffexp_dir:
    - Deriva locus_tag a partir de la columna gene_id (quitando 'ID=' delante).
    - Hace merge con la tabla ko_gene (locus_tag -> ko_id).
    - Guarda un archivo .with_KO.tsv al lado del original.

    Devuelve un diccionario: {ruta_original -> ruta_con_KO}.
    """
    mapping: Dict[Path, Path] = {}

    for tsv in sorted(diffexp_dir.glob(pattern)):
        print(f"[INFO] Anotando KO en: {tsv.name}")
        df = pd.read_csv(tsv, sep="\t")

        if gene_id_col not in df.columns:
            raise ValueError(
                f"El archivo {tsv} no tiene la columna '{gene_id_col}'. "
                "Asegúrate de que los resultados de DESeq2 incluyan esa columna."
            )

        # locus_tag típico: quitar el prefijo 'ID=' de gene_id si existe
        df["locus_tag"] = df[gene_id_col].astype(str).str.replace("^ID=", "", regex=True)

        # Merge con la tabla de KOs
        res_kegg = df.merge(ko_gene, on="locus_tag", how="left")  # agrega columna 'ko_id'

        out_path = tsv.with_suffix(".with_KO.tsv")
        res_kegg.to_csv(out_path, sep="\t", index=False)
        mapping[tsv] = out_path

        print(f"[INFO] Guardado: {out_path.name} (n={res_kegg.shape[0]} filas)")

    if not mapping:
        print(f"[ADVERTENCIA] No se encontraron archivos que cumplan el patrón {pattern} en {diffexp_dir}")

    return mapping


# ============================================================
# 4. Clasificar DEGs y generar listas de KO IDs (UP/DOWN)
# ============================================================

def load_and_classify_kos(
    path: Path,
    padj_cut: float,
    lfc_cut: float,
) -> Tuple[pd.Series, pd.Series]:
    """
    Dado un archivo .with_KO.tsv, genera dos listas de KO IDs:

      - up_kos: KOs de genes UP (padj <= padj_cut, log2FC >= lfc_cut)
      - down_kos: KOs de genes DOWN (padj <= padj_cut, log2FC <= -lfc_cut)

    El archivo debe tener al menos:
      - columna 'ko_id' (puede contener varios KO separados por coma)
      - columnas 'padj' y 'log2FoldChange'
    """
    enrichment = pd.read_csv(path, sep="\t")

    # Nos quedamos únicamente con filas que tengan KOs
    enrichment = enrichment[
        enrichment["ko_id"].notna() & (enrichment["ko_id"].astype(str) != "")
    ].copy()

    # Expandir KOs separados por coma
    enrichment_long = enrichment.assign(
        ko_id=enrichment["ko_id"].astype(str).str.split(",")
    ).explode("ko_id")

    # Asegurar que padj y log2FoldChange existan y no sean NaN
    enrichment_long = enrichment_long.dropna(subset=["padj", "log2FoldChange"])

    up = enrichment_long[
        (enrichment_long["padj"] <= padj_cut)
        & (enrichment_long["log2FoldChange"] >= lfc_cut)
    ]
    down = enrichment_long[
        (enrichment_long["padj"] <= padj_cut)
        & (enrichment_long["log2FoldChange"] <= -lfc_cut)
    ]

    return up["ko_id"], down["ko_id"]


def generate_deg_ko_lists(
    diffexp_dir: Path,
    padj_cut: float,
    lfc_cut: float,
    suffix: str = ".with_KO.tsv",
) -> None:
    """
    Recorre todos los archivos *<suffix> en diffexp_dir (p.ej. *_deseq2_results.with_KO.tsv),
    genera listas de KO IDs para genes UP y DOWN, y las guarda como:

      <nombre>_UP_genes.txt
      <nombre>_DOWN_genes.txt

    (Siguiendo la convención del notebook original, aunque sean KO IDs).
    """
    for path in sorted(diffexp_dir.glob(f"*{suffix}")):
        name = path.stem  # p.ej. 'Nod_vs_Control_deseq2_results.with_KO'

        print(f"[INFO] Generando listas de KO para: {name}")
        up_kos, down_kos = load_and_classify_kos(path, padj_cut, lfc_cut)

        # Diccionario, por si luego quieres usarlo en código en lugar de sólo archivos
        print(f"  - UP:   {len(up_kos)} KOs")
        print(f"  - DOWN: {len(down_kos)} KOs")

        up_out = diffexp_dir / f"{name}_UP_genes.txt"
        down_out = diffexp_dir / f"{name}_DOWN_genes.txt"

        pd.Series(up_kos).drop_duplicates().to_csv(
            up_out,
            index=False,
            header=True,
        )
        pd.Series(down_kos).drop_duplicates().to_csv(
            down_out,
            index=False,
            header=True,
        )

        print(f"  Guardados:\n    {up_out.name}\n    {down_out.name}")


# ============================================================
# 5. main() para uso como script
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Anotación de resultados de DESeq2 con KO IDs y generación de listas "
            "de KO para genes diferencialmente expresados por contraste."
        )
    )

    parser.add_argument(
        "--bioproject",
        default="PRJNA656060",
        help="ID del BioProject (por defecto: PRJNA656060).",
    )
    parser.add_argument(
        "--gff",
        required=True,
        help="Ruta al archivo GFF3 del genoma (con protein_id, locus_tag, etc.).",
    )
    parser.add_argument(
        "--kegg-orthologs",
        required=True,
        help="Ruta al archivo PHAVU_to_KO.tsv (query -> KO IDs).",
    )
    parser.add_argument(
        "--padj-thr",
        type=float,
        default=0.01,
        help="Umbral de padj para genes DE (default: 0.01).",
    )
    parser.add_argument(
        "--lfc-thr",
        type=float,
        default=2.0,
        help="Umbral de |log2FC| para genes DE (default: 2.0).",
    )

    args = parser.parse_args()

    root = project_root()
    diffexp_dir = diffexp_directory(root) / args.bioproject
    diffexp_dir.mkdir(parents=True, exist_ok=True)

    gff_path = Path(args.gff)
    kegg_orthologs_path = Path(args.kegg_orthologs)

    # 1) Construir mapa protein_id -> locus_tag
    prot2locus = build_protein_to_locus_map(gff_path)

    # 2) Construir tabla locus_tag -> KO IDs
    ko_gene = build_gene_to_ko_table(kegg_orthologs_path, prot2locus)

    # (Opcional) Guardar la tabla gene->KO para referencia
    gene_ko_out = diffexp_dir / "gene_to_KO.tsv"
    ko_gene.to_csv(gene_ko_out, sep="\t", index=False)
    print(f"[INFO] Tabla gene_to_KO guardada en: {gene_ko_out}")

    # 3) Anotar todos los resultados de DESeq2 en diffexp_dir con KO IDs
    annotate_deseq_results_with_ko(diffexp_dir, ko_gene)

    # 4) Generar listas de KO IDs para genes DE por contraste
    generate_deg_ko_lists(
        diffexp_dir,
        padj_cut=args.padj_thr,
        lfc_cut=args.lfc_thr,
        suffix=".with_KO.tsv",
    )

    print(f"[OK] Anotación KEGG completa en: {diffexp_dir}")

if __name__ == "__main__":
    main()
