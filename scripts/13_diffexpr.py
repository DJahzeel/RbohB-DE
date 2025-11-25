"""
diffexpr.py

Script para correr análisis de expresión diferencial con PyDESeq2
y generar las gráficas principales (boxplot, densidad, PCA, volcanos),
a partir de una matriz de conteos y un archivo con los nombres de muestra. Pensado primordialmente
para el análisis del BioProject PRJNA482464, pero adaptable a otros datasets similares.

"""

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from paths.pathsval import diffexp_directory, project_root


# ============================================================
# Utilidades básicas
# ============================================================

def get_default_paths(
    bioproject: str,
    counts_path: str | Path | None = None,
    srr_ids_path: str | Path | None = None,
    outdir: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    """
    Construye rutas por defecto a partir de project_root() y diffexp_directory().
    """
    root = project_root()

    if counts_path is None:
        counts_path = root / "results" / bioproject / "All_counts_rbohb.txt"
    else:
        counts_path = Path(counts_path)

    if srr_ids_path is None:
        srr_ids_path = root / "data" / bioproject / "srr_ids.txt"
    else:
        srr_ids_path = Path(srr_ids_path)

    if outdir is None:
        outdir = diffexp_directory(root) / bioproject
    else:
        outdir = Path(outdir)

    outdir.mkdir(parents=True, exist_ok=True)
    return counts_path, srr_ids_path, outdir


def load_counts_and_rename(
    counts_path: Path,
    srr_ids_path: Path,
) -> pd.DataFrame:
    """
    Lee la matriz de conteos y renombra las columnas usando srr_ids.txt,
    manteniendo los nombres de muestra tal como están en SampleName.
    """
    if not counts_path.exists():
        raise FileNotFoundError(f"No se encontró la matriz de conteos: {counts_path}")

    if not srr_ids_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo srr_ids: {srr_ids_path}")

    df_counts = pd.read_csv(counts_path, sep="\t")
    if "gene_id" in df_counts.columns:
        df_counts = df_counts.set_index("gene_id")

    srr_ids = pd.read_csv(srr_ids_path, sep="\t")
    if not {"Run", "SampleName"}.issubset(srr_ids.columns):
        raise ValueError(
            f"El archivo {srr_ids_path} debe contener las columnas 'Run' y 'SampleName'."
        )

    sra_subset = srr_ids[["Run", "SampleName"]]
    samples_dict = dict(zip(sra_subset["Run"], sra_subset["SampleName"]))

    # Renombra columnas SRR -> nombres de muestra
    df_counts_renamed = df_counts.rename(columns=samples_dict)

    return df_counts_renamed


def get_condition_from_sample(sample_name: str) -> str:
    """
    Dado un nombre de muestra tipo 'Control_Nod_1',
    devuelve la condición sin el sufijo de réplica: 'Control_Nod'.
    """
    return re.sub(r"_\d+$", "", sample_name)


def build_conditions_series(df_counts_renamed: pd.DataFrame) -> pd.Series:
    """
    Construye un Series con la condición de cada muestra a partir de los nombres
    de columna de la matriz de conteos renombrada.
    """
    conds = [get_condition_from_sample(s) for s in df_counts_renamed.columns]
    return pd.Series(conds, index=df_counts_renamed.columns, name="condition")


# ============================================================
# Transformaciones y plots exploratorios
# ============================================================

def compute_log2_and_long_format(
    df_counts_renamed: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Devuelve:
      - df_log2: matriz de conteos log2(x+1)
      - df_melted: formato largo para boxplots/densidad
    """
    df_log2 = np.log2(df_counts_renamed + 1)
    df_melted = df_log2.melt(var_name="Sample", value_name="Value")
    df_melted["Cond"] = df_melted["Sample"].str.replace("_.$", "", regex=True)
    return df_log2, df_melted


def plot_boxplot_log2(df_melted: pd.DataFrame, outdir: Path) -> Path:
    """
    Boxplot de conteos log2 por muestra, coloreado por condición.
    """
    colors = ["#009ad6", "#99d6ee", "#00638a", "#e5f4fa", "#90aeeb", "#3267c0"]

    plt.figure(figsize=(10, 8))
    sns.boxplot(
        data=df_melted,
        y="Sample",
        x="Value",
        hue="Cond",
        palette=colors,
    )
    plt.yticks(fontsize=9)
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", title="Cond")
    plt.title("Boxplot de expresiones log₂ normalizadas")
    plt.xlabel("log₂(counts + 1)")
    plt.tight_layout()

    outpath = outdir / "01_boxplot_log2.png"
    plt.savefig(outpath, dpi=300)
    plt.close()
    return outpath


def plot_density_log2(df_melted: pd.DataFrame, outdir: Path) -> Path:
    """
    Densidad de conteos log2 por muestra, facetado por condición.
    """
    graf = sns.displot(
        data=df_melted,
        x="Value",
        hue="Sample",
        col="Cond",
        kind="kde",
        fill=True,
        alpha=0.1,
    )
    graf.figure.suptitle(
        "Estimación de densidad de los conteos log₂ por condición y muestra",
        y=1.02,
    )
    graf.set_axis_labels("log₂(counts + 1)", "Densidad")
    graf.tight_layout()

    outpath = outdir / "02_density_log2.png"
    graf.savefig(outpath, dpi=300)
    plt.close("all")
    return outpath


def run_pca_and_plot(
    df_log2: pd.DataFrame,
    conditions: pd.Series,
    outdir: Path,
) -> tuple[pd.DataFrame, Path]:
    """
    Corre PCA a partir de la matriz log2 (genes x muestras) y genera un scatterplot.
    """
    expression_df = df_log2.T  # filas: muestras, columnas: genes
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(expression_df)

    pca_df = pd.DataFrame(
        data=pcs,
        columns=["PC1", "PC2"],
        index=expression_df.index,
    )
    # Alinea con el índice (muestras)
    pca_df["condition"] = conditions.loc[pca_df.index].values

    # Plot
    plt.figure(figsize=(6, 6))
    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue="condition",
        style="condition",
        s=100,
    )
    plt.title("PCA de expresión génica")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}%)")
    plt.tight_layout()

    outpath = outdir / "03_pca_samples.png"
    plt.savefig(outpath, dpi=300)
    plt.close()

    return pca_df, outpath


# ============================================================
# Filtrado por CPM y preparación de DESeq2
# ============================================================

def filter_low_counts_cpm(
    df_counts_renamed: pd.DataFrame,
    conditions: pd.Series,
    k: int = 10,
) -> pd.DataFrame:
    """
    Aplica el filtro tipo edgeR:
    - Calcula CPM
    - Umbral U = k * 1e6 / Lmin
    - Mantiene genes con CPM >= U en al menos n muestras (n = mínimo # réplicas por condición)
    """
    libsizes = df_counts_renamed.sum(axis=0)
    cpm = df_counts_renamed.div(libsizes, axis=1) * 1e6

    Lmin = libsizes.min()
    n = conditions.value_counts().min()

    U = k * 1e6 / Lmin
    keep = (cpm >= U).sum(axis=1) >= n
    counts_filt = df_counts_renamed.loc[keep].copy()

    return counts_filt


def build_deseq_dataset(
    counts_filt: pd.DataFrame,
    conditions: pd.Series,
) -> tuple[DeseqDataSet, pd.DataFrame]:
    """
    Construye el DeseqDataSet y la metadata a partir de la matriz filtrada
    y el vector de condiciones (por muestra).
    """
    counts_T = counts_filt.T  # filas: muestras, columnas: genes
    metadata = pd.DataFrame(
        {"condition": conditions.loc[counts_T.index].values},
        index=counts_T.index,
    )

    dds = DeseqDataSet(
        counts=counts_T,
        metadata=metadata,
        design_factors="condition",
    )
    return dds, metadata


# ============================================================
# Contrastes, DESeq2 y clasificación de genes DE
# ============================================================

def get_default_contrasts(metadata: pd.DataFrame) -> dict[str, list[str]]:
    """
    Define un conjunto de contrastes por defecto, pero sólo mantiene
    aquellos cuyas condiciones estén presentes en la metadata.
    """
    present = set(metadata["condition"].unique())

    # Ajusta esta tabla según tus condiciones reales
    candidate_contrasts: dict[str, tuple[str, str]] = {
        "Nod_vs_Control": ("Control_Nod", "Control"),
        "Myc_vs_Control": ("Control_Myc", "Control"),
        "RNAi_vs_Control": ("PvRbohB-RNAi", "Control"),
        "RNAiNod_vs_RNAi": ("PvRbohB-RNAi_Nod", "PvRbohB-RNAi"),
        "RNAiMyc_vs_RNAi": ("PvRbohB-RNAi_Myc", "PvRbohB-RNAi"),
    }

    contrasts: dict[str, list[str]] = {}
    for name, (level, ref) in candidate_contrasts.items():
        if level in present and ref in present:
            contrasts[name] = ["condition", level, ref]

    return contrasts


def run_deseq2_with_contrasts(
    dds: DeseqDataSet,
    contrasts: dict[str, list[str]],
) -> dict[str, pd.DataFrame]:
    """
    Ejecuta DESeq2 y calcula estadísticos para cada contraste.
    """
    print("Ejecutando deseq2()...")
    dds.deseq2()
    print("Listo.")

    results: dict[str, pd.DataFrame] = {}

    for name, contrast in contrasts.items():
        print(f"Calculando contraste {name}: {contrast}")
        ds = DeseqStats(dds, contrast=contrast)
        ds.summary()
        results[name] = ds.results_df.copy()

    return results


def classify_de(
    df: pd.DataFrame,
    padj_thr: float = 0.01,
    lfc_thr: float = 2.0,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """
    Añade columna DE_class con categorías:
      - 'UP'   : padj <= padj_thr y log2FC >=  lfc_thr
      - 'DOWN' : padj <= padj_thr y log2FC <= -lfc_thr
      - 'NS'   : resto
    """
    cond_up = (df["padj"] <= padj_thr) & (df["log2FoldChange"] >= lfc_thr)
    cond_down = (df["padj"] <= padj_thr) & (df["log2FoldChange"] <= -lfc_thr)

    df = df.copy()
    df["DE_class"] = np.select(
        [cond_up, cond_down],
        ["UP", "DOWN"],
        default="NS",
    )
    return df, cond_up, cond_down


def classify_all_results(
    results: dict[str, pd.DataFrame],
    padj_thr: float = 0.01,
    lfc_thr: float = 2.0,
) -> dict[str, pd.DataFrame]:
    """
    Aplica classify_de() a todos los contrastes.
    """
    classified: dict[str, pd.DataFrame] = {}
    for name, df in results.items():
        classified_df, _, _ = classify_de(df, padj_thr=padj_thr, lfc_thr=lfc_thr)
        classified[name] = classified_df
    return classified


# ============================================================
# Volcano plots y export de resultados
# ============================================================

def plot_volcano_per_contrast(
    classified_results: dict[str, pd.DataFrame],
    outdir: Path,
    padj_thr: float = 0.01,
    lfc_thr: float = 2.0,
) -> list[Path]:
    """
    Genera volcano plots para cada contraste clasificado.
    """
    volcano_palette = {
        "UP": "red",
        "DOWN": "forestgreen",
        "NS": "darkgray",
    }

    outpaths: list[Path] = []

    for condition, df in classified_results.items():
        df = df.copy()
        df["padj"] = df["padj"].fillna(1.0)
        df["log10Neg"] = -np.log10(df["padj"].clip(lower=1e-300))

        plt.figure(figsize=(4, 6))
        sns.scatterplot(
            data=df,
            x="log2FoldChange",
            y="log10Neg",
            hue="DE_class",
            palette=volcano_palette,
            s=10,
            edgecolor=None,
        )

        plt.axvline(lfc_thr, ls="--", lw=1, c="black")
        plt.axvline(-lfc_thr, ls="--", lw=1, c="black")
        plt.axhline(-np.log10(padj_thr), ls="--", lw=1, c="black")

        plt.xlabel("log2(Fold Change)")
        plt.ylabel("-log10(padj)")
        plt.title(f"Volcano Plot - {condition}")
        plt.grid(True, alpha=0.2)
        plt.legend(title="DE class", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()

        outpath = outdir / f"volcano_{condition}.png"
        plt.savefig(outpath, dpi=300)
        plt.close()

        outpaths.append(outpath)

    return outpaths


def export_results_tsv(
    classified_results: dict[str, pd.DataFrame],
    outdir: Path,
) -> list[Path]:
    """
    Exporta los resultados DESeq2 (ya clasificados) a TSV por contraste.
    """
    outpaths: list[Path] = []
    for name, df in classified_results.items():
        outpath = outdir / f"{name}_deseq2_results.tsv"
        df.to_csv(outpath, sep="\t")
        outpaths.append(outpath)
    return outpaths


def export_background_genes(
    results: dict[str, pd.DataFrame],
    outdir: Path,
) -> Path | None:
    """
    Guarda la tabla de genes del primer contraste como 'background_genes.tsv'
    para usos posteriores (p.ej. enriquecimiento).
    """
    if not results:
        return None

    first_key = list(results.keys())[0]
    background_genes = results[first_key]
    outpath = outdir / "background_genes.tsv"
    background_genes.to_csv(outpath, sep="\t")
    return outpath

# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Análisis de expresión diferencial con PyDESeq2 y generación de "
            "gráficas (boxplot, densidad, PCA, volcanos)."
        )
    )

    parser.add_argument(
        "--bioproject",
        default="482464",
        help="ID del BioProject (por defecto: 482464).",
    )
    parser.add_argument(
        "--counts",
        help=(
            "Ruta a la matriz de conteos. Si no se da, se asume "
            "results/<BIOPROJECT>/All_counts_rbohb.txt"
        ),
    )
    parser.add_argument(
        "--srr-ids",
        help=(
            "Ruta a srr_ids.txt. Si no se da, se asume "
            "data/<BIOPROJECT>/srr_ids.txt"
        ),
    )
    parser.add_argument(
        "--outdir",
        help=(
            "Directorio de salida. Por defecto, diffexp_directory(project_root())/<BIOPROJECT>"
        ),
    )
    parser.add_argument(
        "--padj-thr",
        type=float,
        default=0.01,
        help="Umbral de padj para clasificar genes DE (default: 0.01).",
    )
    parser.add_argument(
        "--lfc-thr",
        type=float,
        default=2.0,
        help="Umbral de |log2FC| para clasificar genes DE (default: 2).",
    )

    args = parser.parse_args()

    counts_path, srr_ids_path, outdir = get_default_paths(
        bioproject=args.bioproject,
        counts_path=args.counts,
        srr_ids_path=args.srr_ids,
        outdir=args.outdir,
    )

    # 1) Leer conteos y renombrar columnas
    df_counts_renamed = load_counts_and_rename(counts_path, srr_ids_path)
    conditions = build_conditions_series(df_counts_renamed)

    # 2) Transformación log2 y plots exploratorios
    df_log2, df_melted = compute_log2_and_long_format(df_counts_renamed)
    plot_boxplot_log2(df_melted, outdir)
    plot_density_log2(df_melted, outdir)
    run_pca_and_plot(df_log2, conditions, outdir)

    # 3) Filtrado CPM
    counts_filt = filter_low_counts_cpm(df_counts_renamed, conditions, k=10)

    # 4) Preparar y correr DESeq2
    dds, metadata = build_deseq_dataset(counts_filt, conditions)
    contrasts = get_default_contrasts(metadata)
    results = run_deseq2_with_contrasts(dds, contrasts)

    # 5) Clasificación, export de tablas y volcános
    classified_results = classify_all_results(
        results,
        padj_thr=args.padj_thr,
        lfc_thr=args.lfc_thr,
    )

    export_results_tsv(classified_results, outdir)
    export_background_genes(results, outdir)
    plot_volcano_per_contrast(
        classified_results,
        outdir,
        padj_thr=args.padj_thr,
        lfc_thr=args.lfc_thr,
    )

    print(f"Análisis completado. Resultados y gráficas en: {outdir}")

if __name__ == "__main__":
    main()
