#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
KEGG KO enrichment desde un archivo TSV/CSV que contiene:
- Una columna con IDs KO
- Una columna con valores de padj (p-valor ajustado)

Uso típico:

    py kegg_enrich_padj.py ^
        -i RNAi_vs_Control_deseq2_results.with_KO.tsv ^
        --ko-column ko_id ^
        --padj-column padj ^
        --alpha 0.05 ^
        -o kegg_enrichment.tsv ^
        --out-plot kegg_enrichment.png
"""

import argparse
import csv
import pandas as pd
import requests
from collections import defaultdict
from math import log10
import time

from scipy.stats import fisher_exact
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt


# ============================================================================
# 1) KEGG API
# ============================================================================

def obtener_rutas_desde_ko(ko_id):
    """
    Devuelve las rutas KEGG asociadas a un KO usando KEGG REST.
    """
    if not ko_id.startswith("ko:"):
        ko_id = f"ko:{ko_id}"

    url = f"http://rest.kegg.jp/link/pathway/{ko_id}"
    r = requests.get(url)

    if r.status_code != 200:
        return []

    pathways = []
    for line in r.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) > 1:
            pathways.append(parts[1])  # ej. 'path:map00010'
    return list(set(pathways))


def obtener_nombre_rutas(pathway_id):
    """
    Obtiene el nombre descriptivo de una vía KEGG.
    """
    if not pathway_id.startswith("path:"):
        pathway_id = f"path:{pathway_id}"

    url = f"http://rest.kegg.jp/get/{pathway_id}"
    r = requests.get(url)

    if r.status_code != 200:
        return "Nombre desconocido"

    for line in r.text.split("\n"):
        if line.startswith("NAME"):
            return line.split("NAME")[1].strip()
    return "Nombre desconocido"


# ============================================================================
# 2) Lectura de KO + padj desde archivo único
# ============================================================================

def leer_kos_y_padj(ruta, ko_col, padj_col, alpha):
    """
    Lee un TSV/CSV y devuelve:
    - lista de KO de background (todos los KO únicos)
    - lista de KO de interés (KO con padj < alpha, únicos)

    Se asume que la tabla tiene al menos:
    - ko_col  : IDs KO (con o sin prefijo 'ko:')
    - padj_col: p-valor ajustado
    """
    # Intentar como TSV, si falla CSV
    try:
        df = pd.read_csv(ruta, sep="\t")
    except Exception:
        df = pd.read_csv(ruta, sep=",")

    cols_lower = {c.lower(): c for c in df.columns}

    # Resolver columnas ignorando mayúsc/minúsculas
    ko_key = ko_col.lower()
    padj_key = padj_col.lower()

    if ko_key not in cols_lower:
        raise ValueError(
            f"No se encontró la columna de KO '{ko_col}' en {ruta}.\n"
            f"Columnas disponibles: {list(df.columns)}"
        )
    if padj_key not in cols_lower:
        raise ValueError(
            f"No se encontró la columna de padj '{padj_col}' en {ruta}.\n"
            f"Columnas disponibles: {list(df.columns)}"
        )

    real_ko_col   = cols_lower[ko_key]
    real_padj_col = cols_lower[padj_key]

    # Limpiar KO
    df = df.copy()
    df[real_ko_col] = (
        df[real_ko_col]
        .astype(str)
        .str.strip()
        .str.replace("ko:", "", regex=False)
    )

    # Quitar KO vacíos / NaN
    df = df[df[real_ko_col].notna() & (df[real_ko_col] != "")]

    # Background = TODOS los KO únicos presentes en la tabla
    kos_bg = sorted(df[real_ko_col].unique().tolist())

    # Interés = KO con padj < alpha (y padj no nulo)
    df_sig = df[df[real_padj_col].notna() & (df[real_padj_col] < alpha)]
    kos_int = sorted(df_sig[real_ko_col].unique().tolist())

    return kos_bg, kos_int


# ============================================================================
# 3) Construcción KO ↔ pathways (batch)
# ============================================================================

def construir_mapas_kos(kos, cache=None, verbose=True, chunk_size=20):
    """
    Construye un dict pathway -> set(KO) usando KEGG, haciendo consultas en lote.

    kos : lista de KO (strings, con o sin prefijo 'ko:')
    cache : dict opcional KO -> [pathways]
    chunk_size : cuántos KO mandar por petición a KEGG (20 es razonable).
    """
    if cache is None:
        cache = {}

    # Normalizar a 'ko:Kxxxxx'
    kos_norm = [f"ko:{k}" if not k.startswith("ko:") else k for k in kos]
    # Quitar duplicados, conservando orden
    kos_norm = list(dict.fromkeys(kos_norm))

    total = len(kos_norm)

    # Separar los que ya están en cache de los que no
    kos_uncached = [ko for ko in kos_norm if ko not in cache]

    if verbose:
        print(f"[INFO] KO totales en background (únicos): {total}")
        print(f"[INFO] KO que requieren consulta a KEGG (no están en cache): {len(kos_uncached)}")
        print(f"[INFO] Usando consultas en lote con chunk_size = {chunk_size}")

    # Consultar KEGG por chunks
    for start in range(0, len(kos_uncached), chunk_size):
        chunk = kos_uncached[start:start + chunk_size]
        if not chunk:
            continue

        # Construimos la URL de KEGG con varios KO
        # Ej: link/pathway/ko:K00001+ko:K00002+...
        ko_str = "+".join(chunk)
        url = f"http://rest.kegg.jp/link/pathway/{ko_str}"

        if verbose:
            print(f"[INFO] Consultando KEGG para KO {start+1}–{start+len(chunk)} "
                  f"de {len(kos_uncached)}...")

        r = requests.get(url)

        # Inicializar en cache como listas vacías, por si algún KO no sale en la respuesta
        for ko in chunk:
            cache.setdefault(ko, [])

        if r.status_code == 200 and r.text.strip():
            for line in r.text.strip().split("\n"):
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                ko_id, pathway = parts[0].strip(), parts[1].strip()
                # Guardamos en cache
                cache.setdefault(ko_id, []).append(pathway)

        # Pausa suave para no spamear KEGG
        time.sleep(0.2)

    # Construir el dict pathway -> set(KO)
    path2ko = defaultdict(set)
    for ko in kos_norm:
        rutas = cache.get(ko, [])
        for p in rutas:
            path2ko[p].add(ko)

    if verbose:
        print(f"[INFO] Construcción KO → vías KEGG terminada.")
        print(f"[INFO] Total de vías distintas: {len(path2ko)}")

    return path2ko, cache


# ============================================================================
# 4) Enriquecimiento KEGG
# ============================================================================

def enriquecimiento_kegg(kos_interes, kos_background, path2ko,
                         min_kos_en_path=5, verbose=True, q_report=0.1):
    """
    ORA (Fisher exact test + FDR) a nivel vía KEGG.
    """
    target = {f"ko:{k}" if not k.startswith("ko:") else k for k in kos_interes}
    bg     = {f"ko:{k}" if not k.startswith("ko:") else k for k in kos_background}

    N = len(bg)
    resultados = []

    if verbose:
        print(f"[INFO] Iniciando enriquecimiento KEGG con {len(path2ko)} vías candidatas...")
        print(f"[INFO] Tamaño universo (background KO únicos): {N}")
        print(f"[INFO] Tamaño conjunto de interés: {len(target)}")

    for j, (pathway, kos) in enumerate(path2ko.items(), start=1):
        kos_bg = kos & bg
        if len(kos_bg) < min_kos_en_path:
            continue

        a = len(target & kos_bg)
        if a == 0:
            continue

        b = len(target) - a
        c = len(kos_bg - target)
        d = N - a - b - c

        tabla = [[a, b], [c, d]]
        _, pval = fisher_exact(tabla, alternative="greater")

        resultados.append((pathway, a, len(kos_bg), pval))

        if verbose and j % 50 == 0:
            print(f"[INFO] Evaluadas ~{j} vías (no todas pasan filtros).")

    if not resultados:
        if verbose:
            print("[WARN] No hubo vías con KO suficientes para testear.")
        return []

    if verbose:
        print(f"[INFO] Vías con datos suficientes para test estadístico: {len(resultados)}")

    pvals = [r[3] for r in resultados]
    _, qvals, _, _ = smm.multipletests(pvals, method="fdr_bh")

    enriched = []
    for (p, overlap, size, pval), q in zip(resultados, qvals):
        enriched.append({
            "pathway_id": p,
            "overlap": overlap,
            "size": size,
            "pval": pval,
            "qval": q,
            "score": -log10(pval) if pval > 0 else 300
        })

    enriched.sort(key=lambda x: x["pval"])

    if verbose:
        n_signif = sum(e["qval"] <= q_report for e in enriched)
        print(f"[INFO] Vías enriquecidas con q <= {q_report}: {n_signif} "
              f"(de un total de {len(enriched)})")

    return enriched


# ============================================================================
# 5) Dot plot
# ============================================================================

def dotplot_kegg(enriquecidos, out=None,
                 top_n=20, q_threshold=0.1,
                 titulo="KEGG enrichment (KO)"):

    """
    Genera un dot plot tipo clusterProfiler.
    """
    enriched = [e for e in enriquecidos if e["qval"] <= q_threshold]
    if not enriched:
        print("[WARN] No hay vías con q <= umbral, no se genera figura.")
        return

    enriched = enriched[:top_n]
    enriched.reverse()

    names  = [obtener_nombre_rutas(e["pathway_id"]) for e in enriched]
    scores = [e["score"] for e in enriched]
    sizes  = [e["overlap"] * 40 for e in enriched]
    colors = [e["qval"] for e in enriched]

    fig, ax = plt.subplots(figsize=(9, 0.5*len(enriched)+2))
    sc = ax.scatter(scores,
                    range(len(enriched)),
                    s=sizes,
                    c=colors,
                    cmap="viridis")

    ax.set_yticks(range(len(enriched)))
    ax.set_yticklabels(names)
    ax.set_xlabel("-log10(p-valor)")
    ax.set_title(titulo)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("q-valor (FDR)")

    plt.tight_layout()

    if out:
        plt.savefig(out, dpi=300)
        plt.close(fig)
        print(f"[OK] Figura guardada en {out}")
    else:
        plt.show()


# ============================================================================
# 6) Argparse / main
# ============================================================================

def parse_args():
    ap = argparse.ArgumentParser(
        description="Enriquecimiento KEGG (KO) a partir de un archivo con KO + padj."
    )
    ap.add_argument(
        "-i", "--input", required=True,
        help="Archivo TSV/CSV con KO y padj."
    )
    ap.add_argument(
        "--ko-column", default="ko_id",
        help="Nombre de la columna donde están los KO (default: ko_id)."
    )
    ap.add_argument(
        "--padj-column", default="padj",
        help="Nombre de la columna con p-valor ajustado (default: padj)."
    )
    ap.add_argument(
        "--alpha", type=float, default=0.05,
        help="Umbral de padj para definir conjunto de interés (default: 0.05)."
    )
    ap.add_argument(
        "--min-kos", type=int, default=5,
        help="Mínimo de KO del background en una vía para considerarla (default: 5)."
    )
    ap.add_argument(
        "-o", "--out-table", default="kegg_enrichment.tsv",
        help="Archivo de salida con la tabla de enriquecimiento (TSV)."
    )
    ap.add_argument(
        "--out-plot", default="kegg_enrichment.png",
        help="Archivo de salida para la figura (PNG/PDF/etc.)."
    )
    ap.add_argument(
        "--top-n", type=int, default=20,
        help="Máximo de vías a mostrar en el dot plot (default: 20)."
    )
    ap.add_argument(
        "--q-threshold", type=float, default=0.1,
        help="Umbral de q-valor para visualizar en el dot plot (default: 0.1)."
    )
    ap.add_argument(
        "--titulo", default="KEGG enrichment (KO)",
        help="Título del dot plot."
    )
    return ap.parse_args()


def main():
    args = parse_args()

    # 1) Leer KO + padj de un solo archivo
    print(f"[INFO] Leyendo tabla desde: {args.input}")
    kos_bg, kos_int = leer_kos_y_padj(
        ruta=args.input,
        ko_col=args.ko_column,
        padj_col=args.padj_column,
        alpha=args.alpha
    )
    print(f"[INFO] KO en background (únicos): {len(kos_bg)}")
    print(f"[INFO] KO en interés (padj < {args.alpha}): {len(kos_int)}")

    if len(kos_int) == 0:
        print("[WARN] No hay KO con padj < alpha. No se realizará enriquecimiento.")
        return

    # 2) Construir mapa KO → vías
    print("[INFO] Construyendo KO → vías KEGG...")
    path2ko, cache = construir_mapas_kos(kos_bg, verbose=True)

    # 3) Enriquecimiento
    print("[INFO] Calculando enriquecimiento...")
    enr = enriquecimiento_kegg(
        kos_interes=kos_int,
        kos_background=kos_bg,
        path2ko=path2ko,
        min_kos_en_path=args.min_kos,
        verbose=True,
        q_report=args.q_threshold
    )

    if not enr:
        print("[WARN] No se encontraron vías enriquecidas (ni siquiera antes de filtrar por q).")
        return

    print(f"[INFO] Vías candidatas (antes de filtrar por q): {len(enr)}")

    # 4) Escribir tabla
    print(f"[INFO] Escribiendo tabla en: {args.out_table}")
    with open(args.out_table, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pathway_id", "name", "overlap", "size",
                    "pval", "qval", "score"])
        for e in enr:
            w.writerow([
                e["pathway_id"],
                obtener_nombre_rutas(e["pathway_id"]),
                e["overlap"],
                e["size"],
                f"{e['pval']:.3e}",
                f"{e['qval']:.3e}",
                f"{e['score']:.4f}"
            ])

    # 5) Figura
    print("[INFO] Generando dot plot...")
    dotplot_kegg(
        enr,
        out=args.out_plot,
        top_n=args.top_n,
        q_threshold=args.q_threshold,
        titulo=args.titulo
    )

if __name__ == "__main__":
    main()