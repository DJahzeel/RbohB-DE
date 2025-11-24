#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 09_bwa_align.sh
# Descripción:
#   Alinea las lecturas trimmed (recortadas) contra un genoma de referencia
#   indexado con BWA, usando bwa mem + samtools sort + samtools index.
#
#   El script:
#     - Lee la lista de muestras desde samples.tsv.
#     - Detecta automáticamente el índice de BWA a partir del único archivo
#       *.bwt que encuentre en data/references/bwa/.
#     - Para cada muestra:
#         * Si tiene R1 y R2 → asume datos paired-end.
#         * Si solo tiene R1 → asume datos single-end.
#       En ambos casos exige que existan los FASTQ trimmed generados por el
#       script de trimming (p.ej. 05_trimmed.sh).
#     - Genera archivos BAM ordenados por coordenadas y su índice (.bai).
#
# Flujo general:
#   1. Detectar PROJECT_ROOT a partir del directorio del script.
#   2. Verificar existencia de:
#        - samples.tsv
#        - directorio con índices BWA: data/references/bwa/
#   3. Detectar automáticamente el prefijo del índice BWA (.bwt único).
#   4. Crear results/<BIOPROJECT_ID>/bam/ si no existe.
#   5. Recorrer samples.tsv:
#        - Construir nombres de FASTQ trimmed esperados:
#            PE: results/<BIOPROJECT_ID>/trimmed/<sample>_clean_1.fastq.gz
#                results/<BIOPROJECT_ID>/trimmed/<sample>_clean_2.fastq.gz
#            SE: results/<BIOPROJECT_ID>/trimmed/<sample>_clean.fastq.gz
#        - Ejecutar:
#            bwa mem ... | samtools sort ... -> sample.bam
#            samtools index sample.bam
#
# Uso:
#   ./09_bwa_align.sh
#   ./09_bwa_align.sh --bioproject 482464
#
# Flags:
#   --bioproject/-b : ID de subcarpeta bajo results/ y data/.
#
# Estructura esperada (minima):
#   PROJECT_ROOT/
#     ├─ samples.tsv
#     ├─ data/
#     │    └─ references/
#     │         └─ bwa/
#     │              ├─ <prefix>.bwt   # índice BWA (y archivos asociados)
#     └─ results/
#          └─ <BIOPROJECT_ID>/
#               ├─ trimmed/            # FASTQ recortados
#               └─ bam/                # salida de este script
#
# Dependencias:
#   - bwa instalado y en el PATH.
#   - samtools instalado y en el PATH.
#   - Índices BWA ya generados (p.ej. con index_bwa.sh).
#   - FASTQ trimmed generados previo al alineamiento.
# -----------------------------------------------------------------------------
set -euo pipefail
# set -e  : termina el script si algún comando devuelve error
# set -u  : error si se usa una variable no definida
# pipefail: si algún comando en un pipe falla, el pipe completo falla

# Detectar ruta del script y raíz del proyecto
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Raíz del proyecto: un nivel arriba del directorio del script
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Valores por defecto. 
THREADS=4
BIOPROJECT_ID="482464" % # Parametro modificable vía flag  

# ---------- Parseo de flags ----------
# Por ahora solo permite cambiar el BIOPROJECT_ID.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bioproject|-b)
      BIOPROJECT_ID="$2"
      shift 2
      ;;
    --help|-h)
      echo "Uso: $0  [--bioproject ID]"
      echo
      echo "  --bioproject/-b ID de subcarpeta bajo results/ (default: 482464)"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

# Rutas fijas basadas en BIOPROJECT_ID
TRIM_OUT="${PROJECT_ROOT}/results/${BIOPROJECT_ID}/trimmed"
BAM_OUT="${PROJECT_ROOT}/results/${BIOPROJECT_ID}/bam"
SAMPLES_TSV="${PROJECT_ROOT}/samples.tsv"
REF_BWA_DIR="${PROJECT_ROOT}/data/references/bwa"

# Comprobar existencia de samples.tsv
if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: No encuentro samples.tsv en ${SAMPLES_TSV}" >&2
  exit 1
fi

# Comprobar existencia de directorio de índices BWA
if [[ ! -d "${REF_BWA_DIR}" ]]; then
  echo "ERROR: No existe el directorio de índices BWA: ${REF_BWA_DIR}" >&2
  echo "       Corre primero el script de indexado (02_index_bwa.sh)." >&2
  exit 1
fi

# Detectar automáticamente el prefijo BWA (único .bwt en data/reference/bwa)
mapfile -t BWT_FILES < <(find "${REF_BWA_DIR}" -maxdepth 1 -type f -name "*.bwt" | sort)

if [[ "${#BWT_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: No se encontraron archivos .bwt en ${REF_BWA_DIR}" >&2
  echo "       Asegúrate de haber creado los índices con 02_index_bwa.sh." >&2
  exit 1
elif [[ "${#BWT_FILES[@]}" -gt 1 ]]; then
  echo "ERROR: Se encontraron múltiples índices BWA (.bwt) en ${REF_BWA_DIR}:" >&2
  printf '  - %s\n' "${BWT_FILES[@]}" >&2
  echo "       Deja solo el índice que quieras usar o ajusta el script." >&2
  exit 1
fi
# Extraer prefijo del archivo .bwt
BWT_FILE="${BWT_FILES[0]}"
BWA_INDEX_PREFIX="${BWT_FILE%.bwt}"

# Crear carpeta de salida si no existe
mkdir -p "${BAM_OUT}"

# Bucle principal: recorrer samples.tsv y alinear cada muestra
while IFS=$'\t' read -r sample r1 r2; do
  # Saltar comentarios o líneas vacías
  [[ "$sample" =~ ^# ]] && continue
  [[ -z "${sample}" || -z "${r1}" ]] && continue

  if [[ -n "${r2:-}" ]]; then
    # ---------------------------------------------------------------
    # Caso PAIR-END: exige que existan los FASTQ trimmed PE
    # ---------------------------------------------------------------
    r1_trim="${TRIM_OUT}/${sample}_clean_1.fastq.gz"
    r2_trim="${TRIM_OUT}/${sample}_clean_2.fastq.gz"

    if [[ ! -f "${r1_trim}" || ! -f "${r2_trim}" ]]; then
      echo "ERROR: No se encontraron FASTQ trimmed PE para sample=${sample}" >&2
      echo "  Esperado:" >&2
      echo "    ${r1_trim}" >&2
      echo "    ${r2_trim}" >&2
      echo "  Corre primero el script de trimming (p.ej. 10_trimmed.sh)." >&2
      exit 1
    fi
    # Alineamiento PE: bwa mem → samtools sort → BAM ordenado
    bwa mem -t "${THREADS}" "${BWA_INDEX_PREFIX}" "${r1_trim}" "${r2_trim}" \
      | samtools sort -@ "${THREADS}" -o "${BAM_OUT}/${sample}.bam" -
  else
    # ---------------------------------------------------------------
    # Caso SINGLE-END: exige que exista el FASTQ trimmed SE
    # ---------------------------------------------------------------
    r1_trim="${TRIM_OUT}/${sample}_clean.fastq.gz"

    if [[ ! -f "${r1_trim}" ]]; then
      echo "ERROR: No se encontró FASTQ trimmed SE para sample=${sample}" >&2
      echo "  Esperado:" >&2
      echo "    ${r1_trim}" >&2
      echo "  Corre primero el script de trimming (p.ej. 10_trimmed.sh)." >&2
      exit 1
    fi
    # Alineamiento SE: bwa mem → samtools sort → BAM ordenado
    bwa mem -t "${THREADS}" "${BWA_INDEX_PREFIX}" "${r1_trim}" \
      | samtools sort -@ "${THREADS}" -o "${BAM_OUT}/${sample}.bam" -
  fi
  # Indexar el BAM resultante
  samtools index -@ "${THREADS}" "${BAM_OUT}/${sample}.bam"
# Se usa awk para leer solo las tres primeras columnas (sample, r1, r2)
done < <(awk 'BEGIN{FS=OFS="\t"} !/^#/{print $1,$2,$3}' "${SAMPLES_TSV}")