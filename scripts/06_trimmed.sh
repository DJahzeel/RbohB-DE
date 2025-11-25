#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 05_trimmed.sh
# Descripción:
#   Realiza el recorte (trimming) de lecturas FASTQ usando cutadapt, tanto
#   para datos single-end (SE) como paired-end (PE), según la información
#   contenida en samples.tsv.
#
#   Aplica parámetros por defecto de calidad y longitud mínima, pensados para
#   datos de tipo NextSeq (opción --nextseq-trim=20) y un umbral de calidad
#   Q20, descartando lecturas recortadas por debajo de longitud 20.
#
# Flujo general:
#   1. Detecta la raíz del proyecto (PROJECT_ROOT) desde la ubicación del script.
#   2. Lee el archivo samples.tsv para obtener:
#        - ID de muestra (sample)
#        - nombre de archivo de R1
#        - nombre de archivo de R2 (opcional)
#   3. Verifica que existan los archivos FASTQ en:
#        data/<BIOPROJECT_ID>/fastq/
#   4. Ejecuta cutadapt:
#        - Si hay R2 → trimming PE, genera *_clean_1.fastq.gz y *_clean_2.fastq.gz
#        - Si no hay R2 → trimming SE, genera *_clean.fastq.gz
#   5. Deja todos los FASTQ recortados en:
#        results/<BIOPROJECT_ID>/trimmed/
#
# Uso:
#   ./05_trimmed.sh
#   ./05_trimmed.sh --threads 8
#   ./05_trimmed.sh --bioproject 482464
#   ./05_trimmed.sh --se "--nextseq-trim=15 -q 20 -m 30"
#   ./05_trimmed.sh --pe "--nextseq-trim=20 -q 20 -m 20"
#
# Flags:
#   --threads/-t   : número de hilos para cutadapt
#   --se           : argumentos para cutadapt en modo single-end
#   --pe           : argumentos para cutadapt en modo paired-end
#   --bioproject/-b: ID de subcarpeta en data/ y results/
#
# Estructura esperada:
#   PROJECT_ROOT/
#     ├─ samples.tsv
#     ├─ data/
#     │    └─ <BIOPROJECT_ID>/
#     │          └─ fastq/
#     │               ├─ *.fastq.gz
#     └─ results/
#          └─ <BIOPROJECT_ID>/
#               └─ trimmed/      # salida de este script
#
# Dependencias:
#   - cutadapt instalado y en el PATH.
# -----------------------------------------------------------------------------

set -euo pipefail
# set -e  : termina el script si algún comando devuelve un error
# set -u  : error si se usa una variable no definida
# pipefail: si un comando en un pipe falla, el pipe completo falla

# Directorio donde está este script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Raíz del proyecto: un nivel arriba de SCRIPT_DIR
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

SAMPLES_TSV="${PROJECT_ROOT}/samples.tsv"

# Valores por defecto. Parámetros modificables vía flags.
THREADS=4
CUTADAPT_ARGS_SE="--nextseq-trim=20 -q 20 -m 20"
CUTADAPT_ARGS_PE="--nextseq-trim=20 -q 20 -m 20"
BIOPROJECT_ID="482464"   # Default

# ---------- Parseo de flags ----------
# Permite sobreescribir:
#   - número de hilos (--threads/-t)
#   - parámetros de cutadapt para SE (--se)
#   - parámetros de cutadapt para PE (--pe)
#   - bioproyecto (--bioproject/-b)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads|-t)
      THREADS="$2"
      shift 2
      ;;
    --se)
      CUTADAPT_ARGS_SE="$2"
      shift 2
      ;;
    --pe)
      CUTADAPT_ARGS_PE="$2"
      shift 2
      ;;
    --bioproject|-b)
      BIOPROJECT_ID="$2"
      shift 2
      ;;
    --help|-h)
      echo "Uso: $0 [--threads N] [--se \"ARGS_SE\"] [--pe \"ARGS_PE\"] [--bioproject ID]"
      echo
      echo "Defaults:"
      echo "  --threads ${THREADS}"
      echo "  --se  \"${CUTADAPT_ARGS_SE}\""
      echo "  --pe  \"${CUTADAPT_ARGS_PE}\""
      echo "  --bioproject 482464"
      echo
      echo "Se espera que r1 y r2 en samples.tsv sean solo nombres de archivo,"
      echo "y se buscarán en: data/ID/fastq/<archivo.fastq.gz>"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

# Directorio donde se guardarán los FASTQ recortados
TRIM_OUT="${PROJECT_ROOT}/results/${BIOPROJECT_ID}/trimmed"
# Directorio donde se encuentran los FASTQ crudos
FASTQ_BASE="${PROJECT_ROOT}/data/${BIOPROJECT_ID}/fastq"

# Resumen de parámetros para el usuario
echo "[INFO] 10_trimmed.sh"
echo "  - PROJECT_ROOT:     ${PROJECT_ROOT}"
echo "  - SAMPLES_TSV:      ${SAMPLES_TSV}"
echo "  - FASTQ_BASE:       ${FASTQ_BASE}"
echo "  - TRIM_OUT:         ${TRIM_OUT}"
echo "  - THREADS:          ${THREADS}"
echo "  - CUTADAPT_ARGS_SE: ${CUTADAPT_ARGS_SE}"
echo "  - CUTADAPT_ARGS_PE: ${CUTADAPT_ARGS_PE}"

# Verificación de existencia de samples.tsv
if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: No encuentro samples.tsv en ${SAMPLES_TSV}" >&2
  exit 1
fi

# Crea el directorio de salida si no existe
mkdir -p "${TRIM_OUT}"

# Lectura línea a línea de samples.tsv
# Se espera: sample_id \t r1 \t r2(opcional)
while IFS=$'\t' read -r sample r1 r2; do
  # Saltar líneas de comentario o encabezados
  [[ "$sample" =~ ^# ]] && continue
  [[ "$sample" == "sample" || "$sample" == "sample_id" ]] && continue
  # Ignorar líneas vacías o sin R1
  [[ -z "${sample}" || -z "${r1}" ]] && continue
  # Ruta al FASTQ de R1
  in1="${FASTQ_BASE}/${r1}"
  if [[ ! -f "${in1}" ]]; then
    echo "ERROR: No encuentro FASTQ: ${in1}" >&2
    exit 1
  fi
  # Base de nombre para archivos de salida
  base="${TRIM_OUT}/${sample}"

  # Verificar si hay R2 para decidir entre SE o PE
  if [[ -n "${r2:-}" ]]; then
    in2="${FASTQ_BASE}/${r2}"
    if [[ ! -f "${in2}" ]]; then
      echo "ERROR: No encuentro FASTQ: ${in2}" >&2
      exit 1
    fi

    echo "[INFO] PE sample=${sample}"
    cutadapt ${CUTADAPT_ARGS_PE} -j "${THREADS}" \
      -o "${base}_clean_1.fastq.gz" -p "${base}_clean_2.fastq.gz" \
      "${in1}" "${in2}"
  # Si no hay R2, es single-end
  else
    echo "[INFO] SE sample=${sample}"
    cutadapt ${CUTADAPT_ARGS_SE} -j "${THREADS}" \
      -o "${base}_clean.fastq.gz" \
      "${in1}"
  fi

done < "${SAMPLES_TSV}"

echo "[OK] Trimming completado. Archivos en: ${TRIM_OUT}"