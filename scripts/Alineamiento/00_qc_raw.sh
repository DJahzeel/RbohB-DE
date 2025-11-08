#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

SAMPLES_TSV="${PROJECT_ROOT}/samples.tsv"
RESULTS_DIR="${PROJECT_ROOT}/results"
QC_DIR="${RESULTS_DIR}/qc"

THREADS=4
BIOPROJECT_ID="482464"   # Default

# ---------- Parseo de flags ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads|-t)
      THREADS="$2"
      shift 2
      ;;
    --bioproject|-b)
      BIOPROJECT_ID="$2"
      shift 2
      ;;
    --help|-h)
      echo "Uso: $0 [--threads N] [--bioproject ID]"
      echo
      echo "  --threads/-t    N hilos para FastQC (default: ${THREADS})"
      echo "  --bioproject/-b ID de subcarpeta bajo data/ (default: 482464)"
      echo "                  FASTQ esperados en: data/ID/fastq/<archivo.fastq.gz>"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

FASTQ_BASE="${PROJECT_ROOT}/data/${BIOPROJECT_ID}/fastq"

echo "[INFO] 00_qc_raw.sh"
echo "  - PROJECT_ROOT: ${PROJECT_ROOT}"
echo "  - SAMPLES_TSV:  ${SAMPLES_TSV}"
echo "  - FASTQ_BASE:   ${FASTQ_BASE}"
echo "  - QC_DIR:       ${QC_DIR}"
echo "  - THREADS:      ${THREADS}"

if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: No encuentro samples.tsv en ${SAMPLES_TSV}" >&2
  exit 1
fi

mkdir -p "${RESULTS_DIR}" "${QC_DIR}"

FQS=()
while IFS=$'\t' read -r sample r1 r2; do
  # Saltar comentarios o encabezado
  [[ "$sample" =~ ^# ]] && continue
  [[ "$sample" == "sample" || "$sample" == "sample_id" ]] && continue
  [[ -z "${sample}" || -z "${r1}" ]] && continue

  fq1="${FASTQ_BASE}/${r1}"
  if [[ ! -f "${fq1}" ]]; then
    echo "ERROR: No encuentro FASTQ: ${fq1}" >&2
    exit 1
  fi
  FQS+=( "${fq1}" )

  if [[ -n "${r2:-}" ]]; then
    fq2="${FASTQ_BASE}/${r2}"
    if [[ ! -f "${fq2}" ]]; then
      echo "ERROR: No encuentro FASTQ: ${fq2}" >&2
      exit 1
    fi
    FQS+=( "${fq2}" )
  fi
done < "${SAMPLES_TSV}"

if [[ "${#FQS[@]}" -eq 0 ]]; then
  echo "ERROR: No se encontraron FASTQ a partir de samples.tsv" >&2
  exit 1
fi

echo "[INFO] FASTQ a analizar:"
printf '  - %s\n' "${FQS[@]}"

fastqc -t "${THREADS}" -o "${QC_DIR}" "${FQS[@]}"

echo "[OK] QC RAW completado en: ${QC_DIR}"
