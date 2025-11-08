#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

SAMPLES_TSV="${PROJECT_ROOT}/samples.tsv"
TRIM_OUT="${PROJECT_ROOT}/results/trimmed"

THREADS=4
CUTADAPT_ARGS_SE="--nextseq-trim=20 -q 20 -m 20"
CUTADAPT_ARGS_PE="--nextseq-trim=20 -q 20 -m 20"
BIOPROJECT_ID="482464"   # Default

# ---------- Parseo de flags ----------
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

FASTQ_BASE="${PROJECT_ROOT}/data/${BIOPROJECT_ID}/fastq"

echo "[INFO] 10_trimmed.sh"
echo "  - PROJECT_ROOT:     ${PROJECT_ROOT}"
echo "  - SAMPLES_TSV:      ${SAMPLES_TSV}"
echo "  - FASTQ_BASE:       ${FASTQ_BASE}"
echo "  - TRIM_OUT:         ${TRIM_OUT}"
echo "  - THREADS:          ${THREADS}"
echo "  - CUTADAPT_ARGS_SE: ${CUTADAPT_ARGS_SE}"
echo "  - CUTADAPT_ARGS_PE: ${CUTADAPT_ARGS_PE}"

if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: No encuentro samples.tsv en ${SAMPLES_TSV}" >&2
  exit 1
fi

mkdir -p "${TRIM_OUT}"

while IFS=$'\t' read -r sample r1 r2; do
  [[ "$sample" =~ ^# ]] && continue
  [[ "$sample" == "sample" || "$sample" == "sample_id" ]] && continue
  [[ -z "${sample}" || -z "${r1}" ]] && continue

  in1="${FASTQ_BASE}/${r1}"
  if [[ ! -f "${in1}" ]]; then
    echo "ERROR: No encuentro FASTQ: ${in1}" >&2
    exit 1
  fi

  base="${TRIM_OUT}/${sample}"

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

  else
    echo "[INFO] SE sample=${sample}"
    cutadapt ${CUTADAPT_ARGS_SE} -j "${THREADS}" \
      -o "${base}_clean.fastq.gz" \
      "${in1}"
  fi

done < "${SAMPLES_TSV}"

echo "[OK] Trimming completado. Archivos en: ${TRIM_OUT}"
