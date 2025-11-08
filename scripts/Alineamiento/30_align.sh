#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

THREADS=4

TRIM_OUT="${PROJECT_ROOT}/results/trimmed"
BAM_OUT="${PROJECT_ROOT}/results/bam"
SAMPLES_TSV="${PROJECT_ROOT}/samples.tsv"
REF_BWA_DIR="${PROJECT_ROOT}/data/references/bwa"

if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: No encuentro samples.tsv en ${SAMPLES_TSV}" >&2
  exit 1
fi

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

BWT_FILE="${BWT_FILES[0]}"
BWA_INDEX_PREFIX="${BWT_FILE%.bwt}"

# Crear carpeta de salida 
mkdir -p "${BAM_OUT}"

# Bucle principal 
while IFS=$'\t' read -r sample r1 r2; do
  [[ "$sample" =~ ^# ]] && continue
  [[ -z "${sample}" || -z "${r1}" ]] && continue

  if [[ -n "${r2:-}" ]]; then
    # Paired-end: exigir trimmed
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

    bwa mem -t "${THREADS}" "${BWA_INDEX_PREFIX}" "${r1_trim}" "${r2_trim}" \
      | samtools sort -@ "${THREADS}" -o "${BAM_OUT}/${sample}.bam" -
  else
    # Single-end: exigir trimmed
    r1_trim="${TRIM_OUT}/${sample}_clean.fastq.gz"

    if [[ ! -f "${r1_trim}" ]]; then
      echo "ERROR: No se encontró FASTQ trimmed SE para sample=${sample}" >&2
      echo "  Esperado:" >&2
      echo "    ${r1_trim}" >&2
      echo "  Corre primero el script de trimming (p.ej. 10_trimmed.sh)." >&2
      exit 1
    fi

    bwa mem -t "${THREADS}" "${BWA_INDEX_PREFIX}" "${r1_trim}" \
      | samtools sort -@ "${THREADS}" -o "${BAM_OUT}/${sample}.bam" -
  fi

  samtools index -@ "${THREADS}" "${BAM_OUT}/${sample}.bam"
done < <(awk 'BEGIN{FS=OFS="\t"} !/^#/{print $1,$2,$3}' "${SAMPLES_TSV}")
