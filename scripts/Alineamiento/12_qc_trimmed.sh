#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

RESULTS_DIR="${PROJECT_ROOT}/results"
TRIM_DIR="${RESULTS_DIR}/trimmed"
QC_TRIM_DIR="${RESULTS_DIR}/qc_trimmed"

THREADS=4

if [[ ! -d "${TRIM_DIR}" ]]; then
  echo "ERROR: No existe el directorio con FASTQ trimmed: ${TRIM_DIR}" >&2
  echo "       Asegúrate de haber corrido primero trim.sh." >&2
  exit 1
fi

# Crear carpeta de salida si no existe
if [[ ! -d "${QC_TRIM_DIR}" ]]; then
  echo "[INFO] Creando directorio de salida para QC trimmed: ${QC_TRIM_DIR}"
  mkdir -p "${QC_TRIM_DIR}"
fi

# Construir lista de FASTQ trimmed
echo "[INFO] Buscando FASTQ trimmed en ${TRIM_DIR}"
mapfile -t FQS < <(find "${TRIM_DIR}" -type f -name "*.fastq.gz" | sort)

if [[ "${#FQS[@]}" -eq 0 ]]; then
  echo "ERROR: No se encontraron archivos *.fastq.gz en ${TRIM_DIR}" >&2
  exit 1
fi

echo "[INFO] Se encontraron ${#FQS[@]} archivos FASTQ trimmed para QC:"
printf '  - %s\n' "${FQS[@]}"

# 7) Ejecutar FastQC (paralelo dentro del nodo con xargs -P THREADS)
printf "%s\n" "${FQS[@]}" \
  | xargs -n 1 -P "${THREADS}" fastqc -o "${QC_TRIM_DIR}"

echo "QC de trimmed completado. Resultados en: ${QC_TRIM_DIR}"

