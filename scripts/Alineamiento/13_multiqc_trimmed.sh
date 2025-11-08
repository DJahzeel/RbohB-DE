#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

RESULTS_DIR="${PROJECT_ROOT}/results"
QC_TRIM_DIR="${RESULTS_DIR}/qc_trimmed" 
MULTIQC_TRIM_DIR="${RESULTS_DIR}/multiqc_trimmed"
REPORT_NAME="multiqc_report.html"

#  Crear solo la carpeta de salida, si no existe
mkdir -p "${MULTIQC_TRIM_DIR}"

#  Ejecutar MultiQC
multiqc "${QC_TRIM_DIR}" -o "${MULTIQC_TRIM_DIR}" -n "${REPORT_NAME}"

# Verificación del reporte
FULL_REPORT_PATH="${MULTIQC_TRIM_DIR}/${REPORT_NAME}"
if [[ -f "${FULL_REPORT_PATH}" ]]; then
    echo "¡Reporte de MultiQC (trimmed) generado con éxito!"
    echo "Reporte disponible en: ${FULL_REPORT_PATH}"
else
    echo "ERROR: Falló la generación del reporte de MultiQC para trimmed."
    exit 1
fi
