#!/usr/bin/env bash
set -euo pipefail

# Detectar ruta del script y raíz del proyecto
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

# Directorios fijos
RESULTS_DIR="${PROJECT_ROOT}/results"
QC_DIR="${RESULTS_DIR}/qc"               # FastQC de los RAW
MULTIQC_DIR="${RESULTS_DIR}/multiqc"     # MultiQC de RAW
REPORT_NAME="multiqc_report.html"

echo "[INFO] MultiQC (RAW)"
echo "  - Entrada QC:   ${QC_DIR}"
echo "  - Salida HTML:  ${MULTIQC_DIR}/${REPORT_NAME}"

# Crear solo la carpeta de salida, si no existe
mkdir -p "${MULTIQC_DIR}"

# Correr MultiQC leyendo de QC_DIR y escribiendo en MULTIQC_DIR
multiqc "${QC_DIR}" -o "${MULTIQC_DIR}" -n "${REPORT_NAME}"

# Verificación
FULL_REPORT_PATH="${MULTIQC_DIR}/${REPORT_NAME}"
if [[ -f "${FULL_REPORT_PATH}" ]]; then
    echo "¡Reporte de MultiQC generado con éxito!"
    echo "Reporte disponible en: ${FULL_REPORT_PATH}"
else
    echo "ERROR: Falló la generación del reporte de MultiQC."
    exit 1
fi

