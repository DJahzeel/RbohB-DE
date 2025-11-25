#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 07_multiqc_trimmed.sh
# Descripción:
#   Genera un reporte integrado de calidad usando MultiQC a partir de los
#   resultados de FastQC sobre las lecturas recortadas (trimmed).
#
#   Toma como entrada los reportes de FastQC en:
#       results/<BIOPROJECT_ID>/qc_trimmed/
#   y genera un reporte HTML de MultiQC en:
#       results/<BIOPROJECT_ID>/multiqc_trimmed/multiqc_report.html
#
# Flujo general:
#   1. Detecta la raíz del proyecto (PROJECT_ROOT) desde la ubicación del script.
#   2. Usa BIOPROJECT_ID para localizar:
#        - results/<BIOPROJECT_ID>/qc_trimmed/      (entrada)
#        - results/<BIOPROJECT_ID>/multiqc_trimmed/ (salida)
#   3. Crea el directorio de salida si no existe.
#   4. Ejecuta MultiQC sobre el directorio qc_trimmed/.
#   5. Verifica que el archivo HTML final exista.
#
# Uso:
#   ./07_multiqc_trimmed.sh
#   ./07_multiqc_trimmed.sh --bioproject 482464
#   ./07_multiqc_trimmed.sh -b 482464
#
# Flags:
#   --bioproject/-b : ID de subcarpeta bajo results/.
#
# Dependencias:
#   - MultiQC instalado y disponible en el PATH.
#   - Haber corrido previamente FastQC sobre los FASTQ trimmed, de forma que
#     existan archivos de salida en results/<BIOPROJECT_ID>/qc_trimmed/.
#
# Estructura esperada (minima):
#   PROJECT_ROOT/
#     └─ results/
#          └─ <BIOPROJECT_ID>/
#               ├─ qc_trimmed/        # entrada: reportes de FastQC
#               └─ multiqc_trimmed/   # salida: reporte de MultiQC
# -----------------------------------------------------------------------------
set -euo pipefail
# set -e  : termina el script si algún comando devuelve error
# set -u  : error si se usa una variable no definida
# pipefail: hace que un pipe falle si falla alguno de sus comandos internos

# Detectar ruta del script y raíz del proyecto
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Raíz del proyecto: un nivel arriba del directorio del script
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# BIOPROJECT por defecto, puede ser modificado por flag
BIOPROJECT_ID="482464" 

# ---------- Parseo de flags ----------
# Permite definir el ID de bioproyecto:
#   --bioproject/-b
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
      echo "                  Entradas esperadas en: results/ID/qc_trimmed/"
      echo "                  Salida en: results/ID/multiqc_trimmed/"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

# Directorio base de resultados para este bioproyecto
RESULTS_DIR="${PROJECT_ROOT}/results"/${BIOPROJECT_ID}
# Directorio donde se esperan los reportes de FastQC sobre trimmed
QC_TRIM_DIR="${RESULTS_DIR}/qc_trimmed" 
# Directorio donde se guardarán los reportes de MultiQC sobre trimmed
MULTIQC_TRIM_DIR="${RESULTS_DIR}/multiqc_trimmed"
# Nombre del archivo de reporte final
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