#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 04_multiqc_raw.sh
# Descripción:
#   Genera un reporte integrado de calidad usando MultiQC a partir de los
#   reportes individuales de FastQC (RAW), previamente generados en:
#       results/<BIOPROJECT_ID>/qc/
#
# Flujo general:
#   1. Detecta la raíz del proyecto a partir de la ubicación del script.
#   2. Permite elegir el BIOPROJECT_ID mediante flag (--bioproject).
#   3. Define rutas de entrada (FastQC) y salida (MultiQC).
#   4. Ejecuta MultiQC para consolidar todos los reportes.
#   5. Verifica que el HTML final se haya generado correctamente.
#
# Uso:
#   ./04_multiqc_raw.sh
#   ./04_multiqc_raw.sh --bioproject 482464
#   ./04_multiqc_raw.sh -b 482464
#
# Estructura esperada:
#   PROJECT_ROOT/
#     ├─ results/
#     │     └─ <BIOPROJECT_ID>/
#     │           ├─ qc/              # Salida de FastQC
#     │           └─ multiqc/         # Salida de este script
#
# Salida principal:
#     results/<BIOPROJECT_ID>/multiqc/multiqc_report.html
#
# Dependencias:
#   - MultiQC instalado y disponible en el PATH.
#   - Directorio results/<BIOPROJECT_ID>/qc/ con resultados de FastQC.
# -----------------------------------------------------------------------------
set -euo pipefail
# set -e  : aborta si un comando falla
# set -u  : error si se usa una variable no declarada
# set -o pipefail : si falla un comando dentro de un pipe, falla el pipe completo

# Detectar ruta del script y raíz del proyecto
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# BIOPROJECT por defecto
BIOPROJECT_ID="482464"  # Puede ser modificado por flag

# ---------- Parseo de flags ----------
# Flags disponibles:
#   --bioproject/-b  : especifica el ID de bioproyecto
#   --help/-h        : muestra ayuda
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
      echo "                  Entradas esperadas en: results/ID/qc/"
      echo "                  Salida: results/ID/multiqc/multiqc_report.html"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

# Directorios fijos basados en BIOPROJECT_ID
RESULTS_DIR="${PROJECT_ROOT}/results"/${BIOPROJECT_ID}
QC_DIR="${RESULTS_DIR}/qc"               # Entrada: reportes FastQC RAW
MULTIQC_DIR="${RESULTS_DIR}/multiqc"     # Salida: reportes MultiQC
REPORT_NAME="multiqc_report.html"

# Información al usuario
echo "[INFO] MultiQC (RAW)"
echo "  - Entrada QC:   ${QC_DIR}"
echo "  - Salida HTML:  ${MULTIQC_DIR}/${REPORT_NAME}"

# Crear solo la carpeta de salida, si no existe
mkdir -p "${MULTIQC_DIR}"

# Correr MultiQC leyendo de QC_DIR y escribiendo en MULTIQC_DIR
multiqc "${QC_DIR}" -o "${MULTIQC_DIR}" -n "${REPORT_NAME}"

# Verificar que el reporte haya sido creado
FULL_REPORT_PATH="${MULTIQC_DIR}/${REPORT_NAME}"
if [[ -f "${FULL_REPORT_PATH}" ]]; then
    echo "¡Reporte de MultiQC generado con éxito!"
    echo "Reporte disponible en: ${FULL_REPORT_PATH}"
else
    echo "ERROR: Falló la generación del reporte de MultiQC."
    exit 1
fi