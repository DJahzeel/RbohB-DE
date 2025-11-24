#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 06_qc_trimmed.sh
# Descripción:
#   Realiza control de calidad (FastQC) sobre los archivos FASTQ ya recortados
#   (trimmed) producidos por el script de trimming (p.ej. 10_trimmed.sh).
#
# Flujo general:
#   1. Detecta la raíz del proyecto (PROJECT_ROOT) desde la ubicación del script.
#   2. Usa BIOPROJECT_ID para localizar:
#        - results/<BIOPROJECT_ID>/trimmed/       (entrada: FASTQ recortados)
#        - results/<BIOPROJECT_ID>/qc_trimmed/    (salida: reportes FastQC)
#   3. Busca todos los archivos *.fastq.gz en el directorio trimmed/.
#   4. Ejecuta FastQC en paralelo sobre todos esos archivos.
#
# Uso:
#   ./06_qc_trimmed.sh
#   ./06_qc_trimmed.sh --threads 8
#   ./06_qc_trimmed.sh --bioproject 482464
#
# Flags:
#   --threads/-t     : número de hilos para FastQC (ejecutados vía xargs -P)
#   --bioproject/-b  : ID de subcarpeta bajo results/
#
# Dependencias:
#   - FastQC instalado y disponible en el PATH.
#   - Haber corrido antes el script de trimming que genera:
#       results/<BIOPROJECT_ID>/trimmed/*.fastq.gz
#
# Estructura esperada:
#   PROJECT_ROOT/
#     ├─ results/
#     │    └─ <BIOPROJECT_ID>/
#     │          ├─ trimmed/      # entrada de este script
#     │          │     ├─ *_clean.fastq.gz
#     │          │     └─ *_clean_*.fastq.gz
#     │          └─ qc_trimmed/   # salida de este script
# -----------------------------------------------------------------------------
set -euo pipefail
# set -e  : termina el script si algún comando falla
# set -u  : error si se usa una variable no definida
# pipefail: si falla algún comando de un pipe, el pipe completo falla

# Detectar el directorio donde está este script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Subir un nivel para obtener la raíz del proyecto
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Valores por defecto. Parámetros modificables vía flags.
THREADS=4
BIOPROJECT_ID="482464"   # Default

# ---------- Parseo de flags ----------
# Permite modificar:
#   - número de hilos (--threads/-t)
#   - bioproyecto (--bioproject/-b)
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
      echo "  --bioproject/-b ID de subcarpeta bajo data/(default: 482464)"
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

# Directorio base de resultados para este bioproyecto
RESULTS_DIR="${PROJECT_ROOT}/results"/${BIOPROJECT_ID}
# Directorio donde se esperan los FASTQ recortados
TRIM_DIR="${RESULTS_DIR}/trimmed"
# Directorio donde se guardarán los reportes de FastQC sobre trimmed
QC_TRIM_DIR="${RESULTS_DIR}/qc_trimmed"

# Verificar que exista el directorio con FASTQ trimmed
if [[ ! -d "${TRIM_DIR}" ]]; then
  echo "ERROR: No existe el directorio con FASTQ trimmed: ${TRIM_DIR}" >&2
  echo "       Asegúrate de haber corrido primero trim.sh." >&2
  exit 1
fi

# # Crear carpeta de salida para resultados de FastQC (si no existe)
if [[ ! -d "${QC_TRIM_DIR}" ]]; then
  echo "[INFO] Creando directorio de salida para QC trimmed: ${QC_TRIM_DIR}"
  mkdir -p "${QC_TRIM_DIR}"
fi

# Construir lista de archivos FASTQ recortados (*.fastq.gz)
echo "[INFO] Buscando FASTQ trimmed en ${TRIM_DIR}"
mapfile -t FQS < <(find "${TRIM_DIR}" -type f -name "*.fastq.gz" | sort)
# Verificar que se hayan encontrado archivos
if [[ "${#FQS[@]}" -eq 0 ]]; then
  echo "ERROR: No se encontraron archivos *.fastq.gz en ${TRIM_DIR}" >&2
  exit 1
fi

echo "[INFO] Se encontraron ${#FQS[@]} archivos FASTQ trimmed para QC:"
printf '  - %s\n' "${FQS[@]}"

# Ejecutar FastQC en paralelo usando xargs:
#   -n 1 : un archivo por invocación
#   -P THREADS : número de procesos simultáneos
printf "%s\n" "${FQS[@]}" \
  | xargs -n 1 -P "${THREADS}" fastqc -o "${QC_TRIM_DIR}"

echo "QC de trimmed completado. Resultados en: ${QC_TRIM_DIR}"