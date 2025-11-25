#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 03_qc_raw.sh
# Descripción:
#   Ejecuta FastQC sobre los archivos FASTQ crudos (sin recorte) definidos en
#   un archivo samples.tsv, para evaluar la calidad de las lecturas antes de
#   cualquier procesamiento posterior.
#
# Flujo general:
#   1. Localiza la raíz del proyecto (PROJECT_ROOT) a partir de la ubicación
#      del script.
#   2. Define el archivo de muestras (samples.tsv) en la raíz del proyecto.
#   3. Permite modificar el número de hilos y el ID de bioproyecto mediante
#      flags (--threads, --bioproject).
#   4. Construye la ruta base donde se esperan los FASTQ:
#        data/<BIOPROJECT_ID>/fastq/
#   5. Lee samples.tsv (formato tabulado: sample  R1  R2) y verifica que
#      existan los archivos FASTQ correspondientes.
#   6. Llama a FastQC con todos los FASTQ válidos y guarda los resultados
#      en results/<BIOPROJECT_ID>/qc/.
#
# Formato esperado de samples.tsv (tabulado):
#   sample_id <TAB> archivo_R1.fastq.gz <TAB> [archivo_R2.fastq.gz]
#   - Puede contener encabezados "sample" o "sample_id".
#   - Puede contener líneas de comentario que comiencen con "#".
#   - El campo R2 puede estar vacío para datos single-end.
#
# Uso:
#   ./03_qc_raw.sh
#   ./03_qc_raw.sh --threads 8
#   ./03_qc_raw.sh --bioproject 482464
#   ./03_qc_raw.sh -t 8 -b 482464
#
# Requisitos:
#   - FastQC instalado y disponible en el PATH.
#   - Estructura de directorios:
#       PROJECT_ROOT/
#         ├─ samples.tsv
#         ├─ data/
#         │    └─ <BIOPROJECT_ID>/
#         │         └─ fastq/
#         │              ├─ *.fastq.gz
#         └─ results/
#
# Salida:
#   - Archivos de reporte FastQC (*.html, *.zip) en:
#       results/<BIOPROJECT_ID>/qc/
# -----------------------------------------------------------------------------
set -euo pipefail
# set -e  : termina el script si algún comando devuelve estado distinto de 0.
# set -u  : error si se usa una variable no definida.
# set -o pipefail : si algún comando en un pipe falla, el pipe completo falla.

# Determina la ruta absoluta del directorio donde se encuentra este script.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# A partir del directorio del script, sube un nivel para obtener la raíz del proyecto.
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Archivo de metadatos con la lista de muestras y nombres de FASTQ.
SAMPLES_TSV="${PROJECT_ROOT}/samples.tsv"

# Valores por defecto. Parametros modificables via flags.
THREADS=4
BIOPROJECT_ID="482464"   # Default

# ---------- Parseo de flags ----------
# Se procesan parámetros de línea de comando:
#   --threads/-t   : número de hilos para FastQC
#   --bioproject/-b: ID de subcarpeta bajo data/
#   --help/-h      : muestra ayuda y termina
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
      # Cualquier argumento no reconocido se considera un error.
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

# Directorio donde se esperan los FASTQ de este bioproyecto.
FASTQ_BASE="${PROJECT_ROOT}/data/${BIOPROJECT_ID}/fastq"
# Directorio base de resultados para este bioproyecto.
RESULTS_DIR="${PROJECT_ROOT}/results"/${BIOPROJECT_ID}
# Subcarpeta específica para los reportes de calidad (FastQC).
QC_DIR="${RESULTS_DIR}/qc"

# Mensaje informativo con rutas y parámetros efectivos.
echo "[INFO] 00_qc_raw.sh"
echo "  - PROJECT_ROOT: ${PROJECT_ROOT}"
echo "  - SAMPLES_TSV:  ${SAMPLES_TSV}"
echo "  - FASTQ_BASE:   ${FASTQ_BASE}"
echo "  - QC_DIR:       ${QC_DIR}"
echo "  - THREADS:      ${THREADS}"

# Verifica que exista el archivo samples.tsv.
if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: No encuentro samples.tsv en ${SAMPLES_TSV}" >&2
  exit 1
fi

# Crea directorios de salida solo si no existen.
mkdir -p "${RESULTS_DIR}" "${QC_DIR}"

# ---------- Recolección de archivos FASTQ ----------
# Lee samples.tsv y construye una lista de archivos FASTQ a procesar.
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
  # Si hay un archivo R2 (datos paired-end), también se valida y se añade.
  if [[ -n "${r2:-}" ]]; then
    fq2="${FASTQ_BASE}/${r2}"
    if [[ ! -f "${fq2}" ]]; then
      echo "ERROR: No encuentro FASTQ: ${fq2}" >&2
      exit 1
    fi
    FQS+=( "${fq2}" )
  fi
done < "${SAMPLES_TSV}"

# Si no se encontró ningún FASTQ válido, se aborta.
if [[ "${#FQS[@]}" -eq 0 ]]; then
  echo "ERROR: No se encontraron FASTQ a partir de samples.tsv" >&2
  exit 1
fi

# Lista de archivos FASTQ que se enviarán a FastQC.
echo "[INFO] FASTQ a analizar:"
printf '  - %s\n' "${FQS[@]}"

# Llamada a FastQC con todos los archivos detectados.
fastqc -t "${THREADS}" -o "${QC_DIR}" "${FQS[@]}"

echo "[OK] QC RAW completado en: ${QC_DIR}"