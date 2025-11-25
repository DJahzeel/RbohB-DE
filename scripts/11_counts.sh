#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 10_counts.sh 
# Descripción:
#   Genera archivos de conteos por característica usando
#   coverageBed (bedtools), a partir de:
#     - BAMs alineados (uno por muestra)
#     - Un archivo de anotación en formato GFF
#
#   Para cada archivo BAM en:
#       results/<BIOPROJECT_ID>/bam/*.bam
#   se ejecuta:
#       coverageBed -a <anotación.gff> -b <sample.bam>
#   y se guarda la salida en:
#       results/<BIOPROJECT_ID>/counts/<sample>.counts.txt
#
# Flujo general:
#   1. Detecta la raíz del proyecto (PROJECT_ROOT) desde la ubicación del script.
#   2. Localiza:
#        - Directorio de referencias: data/references/
#        - Directorio con BAMs:       results/<BIOPROJECT_ID>/bam/
#        - Directorio de salida:      results/<BIOPROJECT_ID>/counts/
#   3. Determina el archivo GFF:
#        - Si se pasa --gff ARCHIVO_GFF → usa ese.
#        - Si no, intenta autodetectar un único .gff en data/references/.
#   4. Recorre todos los BAM (*.bam) en BAM_DIR y para cada uno:
#        - Obtiene el nombre de muestra a partir del basename del BAM.
#        - Ejecuta coverageBed con ANN_GFF como -a y el BAM como -b.
#        - Escribe la salida en <sample>.counts.txt.
#
# Uso:
#   ./10_counts.sh
#   ./10_counts.sh --bioproject 482464
#   ./10_counts.sh --gff anotacion.gff
#
# Flags:
#   --bioproject/-b ID  : ID de subcarpeta bajo results/ y data/.
#   --gff ARCHIVO_GFF   : nombre de archivo GFF dentro de data/references/
#                         (si no se pasa, se intenta autodetectar un único .gff).
#
# Estructura esperada (minima):
#   PROJECT_ROOT/
#     ├─ data/
#     │    └─ references/
#     │         ├─ *.gff
#     └─ results/
#          └─ <BIOPROJECT_ID>/
#               ├─ bam/           # entrada: archivos .bam
#               └─ counts/        # salida: conteos por muestra
#
# Dependencias:
#   - bedtools instalado y coverageBed disponible en el PATH.
#   - BAMs generados previamente (p.ej. con el script de alineamiento).
# -----------------------------------------------------------------------------
set -euo pipefail
# set -e  : termina el script si algún comando devuelve error
# set -u  : error si se usa una variable no definida
# pipefail: si un comando en un pipe falla, el pipe completo falla

# Determinar rutas
# Obtener la ruta del directorio donde está este script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# La raíz del proyecto es el padre del directorio del script
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
#| Directorio de referencias
REF_DIR="${PROJECT_ROOT}/data/references"

#  Parámetros por defecto, pueden ser sobrescritos por flags
BIOPROJECT_ID="482464" 
#  Permite que el usuario pase el nombre del .gff si quiere.
GFF_BASENAME=""
# (si se deja vacío, intentaremos autodetectarlo)

# Parseo de flags
# [--bioproject ID] [--gff  ARCHIVO_GFF]"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bioproject|-b)
      BIOPROJECT_ID="$2"
      shift 2
      ;;
    --gff)
      GFF_BASENAME="$2"
      shift 2
      ;;
    --help|-h)
      echo "Uso: $0 [--bioproject ID] [--gff  ARCHIVO_GFF]"
      echo
      echo "Si no se pasa --gff, el script intenta detectar un único .gff en:"
      echo "  ${REF_DIR}"
      echo
      echo "Salida:"
      echo "  - BAMs leídos desde:  results/ID/bam"
      echo "  - Counts por muestra: results/ID/counts/<sample>.counts.txt"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done
# Rutas de entrada/salida
BAM_DIR="${PROJECT_ROOT}/results/${BIOPROJECT_ID}/bam"
COUNTS_DIR="${PROJECT_ROOT}/results/${BIOPROJECT_ID}/counts"

# Resolver anotación GFF, ya sea pasada por el usuario o autodetectada
ANN_GFF=""

if [[ -n "${GFF_BASENAME}" ]]; then
  # El usuario dio un nombre explícito, usarlo
  ANN_GFF="${REF_DIR}/${GFF_BASENAME}"
  # Verificar que el archivo exista
  if [[ ! -f "${ANN_GFF}" ]]; then
    echo "ERROR: No se encontró el GFF especificado: ${ANN_GFF}" >&2
    exit 1
  fi
else
  # el usuario no dio un GFF, intentar autodetectar uno único en REF_DIR
  mapfile -t GFF_FILES < <(find "${REF_DIR}" -maxdepth 1 -type f -name "*.gff" | sort)
  if [[ "${#GFF_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No se encontró ningún .gff en ${REF_DIR}" >&2
    echo "       Pasa uno explícito con: 05_counts.sh --gff-name archivo.gff" >&2
    exit 1
  elif [[ "${#GFF_FILES[@]}" -gt 1 ]]; then
    echo "ERROR: Se encontraron múltiples .gff en ${REF_DIR}:" >&2
    printf '  - %s\n' "${GFF_FILES[@]}" >&2
    echo "       Especifica cuál usar con: --gff  archivo.gff" >&2
    exit 1
  else
    ANN_GFF="${GFF_FILES[0]}"
  fi
fi

echo "[DEBUG] Ruta GFF/Genoma que se usará: ${ANN_GFF}"

if [[ ! -d "${BAM_DIR}" ]]; then
  echo "ERROR: No existe el directorio de BAMs: ${BAM_DIR}" >&2
  echo "       Asegúrate de haber corrido 03_align.sh." >&2
  exit 1
fi

# Crear carpeta de salida si no existe
if [[ ! -d "${COUNTS_DIR}" ]]; then
  echo "[INFO] Creando directorio de counts: ${COUNTS_DIR}"
  mkdir -p "${COUNTS_DIR}"
fi

# Recorrer cada BAM y correr coverageBed
# Usar nullglob para que el patrón *.bam no falle si no hay archivos,llene el array con cero elementos 
# para manejar el caso de no encontrar BAMs.
shopt -s nullglob
# Obtener lista de BAMs
BAMS=( "${BAM_DIR}"/*.bam )
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "ERROR: No se encontraron archivos .bam en ${BAM_DIR}" >&2
  exit 1
fi
for bam in "${BAMS[@]}"; do
  sample="$(basename "${bam%.bam}")"
  out="${COUNTS_DIR}/${sample}.counts.txt"

  echo "[INFO] Procesando BAM: ${bam}"
  echo "       -> ${out}"

  # coverageBed (bedtools) con anotación en -a y BAM en -b
  coverageBed -a "${ANN_GFF}" -b "${bam}" > "${out}"
done

echo "Cuantificación completada. Archivos en: ${COUNTS_DIR}"