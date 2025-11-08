#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

REF_DIR="${PROJECT_ROOT}/data/references"
BAM_DIR="${PROJECT_ROOT}/results/bam"
COUNTS_DIR="${PROJECT_ROOT}/results/counts"

#  Parámetros por defecto
#  Permite que el usuario pase el nombre del .gff si quiere.
GFF_BASENAME=""
# (si se deja vacío, intentaremos autodetectarlo)

#   05_counts.sh [--gff-name archivo.gff]
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gff)
      GFF_BASENAME="$2"
      shift 2
      ;;
    --help|-h)
      echo "Uso: $0 [--gff  ARCHIVO_GFF]"
      echo
      echo "Si no se pasa --gff, el script intenta detectar un único .gff en:"
      echo "  ${REF_DIR}"
      echo
      echo "Salida:"
      echo "  - BAMs leídos desde:  ${BAM_DIR}"
      echo "  - Counts por muestra: ${COUNTS_DIR}/<sample>.counts.txt"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

# Resolver anotación GFF
ANN_GFF=""

if [[ -n "${GFF_BASENAME}" ]]; then
  # El usuario dio un nombre explícito
  ANN_GFF="${REF_DIR}/${GFF_BASENAME}"
  if [[ ! -f "${ANN_GFF}" ]]; then
    echo "ERROR: No se encontró el GFF especificado: ${ANN_GFF}" >&2
    exit 1
  fi
else
  # Autodetectar un único .gff en data/reference
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
shopt -s nullglob
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
