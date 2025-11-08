#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

REF_DIR="${PROJECT_ROOT}/data/reference"
BAM_DIR="${PROJECT_ROOT}/results/bam"
COUNTS_DIR="${PROJECT_ROOT}/results/counts"

# ============================
# Defaults estilo artículo
# ============================
# htseq-count defaults (ajustados a tu pipeline):
#   -f bam        (tú usas BAM)
#   -r pos        (tus BAM están ordenados por posición)
#   -s yes        (htseq por defecto asume librerías stranded)
#   -t exon       (se cuentan exones)
#   -i gene_id    (atributo gene_id del GFF)
#
HTSEQ_STRANDED="yes"    # -s yes | no | reverse
HTSEQ_FEATURE="exon"    # -t exon (puede cambiar a gene, CDS, etc.)
HTSEQ_ATTR="gene_id"    # -i gene_id (puede cambiar a ID, locus_tag, ...)
HTSEQ_ORDER="pos"       # -r pos (porque tus BAM están pos-sorted)

# El usuario puede elegir el .gff concreto dentro de data/reference
GFF_BASENAME=""

# ============================
# Parseo de flags
# ============================
# Uso:
#   40_counts.sh [--gff ARCHIVO_GFF] \
#                [--stranded yes|no|reverse] \
#                [--feature exon|gene|CDS] \
#                [--attr gene_id|ID|...] \
#                [--order pos|name]
#
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gff)
      GFF_BASENAME="$2"
      shift 2
      ;;
    --stranded)
      HTSEQ_STRANDED="$2"
      shift 2
      ;;
    --feature)
      HTSEQ_FEATURE="$2"
      shift 2
      ;;
    --attr)
      HTSEQ_ATTR="$2"
      shift 2
      ;;
    --order)
      HTSEQ_ORDER="$2"
      shift 2
      ;;
    --help|-h)
      echo "Uso: $0 [opciones]"
      echo
      echo "Opciones (solo estos parámetros de htseq-count son configurables):"
      echo "  --gff ARCHIVO_GFF      Nombre del .gff dentro de data/reference."
      echo "                         Si se omite, se autodetecta un único .gff."
      echo
      echo "  --stranded VAL         Valor para -s (default: yes)."
      echo "                         Posibles: yes | no | reverse"
      echo
      echo "  --feature TIPO         Valor para -t (default: exon)."
      echo "                         Ejemplos: exon, gene, CDS"
      echo
      echo "  --attr NOMBRE          Valor para -i (default: gene_id)."
      echo "                         Ejemplos: gene_id, ID, locus_tag"
      echo
      echo "  --order VAL            Valor para -r (default: pos)."
      echo "                         Posibles: pos | name"
      echo
      echo "Ejemplos:"
      echo "  $0"
      echo "  $0 --gff GCA_000499845.2_P._vulgaris_v2.0_genomic.gff"
      echo "  $0 --stranded no --feature gene --attr ID"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

# ============================
# Resolver anotación GFF
# ============================
ANN_GFF=""

if [[ -n "${GFF_BASENAME}" ]]; then
  ANN_GFF="${REF_DIR}/${GFF_BASENAME}"
  if [[ ! -f "${ANN_GFF}" ]]; then
    echo "ERROR: No se encontró el GFF especificado: ${ANN_GFF}" >&2
    exit 1
  fi
else
  mapfile -t GFF_FILES < <(find "${REF_DIR}" -maxdepth 1 -type f -name "*.gff" | sort)
  if [[ "${#GFF_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No se encontró ningún .gff en ${REF_DIR}" >&2
    echo "       Pasa uno explícito con: $0 --gff ARCHIVO.gff" >&2
    exit 1
  elif [[ "${#GFF_FILES[@]}" -gt 1 ]]; then
    echo "ERROR: Se encontraron múltiples .gff en ${REF_DIR}:" >&2
    printf '  - %s\n' "${GFF_FILES[@]}" >&2
    echo "       Especifica cuál usar con: --gff ARCHIVO.gff" >&2
    exit 1
  else
    ANN_GFF="${GFF_FILES[0]}"
  fi
fi

echo "[INFO] 40_counts.sh (htseq-count)"
echo "  - PROJECT_ROOT:   ${PROJECT_ROOT}"
echo "  - BAM_DIR:        ${BAM_DIR}"
echo "  - COUNTS_DIR:     ${COUNTS_DIR}"
echo "  - GFF:            ${ANN_GFF}"
echo "  - HTSEQ_STRANDED: ${HTSEQ_STRANDED}"
echo "  - HTSEQ_FEATURE:  ${HTSEQ_FEATURE}"
echo "  - HTSEQ_ATTR:     ${HTSEQ_ATTR}"
echo "  - HTSEQ_ORDER:    ${HTSEQ_ORDER}"

if [[ ! -d "${BAM_DIR}" ]]; then
  echo "ERROR: No existe el directorio de BAMs: ${BAM_DIR}" >&2
  echo "       Asegúrate de haber corrido el script de alineamiento." >&2
  exit 1
fi

# Verificar htseq-count
if ! command -v htseq-count >/dev/null 2>&1; then
  echo "ERROR: No encuentro htseq-count en el PATH." >&2
  echo "       Activa el entorno de conda donde esté instalado." >&2
  exit 1
fi

# Crear carpeta de salida si no existe
if [[ ! -d "${COUNTS_DIR}" ]]; then
  echo "[INFO] Creando directorio de counts: ${COUNTS_DIR}"
  mkdir -p "${COUNTS_DIR}"
fi

# ============================
# Recorrer cada BAM y correr htseq-count
# ============================
shopt -s nullglob
BAMS=( "${BAM_DIR}"/*.bam )
shopt -u nullglob

if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "ERROR: No se encontraron archivos .bam en ${BAM_DIR}" >&2
  exit 1
fi

for bam in "${BAMS[@]}"; do
  sample="$(basename "${bam%.bam}")"
  out="${COUNTS_DIR}/${sample}.counts.txt"

  echo "[INFO] Procesando BAM: ${bam}"
  echo "       -> ${out}"

  htseq-count \
    -f bam \
    -r "${HTSEQ_ORDER}" \
    -s "${HTSEQ_STRANDED}" \
    -t "${HTSEQ_FEATURE}" \
    -i "${HTSEQ_ATTR}" \
    "${bam}" "${ANN_GFF}" \
    > "${out}"
done

echo "[OK] Cuantificación completada. Archivos en: ${COUNTS_DIR}"
