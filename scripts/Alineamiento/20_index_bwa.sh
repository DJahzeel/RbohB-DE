#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

# Directorio base de referencia
REF_DIR="${PROJECT_ROOT}/data/references"

# Valores por defecto
REF_BASENAME="pvulgaris_v2.fa"
PREFIX_BASENAME="pvulgaris_v2"

# Parseo de flags
# Uso:
#   index_bwa.sh [--ref NOMBRE_FA] [--prefix PREFIJO]
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref)
      REF_BASENAME="$2"
      shift 2
      ;;
    --prefix)
      PREFIX_BASENAME="$2"
      shift 2
      ;;
    --help|-h)
      echo "Uso: $0 [--ref NOMBRE_FASTA] [--prefix PREFIJO]"
      echo
      echo "Por defecto:"
      echo "  --ref    ${REF_BASENAME}"
      echo "  --prefix ${PREFIX_BASENAME}"
      echo
      echo "Esto construye:"
      echo "  REF_FA         = \${PROJECT_ROOT}/data/reference/\${ref-name}"
      echo "  BWA_INDEX_PREF = \${PROJECT_ROOT}/data/reference/bwa/\${prefix-name}"
      exit 0
      ;;
    *)
      echo "ERROR: Opción desconocida: $1" >&2
      echo "Usa --help para ver las opciones." >&2
      exit 1
      ;;
  esac
done

#  Construir rutas 
REF_FA="${REF_DIR}/${REF_BASENAME}"
BWA_INDEX_PREFIX="${REF_DIR}/bwa/${PREFIX_BASENAME}"
INDEX_DIR="$(dirname "${BWA_INDEX_PREFIX}")"

if [[ ! -f "${REF_FA}" ]]; then
  echo "ERROR: No se encontró el genoma de referencia en: ${REF_FA}" >&2
  exit 1
fi

# Crear carpeta de índices si no existe
if [[ ! -d "${INDEX_DIR}" ]]; then
  echo "[INFO] Creando directorio para índices BWA: ${INDEX_DIR}"
  mkdir -p "${INDEX_DIR}"
fi

echo "[INFO] Usando referencia:     ${REF_FA}"
echo "[INFO] Prefijo de índices:    ${BWA_INDEX_PREFIX}"

# 7) Construir índices solo si no existen
if [[ -f "${BWA_INDEX_PREFIX}.bwt" ]]; then
  echo "[SKIP] Índices BWA ya existen para este prefijo (${BWA_INDEX_PREFIX})."
else
  echo "[INFO] Construyendo índices BWA..."
  bwa index -p "${BWA_INDEX_PREFIX}" "${REF_FA}"
  echo "[OK] Índices BWA generados."
fi

