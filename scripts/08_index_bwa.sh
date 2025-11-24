#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: 08_index_bwa.sh
# Descripción:
#   Construye (si no existen) los índices de BWA para un genoma de referencia
#   en formato FASTA. La referencia se busca dentro de:
#       PROJECT_ROOT/data/references/
#
#   Los índices se almacenan en:
#       PROJECT_ROOT/data/references/bwa/<PREFIJO>*
#
# Flujo general:
#   1. Detecta la raíz del proyecto (PROJECT_ROOT) a partir de la ubicación
#      del propio script.
#   2. Define un directorio base de referencias: data/references/.
#   3. Permite indicar:
#        - el nombre del archivo FASTA de referencia (--ref)
#        - el prefijo para los índices de BWA (--prefix)
#   4. Verifica que exista el archivo FASTA de referencia.
#   5. Crea el directorio de índices si no existe (data/references/bwa/).
#   6. Si ya existe el archivo <prefijo>.bwt, asume que los índices están
#      construidos y no los regenera.
#   7. Si no existen, llama a:
#         bwa index -p <PREFIJO> <REF_FA>
#
# Uso:
#   ./08_index_bwa.sh
#   ./08_index_bwa.sh --ref pvulgaris_v2.fa
#   ./08_index_bwa.sh --prefix pvulgaris_v2
#   ./08_index_bwa.sh --ref otro_genoma.fa --prefix otro_prefijo
#
# Flags:
#   --ref NOMBRE_FASTA   : nombre del archivo FASTA de referencia dentro de
#                          data/references/ (por defecto pvulgaris_v2.fa)
#   --prefix PREFIJO     : prefijo para los índices BWA bajo data/references/bwa/
#   --help/-h            : muestra este mensaje de ayuda y termina.
#
# Estructura esperada (minima):
#   PROJECT_ROOT/
#     └─ data/
#          └─ references/
#               ├─ pvulgaris_v2.fa           # referencia por defecto
#               └─ bwa/                      # índices generados aquí
#
# Dependencias:
#   - bwa instalado y disponible en el PATH.
# -----------------------------------------------------------------------------
set -euo pipefail
# set -e  : termina el script si cualquier comando devuelve código de error
# set -u  : error si se usa una variable no definida
# pipefail: si un comando en un pipe falla, el pipe completo falla

# Detectar directorio donde está el script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Raíz del proyecto: un nivel arriba del directorio del script
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Directorio base de referencia
REF_DIR="${PROJECT_ROOT}/data/references"

# Valores por defecto, modificables vía flags
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
      echo "  REF_FA         = \${PROJECT_ROOT}/data/references/\${ref-name}"
      echo "  BWA_INDEX_PREF = \${PROJECT_ROOT}/data/references/bwa/\${prefix-name}"
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

# Verificar que exista el archivo FASTA de referencia
if [[ ! -f "${REF_FA}" ]]; then
  echo "ERROR: No se encontró el genoma de referencia en: ${REF_FA}" >&2
  exit 1
fi

# Crear carpeta de índices si no existe
if [[ ! -d "${INDEX_DIR}" ]]; then
  echo "[INFO] Creando directorio para índices BWA: ${INDEX_DIR}"
  mkdir -p "${INDEX_DIR}"
fi

# Información al usuario
echo "[INFO] Usando referencia:     ${REF_FA}"
echo "[INFO] Prefijo de índices:    ${BWA_INDEX_PREFIX}"

# Construir índices solo si no existen
if [[ -f "${BWA_INDEX_PREFIX}.bwt" ]]; then
  echo "[SKIP] Índices BWA ya existen para este prefijo (${BWA_INDEX_PREFIX})."
else
  echo "[INFO] Construyendo índices BWA..."
  bwa index -p "${BWA_INDEX_PREFIX}" "${REF_FA}"
  echo "[OK] Índices BWA generados."
fi