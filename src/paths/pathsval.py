from pathlib import Path

from pathlib import Path
import sys

def project_root() -> Path:
    """
    Raíz del proyecto:
    1) Parte del entrypoint (__main__.__file__) si existe; si no, del archivo actual.
    2) Sube buscando 'pyproject.toml' o '.git' y regresa ahí.
    3) Si no hay marcadores, si ve 'scripts' en el camino, regresa su padre.
    4) Fallback: el padre del directorio base.
    """
    start = Path(getattr(sys.modules.get("__main__"), "__file__", __file__)).resolve()
    base = start if start.is_dir() else start.parent

    # 2) Marcadores de proyecto
    for p in (base,) + tuple(base.parents):
        if (p / "pyproject.toml").exists() or (p / ".git").exists():
            return p

    # 3) Padre de 'scripts' si existe en la cadena
    for p in (base,) + tuple(base.parents):
        if p.name == "scripts":
            return p.parent

    # 4) Fallback
    return base.parent

def dboutput_directory(record: dict, root: Path) -> Path:
    """Get the output directory for a given record."""
    try:
        bioproject = record["bioproject"]
        database = record["database"]
    except KeyError as e:
        raise KeyError(f"Falta clave en record: {e}") from e

    outdir = root / "data" / bioproject / f"{database}_files"

    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        raise PermissionError(f"No hay permisos para crear {outdir}")
    except OSError as e:
        raise OSError(f"Error del sistema de archivos creando {outdir}: {e}") from e

    return outdir

def summary_directory(root: Path) -> Path:
    sum_dir = root / "data" / "" 

    try:
        sum_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError as permission_error:
        raise PermissionError(f"No hay permisos para crear {sum_dir}") from permission_error
    except OSError as e:
        raise OSError(f"Error del sistema de archivos creando {sum_dir}: {e}") from e

    return sum_dir

def references_directory(root: Path) -> Path:
    ref_dir = root / "data" / "references" 

    try:
        ref_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError as permission_error:
        raise PermissionError(f"No hay permisos para crear {ref_dir}") from permission_error
    except OSError as e:
        raise OSError(f"Error del sistema de archivos creando {ref_dir}: {e}") from e

    return ref_dir

def runends_directory(root: Path, bioproject: str) -> Path:
    re_dir = root / "results" / bioproject / "runends" 

    try:
        re_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError as permission_error:
        raise PermissionError(f"No hay permisos para crear {re_dir}") from permission_error
    except OSError as e:
        raise OSError(f"Error del sistema de archivos creando {re_dir}: {e}") from e
    
    return re_dir