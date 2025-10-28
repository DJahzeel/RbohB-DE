import argparse as ap
import time
from pathlib import Path
from paths.pathsval import(
    project_root,
    dboutput_directory,
    references_directory,
)
from exploration.explore_project import (
    get_projects_id,
    get_projects_summary,
    save_project_summaries,
)

from exploration.explore_databases import (
    fetch_bioproject_links,
    geo_links,
    gse_gsm,
    sra_links,
)

from download_data.get_downloads_info import (
    get_assemblies_data,
    get_reference,
    run_end,
)

parser = ap.ArgumentParser(
    description=("Herramienta de línea de comandos para recuperar información asociada a "
            "uno o varios BioProjects de NCBI y preparar los metadatos necesarios "
            "para la descarga de datos de secuenciación.")
)
parser.add_argument("-p", "--projects",
                    nargs="+",
                    required=True,
                    help="Lista de accesiones de BioProject (e.g., PRJNA123456 PRJNA654321)")
parser.add_argument("-m", "--mail",
                    required=True,
                    help="Dirección de correo electrónico para Entrez de NCBI")
parser.add_argument("-d", "--database",
                    nargs="+",
                    default=["GDS", "SRA"],
                    help="Bases de datos a explorar (GEO, SRA)")
parser.add_argument("-o", "--organism",
                    required=False,
                    nargs="+", 
                    help="Nombre del organismo (e.g., Homo sapiens)")

args = parser.parse_args()

def main():
    print("Proyectos:", args.projects)
    print("Correo:", args.mail)
    print("Bases de datos:", args.database)
    print("Organismo:", args.organism)

    base_project_path = project_root()
    data_dir = base_project_path / "data"

    try:
        Path(data_dir).mkdir(exist_ok=True)
    except PermissionError:
        raise PermissionError(f"No hay permisos para crear {data_dir}")
    except OSError as e:
        raise OSError(f"Error del sistema de archivos creando {data_dir}: {e}") from e

    projects_id_dic = get_projects_id(args.projects, args.mail)
    if not projects_id_dic:
        raise RuntimeError("No se encontraron UIDs para los proyectos especificados.")
    
    project_uids = list(projects_id_dic.values())
    summaries = get_projects_summary(project_uids, args.mail)
    if not summaries:
        raise RuntimeError("No se encontraron resúmenes para los proyectos especificados, no se ha guardado ningún archivo.")
    save_project_summaries(summaries, base_project_path) if summaries else []

    links_to_db = fetch_bioproject_links(project_uids, args.database, args.mail)
    if not links_to_db:
        raise RuntimeError("No se encontraron enlaces a bases de datos para los proyectos especificados.")
    
    for record in links_to_db:
        outputdb = dboutput_directory(record, base_project_path)
        if record['database'] == 'gds':
            geo = geo_links(record, args.mail, outputdb)
            print(f"Guardado en {geo} # Añadir una pausa para evitar sobrecargar el servidor de NCBI")
            time.sleep(0.4)
            gsm_sample = gse_gsm(record, args.mail, outputdb)
            print(f"Guardado en {gsm_sample} # Añadir una pausa para evitar sobrecargar el servidor de NCBI")
            time.sleep(0.4)
        elif record['database'] == 'sra':
            table,sra_out = sra_links(args.mail, record, outputdb)
            print(f"Guardado en {sra_out} # Añadir una pausa para evitar sobrecargar el servidor de NCBI")
            time.sleep(0.4)
            single, paired= run_end(table, outputdb)
            print(f"Información de runs guardada en: {single}, {paired}")
        else:
            print(f"Base de datos no soportada: {record['database']}")
    
    if args.organism:
        for arg in args.organism:
            outputref = references_directory(base_project_path)
            assembly_data = get_assemblies_data(arg, args.mail)
            time.sleep(0.4)
            if assembly_data:
                referencegenome = get_reference(assembly_data, outputref, arg)
                print(f"Referencia guardada en: {referencegenome}")
            else:
                print(f"No se encontraron datos de ensamblaje para el organismo: {arg}")
 
if __name__ == "__main__":
    main()

