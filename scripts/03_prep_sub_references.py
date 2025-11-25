from download_data.get_downloads_info import (
    find_references,
    save_ftp_links
)
from paths.pathsval import (
    referencelinks_directory,
    refseq_directory,
    project_root,
)
from utils.job_submission import submit_qsub_job

def main():
    base_project_root = project_root()
    reference_genomes_meta = base_project_root / "data" / "references"

    found_organisms = find_references(reference_genomes_meta)
    if not found_organisms:
        print("No se encontraron genomas de referencia en el directorio especificado.")
        return
    
    print("Organismos encontrados en los archivos de referencia , seleccione un numero o rango de los que desea descargar:")
    only_organism = []
    for i, org in enumerate(found_organisms):
        only_organism.append(org['Organism'])
        print(f"{i+1}. {org['Organism']}")
    selected_indices = input("Ingrese los números separados por comas (por ejemplo, 1,3,5) o un rango (por ejemplo, 1-3): ")
    if '-' in selected_indices:
        start, end = map(int, selected_indices.split('-'))
        selected_indices = list(range(start - 1, end))
    else:
        selected_indices = [int(i) - 1 for i in selected_indices.split(',')]
    
    processing_prefix = f"ref_download_{len(selected_indices)}_items"
    link_dir = referencelinks_directory(base_project_root)
    to_download = save_ftp_links(selected_indices, found_organisms, link_dir, only_organism)

    references_script = base_project_root / "src" / "references_bash" / "refseq_download.jdl"
    refseq_path = refseq_directory(base_project_root)
    submit_qsub_job(f"{processing_prefix}_reference_download", references_script, to_download, refseq_path, wait_for=True)
    print(f"Descarga de genomas de referencia para {processing_prefix} terminada.")
    
if __name__ == "__main__":
    main()  