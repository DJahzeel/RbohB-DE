import subprocess
from pathlib import Path

def submit_fasterqdump_job(job_name: str, script_path: Path, srr_list_path: Path, output_dir: Path, wait_for: bool = False):
    """
    Lanza un trabajo a qsub con el script y la lista de SRR especificados.
    
    Args:
        job_name (str): Nombre descriptivo para los logs.
        script_path (Path): Ruta al script .jdl que se ejecutará.
        srr_list_path (Path): Ruta al archivo .tsv (lista de SRRs) 
                                que se pasará como argumento $1.
    """
    
    # Verificación de que los archivos de entrada existan
    if not srr_list_path or not srr_list_path.exists():
        print(f"No se encontró archivo de SRR para '{job_name}' en {srr_list_path}, no se enviará el trabajo.")
        return
    
    if not script_path.exists():
        print(f"ERROR: No se encontró el script JDL: {script_path}")
        return
    
    if not output_dir.is_dir():
        print(f"ERROR: El directorio de salida no existe: {output_dir}")
        return

    print(f"Enviando trabajo '{job_name}' al clúster con 'qsub'...")
    
    # Construimos el comando (convirtiendo los Paths a strings)
    comando = [
        "qsub", 
        str(script_path), 
        str(srr_list_path),
        str(output_dir)
    ]

    if wait_for:
        comando.extend(["-sync", "y"])

    try:
        resultado = subprocess.run(comando, check=True, capture_output=True, text=True)
        
        print(f"Trabajo '{job_name}' enviado exitosamente.")
        print(f"Salida de qsub: {resultado.stdout.strip()}")
        
    except FileNotFoundError:
        print(f"ERROR: El comando 'qsub' no se encontró.")
        print("Asegúrate de estar en el servidor y que 'qsub' esté en tu PATH.")
    except subprocess.CalledProcessError as e:
        print(f"ERROR al enviar el trabajo '{job_name}':")
        print(f"Comando: {' '.join(comando)}")
        print(f"Error (stderr): {e.stderr.strip()}")
    except Exception as e:
        print(f"Ocurrió un error inesperado al lanzar el subprocess: {e}")