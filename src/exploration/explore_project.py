from Bio import Entrez
from urllib.error import HTTPError, URLError
import os
from pathlib import Path

def get_projects_id(accesion_ids,mail):
    try:

        found_ids = {}
        Entrez.email = mail
        search_str=" OR ".join(accesion_ids)

        handle = Entrez.esearch(db="bioproject", term=search_str)
        search_results = Entrez.read(handle)
        handle.close()

        if search_results['Count'] == '0':
            print("No project found for the given accession IDs.")
            return None
        
        elif search_results['Count'] >= '1':
            project_uid = search_results['IdList']
            uid_str = ",".join(project_uid)
            handle = Entrez.esummary(db="bioproject", id=uid_str)
            summaries = Entrez.read(handle)
            handle.close()

        for summary in summaries['DocumentSummarySet']['DocumentSummary']:
            print(f"Project ID: {summary['Project_Id']}, Accesion ID: {summary['Project_Acc']}")
            found_ids[summary['Project_Acc']] = summary['Project_Id']
        
        for accesion in accesion_ids:
            if accesion not in found_ids:
                print(f"No project found for accession ID: {accesion}")
        return found_ids
    
    except (HTTPError, URLError) as connection_error:
        print(f"Connection error: {connection_error}")
        return None
    except Exception as error:
        print(f"An error occurred: {error}")
        return None

def get_projects_summary(project_uid, mail):
    try:

        Entrez.email = mail
        uid_str = ",".join(project_uid)

        handle = Entrez.esummary(db="bioproject", id=uid_str)
        records = Entrez.read(handle)
        handle.close()
        return records['DocumentSummarySet']['DocumentSummary']
    
    except (HTTPError, URLError) as connection_error:
        raise ConnectionError(f"Connection error: {connection_error}") from connection_error
    except Exception as error:
        raise Exception(f"An error occurred: {error}") from error

def save_project_summaries(records, root: Path):
        
        if not records:
            print("No records to save.")
            return
        
        for record in records:
            project_uid = record.get('Project_Id')

            if not project_uid:
                print("No Project ID found in the record, skipping.")
                continue
            
            sum_dir = root/"data"/ project_uid
            try:
                sum_dir.mkdir(parents=True, exist_ok=True)
            except PermissionError:
                raise PermissionError(f"No permissions to create directory {sum_dir}, skipping.")
            except OSError as e:
                raise OSError(f"Error creating directory {sum_dir}: {e}") from e

            file_path = sum_dir / "project_summary.txt"
            content = "\n".join(f"{k}: {v}" for k, v in record.items()) + "\n"

            try:
                file_path.write_text(content, encoding="utf-8")
                print(f"Project summary saved to {file_path}")
            except PermissionError as e:
                raise PermissionError(f"Sin permisos para escribir {file_path}") from e
            except FileNotFoundError as e:
                raise FileNotFoundError(f"Directorio padre no existe para {file_path}") from e
            except IsADirectoryError as e:
                raise IsADirectoryError(f"{file_path} es un directorio, no un archivo") from e
            except UnicodeEncodeError as e:
                raise UnicodeEncodeError(e.encoding, e.object, e.start, e.end,
                                 f"Error de codificación al escribir {file_path}")
            except OSError as e:
                raise OSError(f"Error del sistema escribiendo {file_path}: {e}") from e
            
        return file_path


