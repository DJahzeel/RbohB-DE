from Bio import Entrez
from urllib.error import HTTPError, URLError
import os

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
        print(f"Connection error: {connection_error}")
        return None
    except Exception as error:
        print(f"An error occurred: {error}")
        return None

def save_project_summaries(records):
    try:
        if not records:
            print("No records to save.")
            return
    
        os.makedirs("data", exist_ok=True)
    
        for record in records:
            project_uid = record.get('Project_Id')

            if not project_uid:
                print("No Project ID found in the record, skipping.")
                continue

            os.makedirs(f"data/{project_uid}", exist_ok=True)
            file_path = f"data/{project_uid}/project_summary.txt"

            try:
                with open(file_path, 'w', encoding='utf-8') as file:
                    for key, value in record.items():
                        file.write(f"{key}: {value}\n")
                        print(f"Project summary saved to {file_path}")
                        
            except IOError as io_error:
                print(f"File error: {io_error}")
                return None

    except Exception as error:
        print(f"An error occurred: {error}")


