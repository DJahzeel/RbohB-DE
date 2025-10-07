
from Bio import Entrez
import argparse as ap

def fetch_project_id(accesion,mail):
    Entrez.email = mail
    handle = Entrez.esearch(db="bioproject", term=accesion)
    search_results = Entrez.read(handle)
    print(search_results)
    handle.close()
    return search_results

def recover_uid(search_results):
    project_uid = search_results['IdList'][0]
    print(f"Project UID: {project_uid}")
    return project_uid

def project_summary(project_uid):
    handle = Entrez.esummary(db="bioproject", id=project_uid)
    proj_record = Entrez.read(handle)
    handle.close()
    proj = proj_record["DocumentSummarySet"]["DocumentSummary"][0]
    for key, value in proj.items():
        print(f"{key}: {value}")


