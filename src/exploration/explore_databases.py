from Bio import Entrez
from urllib.error import HTTPError, URLError
import os
import time
import pandas as pd
from io import StringIO 

def fetch_bioproject_ids(bioprojects, databases, mail):
    try:
        found_ids = []
        for bioproject in bioprojects:
            for db in databases:
                Entrez.email = mail
                handle = Entrez.elink(db=db, dbfrom="bioproject", id=bioproject)
                search_results = Entrez.read(handle)
                handle.close()
                if 'LinkSetDb' in search_results[0] and len(search_results[0]['LinkSetDb']) > 0:
                    dict_ids = {"bioproject": bioproject, "database": db, "ids":search_results[0]['LinkSetDb'][0]['Link']}
                    found_ids.append(dict_ids)
                    print(f"Found results for Bioproject {bioproject} in database {db}: {dict_ids['ids']}")
                else:
                    print(f"No results found for Bioproject {bioproject} in database {db}.")
                    continue

                time.sleep(0.4)
                
        return found_ids
    
    except (HTTPError, URLError) as e:
        print(f"An error occurred: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None

def geo_links(record, mail):
    
    try:
        id_string = ','.join([id_dict['Id'] for id_dict in record['ids']])
        print(f"Fetching summaries for IDs: {id_string}")
        Entrez.email = mail
        handle = Entrez.efetch(db="gds", id=id_string, rettype="full", retmode="text")
        summaries = handle.read()
        handle.close()
        print(summaries)

        try:
            os.makedirs("data", exist_ok=True)
            os.makedirs(f"data/{record['bioproject']}/{record['database']}_files", exist_ok=True)
            file_path = f"data/{record['bioproject']}/{record['database']}_files/geo_links.txt"
            with open(file_path, "w", encoding="utf-8") as file:
                file.write(summaries)

        except FileNotFoundError as not_found_error:
            print(f"File not found error occurred: {not_found_error}")
            exit(1)
        except PermissionError as perm_error:
            print(f"Permission error occurred: {perm_error}")
            exit(1)
        except OSError as os_error:      
            print(f"OS error occurred: {os_error}")
            exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            exit(1)

    except (HTTPError, URLError) as e:
        print(f"An error occurred: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    
def gse_gsm(mail, record):
    try:
        Entrez.email = mail
        id_string = ','.join([id_dict['Id'] for id_dict in record['ids']])
        handle = Entrez.esummary(db="gds", id=id_string)
        summaries = Entrez.read(handle)
        handle.close()

        rows = []
        for summary in summaries:
            gse = summary.get('GSE') or summary.get('Accession')
            if gse and isinstance(gse, str) and not gse.startswith("GSE"):
                gse = f"GSE{gse}"
                        
                samples = summary.get('Samples') or []
                for sample in samples:
                    gsm  = sample.get('Accession')
                    if gsm:
                        rows.append({'GSE': gse, 'GSM': gsm})
        try:
            os.makedirs("data", exist_ok=True)
            os.makedirs(f"data/{record['bioproject']}/{record['database']}_files", exist_ok=True)
            file_path = f"data/{record['bioproject']}/{record['database']}_files/gse_gsm.tsv"
            with open(file_path, "w", encoding="utf-8") as file:
                df = pd.DataFrame(rows)
                df.to_csv(file, sep="\t", index=False)

          
        except FileNotFoundError as not_found_error:
            print(f"File not found error occurred: {not_found_error}")
            exit(1)
        except PermissionError as perm_error:
            print(f"Permission error occurred: {perm_error}")
            exit(1)
        except OSError as os_error:
            print(f"OS error occurred: {os_error}")
            exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            exit(1)

    except Exception as e:
        print(f"An error occurred while processing summaries: {e}")
    except (HTTPError, URLError) as e:
        print(f"An error occurred while fetching summaries: {e}")
    
    return rows

def sra_links(mail, record):
    try:
        id_string = ','.join([id_dict['Id'] for id_dict in record['ids']])
        Entrez.email = mail
        handle = Entrez.efetch(db="sra", id=id_string, rettype="runinfo", retmode="text")
        runinfo = handle.read()
        handle.close()
            
        if isinstance(runinfo, bytes):
            runinfo = runinfo.decode("utf-8")
        else:
            runinfo = str(runinfo)

        try:
            os.makedirs("data", exist_ok=True)
            os.makedirs(f"data/{record['bioproject']}/{record['database']}_files", exist_ok=True)
            file_path = f"data/{record['bioproject']}/{record['database']}_files/sra_runinfo.tsv"

            table = pd.read_csv(StringIO(runinfo))
            columns = {
                'Run': table.get('Run'),
                'Experiment': table.get('Experiment'),
                'Platform': table.get('Platform'),
                'LibraryName': table.get('LibraryName'),
                'LibraryLayout': table.get('LibraryLayout'),
                'Sample': table.get('Sample'),
                'ScientificName': table.get('ScientificName'),
                'SampleName': table.get('SampleName'),
            }

            df = pd.DataFrame(columns)
            df.to_csv(file_path, index=False, sep='\t')

        except FileNotFoundError as not_found_error:
            print(f"File not found error occurred: {not_found_error}")
            exit(1)
        except PermissionError as perm_error:
            print(f"Permission error occurred: {perm_error}")
            exit(1)
        except OSError as os_error:
            print(f"OS error occurred: {os_error}")
            exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            exit(1)
        return df

    except (HTTPError, URLError) as e:
        print(f"An error occurred: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")