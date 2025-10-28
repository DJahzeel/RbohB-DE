from Bio import Entrez
from urllib.error import HTTPError, URLError
import time
import pandas as pd
from io import StringIO
from pathlib import Path

def fetch_bioproject_links(bioprojects, databases, mail):
    found_ids = []
    try:
        for bioproject in bioprojects:
            for db in databases:
                Entrez.email = mail
                handle = Entrez.elink(db=db, dbfrom="bioproject", id=bioproject)
                search_results = Entrez.read(handle)
                handle.close()
                if 'LinkSetDb' in search_results[0] and len(search_results[0]['LinkSetDb']) > 0:
                    dict_ids = {"bioproject": bioproject, "database": db, "ids":search_results[0]['LinkSetDb'][0]['Link']}
                    found_ids.append(dict_ids)
                else:
                    print(f"No results found for Bioproject {bioproject} in database {db}.")
                    continue

                time.sleep(0.4)
                
    except (HTTPError, URLError) as e:
       raise ConnectionError(f"An error occurred: {e}") from e
    except Exception as e:
       raise Exception(f"An unexpected error occurred: {e}") from e

    return found_ids

def geo_links(record, mail, output):

    try:

        id_string = ','.join([id_dict['Id'] for id_dict in record['ids']])
        Entrez.email = mail
        handle = Entrez.efetch(db="gds", id=id_string, rettype="full", retmode="text")
        summaries = handle.read()
        handle.close()

    except (HTTPError, URLError) as e:
        raise ConnectionError(f"An error occurred: {e}")
    except Exception as e:
        raise Exception(f"An unexpected error occurred: {e}")

    outdir = Path(output)
    if outdir.exists() and not outdir.is_dir():
        raise NotADirectoryError(f"{outdir} existe pero no es un directorio")
    
    file_path = outdir / "geo_links.txt"

    try:
        file_path.write_text(summaries, encoding="utf-8")
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

def gse_gsm(record, mail, output):

    try:
        Entrez.email = mail
        id_string = ','.join([id_dict['Id'] for id_dict in record['ids']])
        handle = Entrez.esummary(db="gds", id=id_string)
        summaries = Entrez.read(handle)
        handle.close()

    except Exception as e:
        raise Exception(f"An error occurred while processing summaries: {e}") from e
    except (HTTPError, URLError) as connection_error:
        raise ConnectionError(f"An error occurred while fetching summaries: {connection_error}") from connection_error

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

    outdir = Path(output)
    if outdir.exists() and not outdir.is_dir():
        raise NotADirectoryError(f"{outdir} existe pero no es un directorio")

    file_path = outdir / "gse_gsm.tsv"
    try:
        df = pd.DataFrame(rows)
        df.to_csv(file_path, index=False, encoding="utf-8")

    except PermissionError as perm_error:
        raise PermissionError(f"Permission error occurred: {perm_error}") from perm_error
    except FileNotFoundError as not_found_error:
        raise FileNotFoundError(f"File not found error occurred: {not_found_error}") from not_found_error
    except IsADirectoryError as dir_error:
        raise IsADirectoryError(f"Is a directory error occurred: {dir_error}") from dir_error
    except UnicodeEncodeError as encode_error:  
        raise UnicodeEncodeError(encode_error.encoding, encode_error.object, encode_error.start, encode_error.end,
                                 f"Encoding error occurred: {encode_error}") from encode_error
    except OSError as os_error:
        raise OSError(f"OS error occurred: {os_error}") from os_error
    
    return file_path

def sra_links(mail, record, output):

    try:
        id_string = ','.join([id_dict['Id'] for id_dict in record['ids']])
        Entrez.email = mail
        handle = Entrez.efetch(db="sra", id=id_string, rettype="runinfo", retmode="text")
        runinfo = handle.read()
        handle.close()

    except (HTTPError, URLError) as e:
        raise ConnectionError(f"An error occurred: {e}") from e
    except Exception as e:
        raise Exception(f"An unexpected error occurred: {e}") from e
            
    if isinstance(runinfo, bytes):
        runinfo = runinfo.decode("utf-8")
    else:
        runinfo = str(runinfo)

    outdir = Path(output)
    if outdir.exists() and not outdir.is_dir():
        raise NotADirectoryError(f"{outdir} existe pero no es un directorio")

    file_path = outdir / "sra_runinfo.tsv"
    try:
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

    except PermissionError as perm_error:
        raise PermissionError(f"Permission error occurred: {perm_error}") from perm_error
    except FileNotFoundError as not_found_error:
        raise FileNotFoundError(f"File not found error occurred: {not_found_error}") from not_found_error
    except IsADirectoryError as dir_error:
        raise IsADirectoryError(f"Is a directory error occurred: {dir_error}") from dir_error
    except UnicodeEncodeError as encode_error:  
        raise UnicodeEncodeError(encode_error.encoding, encode_error.object, encode_error.start, encode_error.end,
                                 f"Encoding error occurred: {encode_error}") from encode_error
    except OSError as os_error:
        raise OSError(f"OS error occurred: {os_error}") from os_error
    
    return df, file_path
