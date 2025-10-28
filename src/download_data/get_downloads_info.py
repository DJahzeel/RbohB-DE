from Bio import Entrez
from pathlib import Path
from urllib.error import HTTPError, URLError
import pandas as pd
import time


def run_end(table, output):

    run_ends = {
        'SRR': table['Run'],
        'LibraryLayout': table['LibraryLayout']
    }

    filter = pd.DataFrame(run_ends)
    p = filter['LibraryLayout'] == 'PAIRED'
    s = filter['LibraryLayout'] == 'SINGLE'

    paired_end = pd.DataFrame(filter['SRR'][p])
    single_end = pd.DataFrame(filter['SRR'][s])

    print("Single-end reads:", single_end)
    print("Paired-end reads:", paired_end)
    print("Total Single-end:", len(single_end))
    print("Total Paired-end:", len(paired_end))
    
    outdir = Path(output)
    if outdir.exists() and not outdir.is_dir():
        raise NotADirectoryError(f"{outdir} existe pero no es un directorio")
    
    paired_file_path = outdir / "paired_end_runs.tsv"
    single_file_path = outdir / "single_end_runs.tsv"

    try:
        paired_end.to_csv(paired_file_path, index=False, sep='\t')
        single_end.to_csv(single_file_path, index=False, sep='\t')

    except PermissionError as e:
        raise PermissionError(f"Sin permisos para escribir ") from e
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Directorio padre no existe para ") from e
    except IsADirectoryError as e:
        raise IsADirectoryError(f" es un directorio, no un archivo") from e
    except UnicodeEncodeError as e:
        raise UnicodeEncodeError(e.encoding, e.object, e.start, e.end,
                                 f"Error de codificación al escribir ")
    except OSError as e:
        raise OSError(f"Error del sistema escribiendo: {e}") from e
    
    return paired_file_path, single_file_path

def get_assemblies_data(organism, email):
    try:
        Entrez.email = email
        handle = Entrez.esearch(db="assembly", term=organism, retmax=150)
        search_results = Entrez.read(handle)
        handle.close()
        assembly_uids = search_results["IdList"]

    except (HTTPError, URLError) as e:
        raise ConnectionError(f"An error occurred: {e}") from e
    except Exception as e:
        raise Exception(f"An unexpected error occurred: {e}") from e

    assembly_data = []
    for uid in assembly_uids:
        try:
            Entrez.email = email
            handle = Entrez.esummary(db="assembly", id=uid)
            summary = Entrez.read(handle)
            handle.close()
        except (HTTPError, URLError) as e:
            raise ConnectionError(f"An error occurred: {e}") from e
        except Exception as e:
            raise Exception(f"An unexpected error occurred: {e}") from e

        time.sleep(0.4)

        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        assembly_data.append(docsum)
    
    return assembly_data

def get_reference(assembly_data, output, organism):

    df_assembly = pd.DataFrame(assembly_data)
    df_selected = df_assembly[[
    "Organism",
    "AssemblyName",
    "AssemblyStatus",
    "AssemblyAccession",
    "RefSeq_category",
    "FtpPath_GenBank"
    ]]
    
    df_filtered = df_selected[df_selected["RefSeq_category"] != "na"].copy()
    filename_prefix = df_filtered['FtpPath_GenBank'].str.split('/').str[-1] 
    base_path = df_filtered['FtpPath_GenBank']
    fasta_url = base_path + '/' + filename_prefix + '_genomic.fna.gz'
    gff_url = base_path + '/' + filename_prefix + '_genomic.gff.gz'
    df_filtered.loc[:, 'FASTA FTP'] = fasta_url
    df_filtered.loc[:, 'GFF FTP'] = gff_url

    outdir = Path(output)
    if outdir.exists() and not outdir.is_dir():
        raise NotADirectoryError(f"{outdir} existe pero no es un directorio")
    
    file_path = outdir / f"{organism}_reference_genome.tsv"
    try:
        df_filtered.to_csv(file_path, sep='\t', index=False)
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