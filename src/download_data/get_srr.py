from Bio import Entrez
from io import StringIO
import pandas as pd

def srr_ids(sra_ids):

    sra_id_list = ",".join(sra_ids)
    handle = Entrez.efetch(db="sra", id=sra_id_list, rettype="runinfo", retmode="text")
    run_info= handle.read().decode('utf-8')
    handle.close()
    
    return run_info

def run_end(run_info):
    reads_table = pd.read_csv(StringIO(run_info))
    run_ends = {
        'SRR': reads_table['Run'],
        'LibraryLayout': reads_table['LibraryLayout']
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

    return paired_end, single_end

def save_srr(paired_end, single_end, file_name):
    paired_end.to_csv(f'{file_name}_paired.tsv', sep='\t', index=False, header=False)
    single_end.to_csv(f'{file_name}_single.tsv', sep='\t', index=False, header=False)
