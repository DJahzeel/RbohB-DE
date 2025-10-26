from Bio import Entrez
import os
import pandas as pd

def run_end(table, record):

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
    
    
    os.makedirs("data", exist_ok=True)
    os.makedirs(f"data/{record['bioproject']}/{record['database']}_files", exist_ok=True)

    paired_file_path = f"data/{record['bioproject']}/{record['database']}_files/paired_end_runs.tsv"
    single_file_path = f"data/{record['bioproject']}/{record['database']}_files/single_end_runs.tsv"

    paired_end.to_csv(paired_file_path, index=False, sep='\t')
    single_end.to_csv(single_file_path, index=False, sep='\t')
    
    return paired_end, single_end

def get_assemblies_data(organism):
    assembly_data = []
    handle = Entrez.esearch(db="assembly", term=organism, retmax=150)
    search_results = Entrez.read(handle)
    handle.close()
    assembly_uids = search_results["IdList"]

    for uid in assembly_uids:
        handle = Entrez.esummary(db="assembly", id=uid)
        summary = Entrez.read(handle)
        handle.close()

        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        assembly_data.append(docsum)
    
    return assembly_data

def get_reference(assembly_data, organism):

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

    os.makedirs("data/references", exist_ok=True)
    df_filtered.to_csv(f'data/{organism}_reference_genome.tsv', sep='\t', index=False)

        
    return df_filtered