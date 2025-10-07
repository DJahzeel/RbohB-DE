from Bio import Entrez
import pandas as pd

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

def get_reference(assembly_data):
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

    print(" Bash commands to download the reference genome and annotation files:")
    print("wget -c", fasta_url)
    print("wget -c", gff_url)
    
    return df_filtered