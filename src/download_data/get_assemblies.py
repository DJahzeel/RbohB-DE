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

    df_filtered = df_selected[df_selected["RefSeq_category"] != "na"]
    df_filtered['FASTA FTP'] = df_filtered['FtpPath_GenBank'] + '/' + '_genomic.fna.gz'
    df_filtered['GFF FTP'] = df_filtered['FtpPath_GenBank'] + '/' + '_genomic.gff.gz'
    return df_filtered