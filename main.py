import argparse as ap
import os

from entrez_queries import explore_project as ep 
from entrez_queries import explore_databases as ed
from download_data import get_srr as srr
from download_data import get_assemblies as ga

    
parser = ap.ArgumentParser(description="Fetch project data from NCBI Bioproject database")
parser.add_argument('-a','--accession',
                    type=str,
                    required=True,
                    help="Project accession number (e.g., PRJNA123456)")
parser.add_argument('-m','--mail',
                    type=str,
                    required=True,
                    help="Email address for NCBI Entrez")
parser.add_argument('-e','--explore',
                    type=str,
                    required=False,
                    help="Explore type (database or only project)")
parser.add_argument('-d','--database',
                    type=str,
                    required=False,
                    default='GEO',
                    help="Database to explore (GEO or SRA)")
parser.add_argument('--srr',
                    type=bool,
                    required=False,
                    default=False,
                    help="Get SRR ids from SRA database")
parser.add_argument('--file',
                    type=bool,
                    required=False,
                    default=False,
                    help="File to save results")
parser.add_argument('--reference',
                    type=str,
                    required=False,
                    help="Reference genome")

args = parser.parse_args()

def project():
    print("============================== PROJECT FETCH RESULTS ==============================")
    fetch_results = ep.fetch_project_id(args.accession, args.mail)
    print("============================== PROJECT UID =========================================")
    project_uid = ep.recover_uid(fetch_results)
    return project_uid

def databases(project_uid):
    geo_ids = None
    sra_ids = None
    print("============================== EXPLORE DATABASES ===================================")
    if args.database == 'GEO' or args.database == 'geo':
        linksBioProj = ed.explore_databases(project_uid, "gds")
        geo_ids = ed.geo(linksBioProj, args.accession)
        if geo_ids is None:
            print("\nEl proyecto", args.accession, " NO tiene enlaces directos a GEO")
    print("============================== SRA DATABASE =========================================")
    if geo_ids is None or args.database == 'SRA' or args.database == 'sra':
        linksBioProj = ed.explore_databases(project_uid, "sra")
        sra_ids = ed.sra(linksBioProj, args.accession)
        if sra_ids is None:
            print("\nEl proyecto", args.accession, " NO tiene enlaces directos a SRA")
    if geo_ids is None and sra_ids is None:
        print("\nEl proyecto", args.accession, " NO tiene enlaces directos a GEO ni a SRA")
    return geo_ids, sra_ids

def main():
    project_uid = project()
    if not args.explore or args.explore.lower() == 'project':
        print("============================== PROJECT SUMMARY =====================================")
        ep.project_summary(project_uid)
    if not args.explore or args.explore.lower() == 'database':
        geo_ids, sra_ids = databases(project_uid)
    
    if args.srr and (sra_ids or geo_ids):
        ids = sra_ids if sra_ids else geo_ids
        srr_ids = srr.srr_ids(ids)
        paired, single = srr.run_end(srr_ids)
        if args.file:
            srr.save_srr(paired, single, args.accession)
    if args.reference:
        reference_data = ga.get_assemblies_data(args.reference)
        if reference_data is not None:
            reference_genome = ga.get_reference(reference_data)
            if args.file:
                os.makedirs('data', exist_ok=True)
                reference_genome.to_csv(f'data/{args.reference}_reference_genome.tsv', sep='\t', index=False)
        else :
            print(f"No assemblies/reference genome found for organism: {args.reference}")

if __name__ == "__main__":
    main()