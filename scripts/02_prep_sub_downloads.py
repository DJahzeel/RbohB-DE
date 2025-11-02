import argparse as ap
from download_data.get_downloads_info import save_run_ends
from paths.pathsval import project_root, runends_directory
from utils.job_submission import submit_fasterqdump_job
from download_data.get_downloads_info import run_end
import pandas as pd

parser = ap.ArgumentParser(
    description=("Prepara y envía trabajos de descarga de datos de secuenciación "
                 "utilizando fasterq-dump para un BioProject específico.")
)
parser.add_argument("-b", "--bioproject",
                    required=True,
                    help="Accesión del BioProject (UID Entrez del Bioproject e.g., 168994)")
parser.add_argument("-t", "--type",
                    required=True,
                    choices=["GEO", "SRA"],
                    help="Descargar datos asociados a GSEs específicos(muestras)(si existe relacion) " \
                    "o a todo el SRA del BioProject.")
parser.add_argument("-g", "--gse",
                    required=False,
                    nargs="+",
                    help="Lista de GSEs separados por comas (e.g., GSE12345,GSE67890). " \
                    "Solo requerido si el tipo es 'GEO'.")
args = parser.parse_args()

def main():
    print("BioProject:", args.bioproject)
    print("Tipo de descarga:", args.type)
    if args.gse:
        print("GSEs:", args.gse)
        gse_list = ",".join(args.gse)

    base_project_path = project_root()
    bpdir = base_project_path / "data"  / args.bioproject
    sra_files_dir = bpdir / "sra_files"
    gds_files_dir = bpdir / "gds_files"

    sra_info_path = sra_files_dir / "sra_runinfo.tsv"
    gds_info_path = gds_files_dir / "gse_gsm.tsv"

    if not sra_info_path.is_file():
        raise FileNotFoundError(f"SRA runinfo file not found at {sra_info_path}")

    print (f"Reading SRA runinfo from {sra_info_path}")
    

    if args.type == "GSE":
        try:
            sra_df = pd.read_csv(sra_info_path, sep="\t")
            if sra_df.empty:
                raise ValueError("SRA runinfo file is empty.")
            if 'Run' not in sra_df.columns or 'LibraryLayout' not in sra_df.columns:
                raise ValueError("SRA runinfo file does not contain 'Run' or 'LibraryLayout' column.")
        except Exception as e:
            raise ValueError(f"Error reading SRA runinfo file: {e}")

        try:
            df_gse_gsm = pd.read_csv(gds_info_path, sep=",")
            if df_gse_gsm.empty:
                raise ValueError("GDS info file is empty.")
            if 'GSE' not in df_gse_gsm.columns or 'GSM' not in df_gse_gsm.columns:
                raise ValueError("GDS info file does not contain 'GSE' or 'GSM' column.")
        except Exception as e:
            raise ValueError(f"Error reading GDS info file: {e}")

   
        for gse in gse_list.split(","):
            gse = gse.strip()
            gsm_list = df_gse_gsm[df_gse_gsm['GSE'] == gse]['GSM'].tolist()
            sra_subset = sra_df[sra_df['SampleName'].isin(gsm_list)]
            if sra_subset.empty:
                print(f"No matching GSM entries found in SRA runinfo for {gse}. Skipping.")
                continue

            single, paired = run_end(sra_subset)
            output_path = runends_directory(base_project_path, args.bioproject)
            paired_path, single_path = save_run_ends(paired, single, output_path, gse)

            print(f"Processed {gse}:")
            print(f"  Paired-end runs saved to: {paired_path}")
            print(f"  Single-end runs saved to: {single_path}")

    elif args.type == "SRA":
        try:
            sra_df = pd.read_csv(sra_info_path, sep="\t")
            if sra_df.empty:
                raise ValueError("SRA runinfo file is empty.")
            if 'Run' not in sra_df.columns or 'LibraryLayout' not in sra_df.columns:
                raise ValueError("SRA runinfo file does not contain 'Run' or 'LibraryLayout' column.")
        except Exception as e:
            raise ValueError(f"Error reading SRA runinfo file: {e}")
         
        single, paired = run_end(sra_df)
        output_path = runends_directory(base_project_path, args.bioproject)
        paired_path, single_path = save_run_ends(paired, single, output_path, args.bioproject)

        print(f"Processed SRA project {args.bioproject}:")
        print(f" Paired-end runs saved to: {paired_path}")
        print(f" Single-end runs saved to: {single_path}")

    paired_script = base_project_path / "src" / "fasterqd_bash" /"sra_paired.jdl"
    single_script = base_project_path / "src" / "fasterqd_bash" /"sra_single.jdl"

    if not paired_path.exists() and not single_path.exists():
        raise RuntimeError("No se crearon archivos de salida. No se mandaron trabajos de descarga.")
    
    if paired_path.exists() and not single_path.exists():
        print("Solo se crearon archivos de pares. Se mandarán trabajos de descarga para pares.")
        submit_fasterqdump_job(f"{args.bioproject}_paired_download", paired_script, paired_path, wait_for=False)

    if single_path.exists() and not paired_path.exists():
        print("Solo se crearon archivos de simples. Se mandarán trabajos de descarga para simples.")
        submit_fasterqdump_job(f"{args.bioproject}_single_download", single_script, single_path, wait_for=False)

    if paired_path.exists() and single_path.exists():
        print("Se crearon archivos de pares y simples. Se mandarán trabajos de descarga para ambos secuencialmente.")
        submit_fasterqdump_job(f"{args.bioproject}_single_download", single_script, single_path, wait_for=True)
        submit_fasterqdump_job(f"{args.bioproject}_paired_download", paired_script, paired_path, wait_for=False)

