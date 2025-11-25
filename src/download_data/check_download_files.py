import pandas as pd

'''
Este módulo contiene funciones para verificar la existencia de archivos FASTQ descargados
y generar una hoja de muestras basada en los archivos presentes.
'''
def check_paired_files(paired_end, fastqdir):
    verified_rows = []
    for srr in (list(paired_end)):
        file1 = fastqdir / f"{srr}_1.fastq.gz"
        file2 = fastqdir / f"{srr}_2.fastq.gz"
        if file1.exists() and file2.exists():
            verified_rows.append({
                'sample_id': srr,
                'file_1': file1.name,
                'file_2': file2.name
                })
        elif file1.exists() and not file2.exists():
            verified_rows.append({
                'sample_id': srr,
                'file_1': file1.name,
                'file_2': "MISSING"
                })
        elif not file1.exists() and file2.exists():
            verified_rows.append({
                'sample_id': srr,
                'file_1': "MISSING",
                'file_2': str(file2.stem)
            })
        else:
            print(f"Advertencia: Ambos archivos faltan para el SRR(paired end) {srr}. No se incluirá en el archivo de salida.")
    return verified_rows

def check_single_files(single_end, fastqdir):
    verified_rows = []
    for srr in (list(single_end)):
        file1 = fastqdir / f"{srr}.fastq.gz"
        if file1.exists():
            verified_rows.append({
                'sample_id': srr,
                'file_1': file1.name,
                'file_2': ""
                })
        else:
            print(f"Advertencia: El archivo falta para el SRR(single end) {srr}. No se incluirá en el archivo de salida.")
    return verified_rows

def save_sample_sheet(verified_rows, output, processing):

    file_path = output / f"{processing}_samples.tsv"
    try:
        df_verified = pd.DataFrame(verified_rows, columns=['sample_id', 'file_1', 'file_2'])
        df_verified.to_csv(file_path, sep='\t', index=False)
        return file_path
    except PermissionError as e:
        raise PermissionError(f"Sin permisos para escribir en {file_path}") from e 
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Directorio padre no existe para {file_path}") from e
    except OSError as e:
        raise OSError(f"Error al escribir el archivo {file_path}: {e}") from e
