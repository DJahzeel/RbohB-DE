# Análisis transcriptómico de raíces de Phaseolus vulgaris (PvRbohB) bajo simbiosis rizobiana y micorrízica arbuscular
**Fecha:** 6/10/2025

**Autores:** 
- Pedro Daniel Pineda Martinez
email: pedropm@lcg.unam.mx
- Dara Jahzeel Palafox Collado
email: darapc@lcg.unam.mx

## Descripción
Replicación reproducible del análisis transcriptómico de raíces de *Phaseolus vulgaris* a **7 días post-inoculación** (dpi), comparando colonización por **Rhizobium tropici** (nodulación) y **Rhizophagus irregularis** (micorrización arbuscular, AM) en líneas **WT** y **PvRbohB-RNAi**, con énfasis en redes de **ROS**, **pared celular** y **fitohormonas**.

## Objetivos Clave
El objetivo principal es **replicar** de forma reproducible el análisis transcriptómico de raíces de *Phaseolus vulgaris*  del estudio base (Fonseca-García et al., 2019), recuperando los contrastes y las salidas clave reportadas, y documentando cualquier desviación técnica o biológica.
### Especificos :
- Generar alineamientos de calidad y una matriz de conteos integrada de raíces de *P. vulgaris* (7 dpi), con métricas de QC documentadas.
- Comparar perfiles de expresión entre nodulación (R. tropici) y micorrización (R. irregularis) en WT y PvRbohB-RNAi mediante DESeq2.
- Identificar genes diferencialmente expresados por contraste (padj ≤ 0.05) y entregar tablas y figuras canónicas (MA, volcán, PCA) con anotación mínima.
- Producir recursos reproducibles (estructura del repo, comandos y versiones) y un breve resumen de concordancia con el estudio base

## Preguntas de Investigación
- ¿Qué diferencias transcripcionales existen entre las condiciones (Control, Nodulación, AM) a 7 dpi en raíces de P. vulgaris?

- ¿Cómo modula la deleción de PvRbohB-RNAi estas respuestas transcripcionales respecto a WT, particularmente en los ejes de ROS, pared celular y fitohormonas?

- ¿En qué medida los resultados replicados concuerdan con el estudio base y qué variaciones podrían explicar las posibles desviaciones?

## Datos 

| Característica | Detalle |
| :--- | :--- |
| **BioProject** | PRJNA482464 (Datos públicos del estudio base) |
| **Organismo** | *Phaseolus vulgaris* (cv. Negro Jamapa) |
| **Referencia** | *P. vulgaris* (G19833), **NCBI GCA\_000499845.1** |
| **Diseño** | **2 genotipos** ($\text{WT}, \text{PvRbohB-RNAi}$) $\times$ **3 condiciones** ($\text{Control}, \text{Rhizobium tropici}, \text{Rhizophagus irregularis} \text{ [AM]}$) $\times$ **3 réplicas** $= \mathbf{18}$ corridas. |
| **Tipo de datos** | RNA-seq (Illumina); *Layout* **mixto** (*PAIRED* en la mayoría, *SINGLE* solo en PvRbohB-RNAi + Rhizobium). |

### Accesiones por Condición (SRR)

| Condición | Genotipo | *Layout* | Accesiones SRR |
| :--- | :--- | :--- | :--- |
| **Control** | WT | PAIRED | 7693915, 7693916, 7693917 |
| **$R. tropici$** | WT | PAIRED | 7696192, 7696193, 7696194 |
| **$R. irregularis$** | WT | PAIRED | 7696200, 7696201, 7696202 |
| **Control** | PvRbohB-RNAi | PAIRED | 7696204, 7696205, 7696206 |
| **$R. tropici$** | PvRbohB-RNAi | SINGLE | 7696589, 7696590, 7696591 |
| **$R. irregularis$** | PvRbohB-RNAi | PAIRED | 7696208, 7696209, 7696210 |


## Flujo de trabajo

```
FASTQ ──► FastQC ──► BWA+samtools ──► htseq-count ──► DESeq2 ──► (anotación opcional)
         QC           BAM/BAl          matriz de conteos        DEGs + figuras
```

### Pasos de alto nivel
1. **Identificación y obtención de datos**: RNA-seq (Illumina) del BioProject PRJNA482464 (SRR por condición/genotipo).
2. **Preparación**: descarga (SRA Toolkit), organización de carpetas y metadatos (`samples.tsv`).
3. **Control de calidad**: FastQC (opcional: trimming si aplica).
4. **Referencia**: genoma/anotación de *P. vulgaris* (G19833; GCA_000499845.1) e índices BWA.
5. **Alineamiento**: `bwa mem` + `samtools sort/index` (respetando **PAIRED/SINGLE** según grupo).
6. **Cuantificación**: `htseq-count` por gen (definiendo `-s`, `-t gene`, `-i ID`).
7. **DE**: DESeq2 por contrastes (WT y PvRbohB-RNAi: **Rhiz vs CTRL**, **AM vs CTRL**).
8. **Anotación mínima**: mapeo `gene_id → descripción` a partir del GFF (GO/KEGG opcional).
9. **Síntesis**: figuras (MA/volcán/PCA), tablas de DEGs y nota de concordancia con el estudio base.

## Requisitos funcionales
- Procesar lecturas RNA-seq y **generar alineamientos** reproducibles.
- **Cuantificar** expresión por gen y realizar **DE** por contraste a 7 dpi.
- Producir **tablas/figuras** (DEGs, MA, volcán, PCA) y una **anotación mínima** de genes.
- Mantener **trazabilidad** de parámetros y versiones (README + `VERSIONS.txt`).


##  Resultados y entregables esperados
| Resultados Esperados | Entregables (Archivos) |
| :--- | :--- |
| **Matriz de expresión** de alta calidad y métricas de mapeo documentadas. | `results/counts/all_counts.tsv`, `results/qc/` (métricas). |
| Listas de **Genes Diferencialmente Expresados (DEGs)** ($\text{padj} \le 0.05$). | `results/DE/DEGs.tsv`. |
| Tendencias biológicas coherentes en ejes de ROS, pared celular y fitohormonas. | Figuras: `results/DE/MA.png`, `volcano.png`, `PCA.png`. |
| Síntesis comparativa WT vs *PvRbohB-RNAi* y **documentación de desvíos**. | `results/annotation/annotation.tsv`, `README.md` final, `VERSIONS.txt`. |



## Estructura del repositorio 
````
.
├── README.md 
├── .gitignore  # Archivos que no deben de ser versionados
├── data/       # FASTQs, metadatos, datos de entrada
├── results/    # QC, alineamientos, conteos, DE, figuras
├── scripts/    # Scripts de Bash/R/Python que controlan el flujo del pipeline
├── src/        # Paquetes y modulos
├── docs/       # Documentación adicional o notas de desvío
├── pyproject.toml # Archivo de configuracion para el proyecto     
└── LICENSE                         
````

## Referencias
Fonseca-García, C., Zayas, A. E., Montiel, J., Nava, N., Sánchez, F., & Quinto, C. (2019). Transcriptome analysis of the differential effect of the NADPH oxidase gene RbohB in Phaseolus vulgaris roots following Rhizobium tropici and Rhizophagus irregularis inoculation. BMC Genomics, 20, 800. https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6162-7

                
