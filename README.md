# Análisis transcriptómico de raíces de Phaseolus vulgaris (PvRbohB) bajo simbiosis rizobiana y micorrízica arbuscular
**Fecha:** 6/10/2025

**Autores:** 
- Pedro Daniel Pineda Martinez
email: pedropm@lcg.unam.mx
- Dara Jahzeel Palafox Collado
email: darapc@lcg.unam.mx

## Descripción
Este repositorio contiene un pipeline reproducible de análisis transcriptómico (RNA-seq) aplicado a raíces de Phaseolus vulgaris a 7 días post-inoculación (dpi), bajo colonización por Rhizobium tropici (nodulación) y Rhizophagus irregularis (micorrización arbuscular, AM), en líneas WT y PvRbohB-RNAi.

El proyecto se basa en los datos públicos del BioProject PRJNA482464 y en la información biológica reportada por Fonseca-García et al. (2019), pero no busca replicar de manera exacta todos los análisis del estudio original.
En cambio, el objetivo es:

- implementar un pipeline de RNA-seq con las herramientas revisadas en clase, y
- explorar si los resultados generales son coherentes con lo descrito para PvRbohB, ROS, pared celular y fitohormonas.

## Objetivos Clave
Desarrollar e implementar un pipeline modular y reproducible para el análisis transcriptómico de raíces de P. vulgaris (7 dpi), utilizando los datos de PRJNA482464 como caso de estudio, para obtener una matriz de conteos y contrastes de expresión diferencial que permitan explorar el papel de PvRbohB en contextos simbióticos.

**Objetivos específicos**
- Generar alineamientos de calidad y una matriz de conteos integrada de raíces de P. vulgaris (7 dpi), documentando métricas de QC en cada etapa del pipeline.
- Comparar perfiles de expresión entre nodulación (R. tropici) y micorrización (R. irregularis) en WT y PvRbohB-RNAi, usando DESeq2 (u otra herramienta de DE) a partir de la matriz de conteos generada.
- Identificar genes diferencialmente expresados por contraste (p.ej. padj ≤ 0.05) y producir tablas y figuras canónicas (MA, volcán, PCA) con una anotación mínima.
- Documentar la estructura del repositorio, los comandos y versiones de software, así como una comparación cualitativa entre las tendencias observadas en este proyecto y las reportadas en el estudio de referencia.

## Preguntas de Investigación
- ¿Qué diferencias transcripcionales existen entre las condiciones (Control, Nodulación, AM) a 7 dpi en raíces de P. vulgaris?
- ¿Cómo modula la deleción de PvRbohB-RNAi estas respuestas transcripcionales respecto a WT, particularmente en los ejes de ROS, pared celular y fitohormonas?
- ¿En qué medida los resultados replicados concuerdan con el estudio base y qué variaciones podrían explicar las posibles desviaciones?

## Datos 

| Característica | Detalle |
| :--- | :--- |
| **BioProject** | PRJNA482464 (Datos públicos del estudio base) |
| **Organismo** | *Phaseolus vulgaris* (cv. Negro Jamapa) |
| **Referencia** | *P. vulgaris* (G19833), **NCBI GCA\_000499845.2** |
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
FASTQ ──► FastQC / MultiQC ──► (opcional) Cutadapt ──► BWA + samtools ──► coverageBed + matriz de conteos ──► DESeq2 ──► anotación funcional
 lecturas crudas       control de calidad        trimming         alineamiento BAM/BAI       conteos por gen               DEGs + figuras          interpretación


```

### Pasos de alto nivel
1. **Identificación y obtención de datos**: RNA-seq (Illumina) del BioProject PRJNA482464 (SRR por condición/genotipo).
2. **Preparación**: descarga (SRA Toolkit) y genoma de referneica *P. vulgaris* (G19833; GCA_000499845.1) , organización de carpetas y metadatos (`samples.tsv`).
4. **Control de calidad**:   Evaluación de la calidad de los FASTQ crudos con **FastQC** y resumen global con **MultiQC**.(Opcional) Aplicar **trimming**  y repetir QC sobre las lecturas recortadas.
5. . **Alineamiento**: `bwa mem` + `samtools sort/index` (respetando **PAIRED/SINGLE** según grupo).
6. . **Cuantificación**:  Uso de **coverageBed** (bedtools) con la anotación GFF para obtener conteos por gen y muestra (`<sample>.counts.txt`), seguido de la integración de todos los archivos en una **matriz de conteos** (genes × muestras) mediante un script en Python..
7.  **DE**: DESeq2 por contrastes (WT y PvRbohB-RNAi: **Rhiz vs CTRL**, **AM vs CTRL**).
8.   **Anotación mínima**: mapeo `gene_id → descripción` a partir del GFF (GO/KEGG opcional).
9.    **Síntesis**: figuras (MA/volcán/PCA), tablas de DEGs y nota de concordancia con el estudio base.

## Requisitos funcionales
- Procesar lecturas RNA-seq y **generar alineamientos** reproducibles.
- Cuantificar expresión por gen a partir de BAM y anotación GFF.
- Generar una matriz de conteos lista para análisis de expresión diferencial.
- Producir tablas y figuras básicas (DEGs, MA, volcán, PCA).
- Mantener trazabilidad de parámetros y versiones (README).


##  Resultados y entregables esperados
| Resultados Esperados | Entregables (Archivos) |
| :--- | :--- |
| **Matriz de expresión** de alta calidad y métricas de mapeo documentadas. | `results/counts/all_counts.tsv`, `results/qc/` (métricas). |
| Listas de **Genes Diferencialmente Expresados (DEGs)** ($\text{padj} \le 0.05$). | `results/DE/DEGs.tsv`. |
| Tendencias biológicas coherentes en ejes de ROS, pared celular y fitohormonas. | Figuras: `results/DE/MA.png`, `volcano.png`, `PCA.png`. |
| Síntesis comparativa WT vs *PvRbohB-RNAi*. | Notas en `docs/`, `README.md` final,Archivos de anotación. |



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

                
