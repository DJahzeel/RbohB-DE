## Resumen del artículo / Proyecto basado en el estudio

### Organismo de referencia
El trabajo se centra en *Phaseolus vulgaris* (frijol común), variedad Mesoamericana cv. Negro Jamapa, un cultivo de interés agrícola y modelo para estudiar interacciones simbióticas. Se analizan raíces inoculadas con el hongo micorrízico arbuscular (*Rhizophagus irregularis*) y con la bacteria fijadora de nitrógeno (*Rhizobium tropici*).

### Objetivos  
Con base en el artículo citado, pero definiendo objetivos propios de replicación computacional, planteamos lo siguiente:
### Objetivo general
Establecer y ejecutar un flujo de trabajo bioinformático robusto y reproducible para RNA-seq en Phaseolus vulgaris, a partir de los datos públicos del estudio de referencia.

### Objetivos especificos
\-Descarga de datos: Obtener los archivos FASTQ de las muestras y el genoma de referencia de Phaseolus vulgaris.

\-Alineamiento: Mapear las lecturas contra el genoma de referencia y generar archivos BAM (ordenados e indexados).

\-Conteo de lecturas: Construir la matriz de conteos por gen a partir de los alineamientos.

\-Expresión diferencial: Identificar genes diferencialmente expresados entre condiciones usando DESeq2.

\-Visualización: Representar los resultados mediante PCA, volcano plots, heatmaps y gráficas de expresión de RbohB.

\-Anotación funcional: Obtener información adicional de los genes mediante búsquedas en bases de datos como BLAST.

\-Interpretación: Analizar e integrar los resultados en el contexto biológico del estudio.

### Resultados esperados
\-Datos & QC: FASTQ y referencia descargados/validados; reporte de integridad y legibilidad.

\-Alineamiento: BAM ordenados e indexados; tasas de mapeo en rangos esperados (flagstat/qualimap).

\-Conteos: Matriz genes×muestras con alta correlación entre réplicas.

\-DESeq2: Tablas de DEGs (FDR < 0.05); respuesta transcriptómica mayor en Rhizobium que en AM.

\-Figuras: PCA con separación clara por condición; volcano/heatmaps consistentes; perfil de PvRbohB acorde a los contrastes.

\-Funcional: Enriquecimientos GO/KEGG en ROS, pared celular y vías hormonales (auxina, citoquinina, etileno); BLAST para genes clave.

\-Reproducibilidad: Repo ejecutable end-to-end (README, environment, scripts/logs).

En conjunto, se anticipa demostrar que RbohB promueve la infección y nodulación por R. tropici pero limita la colonización micorrízica de R. irregularis.

### BioProject ID  
Los datos generados por el proyecto (lecturas de RNA-Seq, muestras biológicas y metadatos) fueron depositados en la base de datos NCBI bajo el BioProject ID **PRJNA482464**
