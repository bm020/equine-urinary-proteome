# Defining the Equine Urinary Proteome: A Reference Baseline for Biomarker Discovery
Overview

This repository contains the R scripts used for the bioinformatics, statistical analyses and figure generation accompanying the manuscript:

Defining the Equine Urinary Proteome: A Reference Baseline for Biomarker Discovery

The study established a reference equine urinary proteome using high-resolution LC–MS/MS and demonstrates a workflow for urinary proteomics and biomarker discovery in equine grass sickness (EGS).

The repository is intended to promote transparency and reproducibility of the computational analyses presented in the manuscript.

Repository contents

The repository includes scripts for:

Data import and preprocessing
Protein filtering and quality control
Principal component analysis (PCA)
Protein abundance visualisation
Differential abundance analysis using limma
Volcano plot generation
Heatmap generation
Gene Ontology enrichment analysis
Reactome pathway enrichment analysis
Venn diagram and protein overlap analyses
Figure generation for the manuscript
Software requirements

Analyses were performed using:

R (version 4.4.2)
Bioconductor (version 3.21)

Key packages include:

tidyverse
limma
ggplot2
ggrepel
ComplexHeatmap
pheatmap
STRINGdb
UpSetR
circlize
viridis
readxl
svglite

Additional package versions are provided within the scripts where applicable.

Data availability

To protect unpublished data and comply with journal submission requirements, raw mass spectrometry files and processed datasets are not included in this repository.

The proteomics data are available through the PRIDE repository under the accession reported in the manuscript.

Reproducibility

Each script is designed to reproduce a specific analysis or figure presented in the manuscript. Users should modify file paths to match the location of their own input files before execution.

Citation

If you use this code in your research, please cite:

Monteiro Moita B., Frampas C., Subbannayya Y., Moulik S., Wells B., Harte T., Proudman C. J., Pinto S. M. Defining the Equine Urinary Proteome: A Reference Baseline for Biomarker Discovery. Journal of Proteome Research (submitted).

Licence

This repository is released under the MIT License.
