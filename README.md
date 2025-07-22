# **scMicropipe**

A modular Nextflow pipeline for single-cell (scRNA-seq) data analysis, structured in independent sub-workflows for flexibility and reproducibility.

### 1. Pipeline structure

The NanoVI pipeline follows a modular architecture implemented in **Nextflow**.

- The main script `main.nf` acts as an entry point and dispatches the execution based on the `--cmd` parameter.
- Each subcommand (`DGE_matrix`, `Data_processing`, `Celltype_annotation`, `Integration`, `Functionality_analysis`, `Coexpression_modules` ) is implemented as an independent workflow in the `modules/` directory.
- The workflow is further divided into modular Nextflow processes, allowing clear separation of steps and easier maintenance.
- Helper Python scripts used during processing are located in the `bin/` directory.
- A `test/` directory is provided with example data to test the pipeline.


This structure ensures that the pipeline is:

- Easy to extend with new subcommands or modules.
- Highly reproducible and portable through the use of **Nextflow** and **Docker**.
- User-friendly: each subcommand can be run independently by specifying the appropriate parameters.

Below is the current directory structure:

```
├── main.nf # Main Nextflow script: dispatches to sub-workflow modules
├── config/
│ ├── containers.config # Docker/Singularity container configuration
├── modules/ # Modular Nextflow workflows by subcommand
│ ├── DGE_matrix/
│ │ └── main.nf # Generation of digital gene expression matrices
│ ├── Data_processing/
│ │ └── main.nf # scRNA-seq data filtering and processing
│ ├── Celltype_annotation/
│ │ └── main.nf # Cell type annotation and classification
│ ├── Integration/
│ │ └── main.nf # Dataset integration and batch effect correction
│ ├── Functionality_analysis/
│ │ └── main.nf # Functional analysis (e.g., GSEA, pathway analysis)
│ └── Coexpression_modules/
│ └── main.nf # Construction and analysis of coexpression modules
├── bin/ # Helper scripts in Python/R
├── test/ 
├── README.md 
└── nextflow.config # Global configuration and default parameters
```
