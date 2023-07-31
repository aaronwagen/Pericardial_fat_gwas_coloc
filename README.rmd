---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# Background
This repository contains scripts used to perform a colocalisation analysis using the pericardial fat GWAS and tissue eQTLs from GTEx v8. 
# Citation
If you use any of the code or data from this repository, please cite our [paper](). Further, if you use any of the software used within this repository (e.g. `coloc`, `colochelpR`, etc.) please make sure to cite the software appropriately. 

# License
The code in this repository is released under an MIT license. This repository is distributed in the hope that it will be useful to the wider community, but without any warranty of any kind. Please see the [LICENSE](LICENSE) file for more details. 

# Code contents

The workflow and analysis, including utilisation of any scripts, is described in analyse_pericard_fat_gwas_coloc.Rmd.

It can be view interactively at: [https://aaronwagen.github.io/pat_gwas_gtex_coloc/](https://aaronwagen.github.io/pat_gwas_gtex_coloc/)

Within this repository you will otherwise find:

| Directory | Description |
| --------- | --------------------------------------------------------------------------- |
| [docs](docs) | Contains all `.Rmd`s and their corresponding `.html`s describing analyses performed for this project. |
| [logs](logs) | For any scripts that were run outside of an `.Rmd` (e.g. scripts from the [scripts](scripts) directory), a log file was recorded and can be accessed here. |
| [manuscript](manuscript) | Figures and tables produced for the manuscript. |
| [processed_data](processed_data) | Results from all analyses. |
| [scripts](scripts) | Contains analysis scripts. |
