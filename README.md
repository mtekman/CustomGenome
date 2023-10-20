# CustomGenome for R / Bioconductor

A package to append custom FASTA sequences to your genome and
annotation files for downstream mapping.

CustomGenome generates a new Ensembl reference genome sequence and
annotation files with custom FASTA sequences (e.g.  GFP, Tomato) and
produces RNA-seq count matrices with them. Genome reference files are
fetched and validated directly from Ensembl, user-sequences are then
appended to both, and then the new reference files are written to
disk. Downstream analysis is then facilitated via Rsubread function
wrappers to index the new reference, align user FASTA sequences to it,
and then perform quantification for desired features (e.g. exon, utr,
etc.)


## Installation

Install the package from Bioconductor or Gitlab, ensuring correct
dependencies.

#### From Bioconductor

```r
BiocManager::install("CustomGenome")
```

#### From Bioconda

```bash
micromamba install -c bioconda -c conda-forge bioconductor-customgenome
```

#### From Gitlab

```r
library(remotes)
remotes::install_github("mtekman/CustomGenome",
                        repos = BiocManager::repositories())
```


## Getting Started

Load the library

```r
library(CustomGenome)
```

Follow the vignette to learn more about how to get started with this package.