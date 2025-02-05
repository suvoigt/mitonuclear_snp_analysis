# README: Mitochondrial and Nuclear Genome SNP Analysis Pipelines

## Overview

These scripts were created as part of the study *Testing for Age- and Sex-Specific Mitonuclear Epistasis in Drosophila* by Garlovsky et al., 2025.

This repository includes two scripts for SNP analysis:

- **Mitochondrial SNP Analysis** (`mito_mt_seq_snps_diff.sh`): Performs SNP calling, filtering, and genetic differentiation analysis on mitochondrial sequencing data.
- **Nuclear Genome SNP Analysis** (`mito_nuc_seq_snps_var_diff.sh`): Conducts SNP calling, variation analysis, and genetic differentiation on nuclear genome sequencing data.

## Requirements

### Software Dependencies:

- `samtools`
- `PoolSNP` ([https://github.com/capoony/PoolSNP](https://github.com/capoony/PoolSNP))
- `Popoolation2` ([https://sourceforge.net/projects/popoolation2/](https://sourceforge.net/projects/popoolation2/))
- `SNPeff` (doi: 10.4161/fly.19695)
- `RepeatMasker`
- `Python2` (with scripts from `DrosEU_pipeline`)
- `Java`
- `Perl`
- `awk`, `sed`, `gzip`

## References

- **PoolSNP:** doi:10.1093/molbev/msaa120
- **SNPeff:** doi: 10.4161/fly.19695
- **Popoolation2:** [https://sourceforge.net/projects/popoolation2/](https://sourceforge.net/projects/popoolation2/)
- **DrosEU\_pipeline:** [https://github.com/capoony/DrosEU\_pipeline](https://github.com/capoony/DrosEU_pipeline)
- **RepeatMasker:** [https://www.repeatmasker.org/](https://www.repeatmasker.org/)

