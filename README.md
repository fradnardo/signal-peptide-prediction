# Signal Peptide Prediction

End-to-end bioinformatics and machine-learning project for the prediction
of secretory signal peptides from protein sequences.

## Overview

Secretory signal peptides are short N-terminal sequences that direct
proteins toward the secretory pathway.

This project investigates computational approaches for signal-peptide
prediction and compares:

- a position-specific statistical approach based on the von Heijne method;
- a feature-based Support Vector Machine classifier.

The workflow includes:

1. biological dataset construction;
2. data preprocessing and exploratory analysis;
3. feature engineering;
4. model development;
5. cross-validation and independent testing;
6. performance comparison;
7. biological interpretation of the results.

## Biological background

A typical signal peptide contains three main regions:

- **n-region**: enriched in positively charged residues;
- **h-region**: hydrophobic core;
- **c-region**: contains the cleavage-site region.

The main computational task is to distinguish proteins containing
secretory signal peptides from proteins without them, while avoiding
confusion with N-terminal transmembrane helices.

## UniProt data collection

The collector starts from reviewed eukaryotic UniProtKB entries. Both classes
require a non-fragment sequence of at least 40 residues and evidence at the
protein level.

- **Positives** need a Signal feature supported by `ECO:0000269`, an exact
  cleavage site, and a signal peptide of at least 14 residues. The
  `cleavage_after` column is the one-based index of the last signal peptide
  residue. A separate audit row counts exclusions for multiple eligible cuts,
  since the output needs one coordinate per protein.
- **Negatives** need zero Signal annotations and at least one experimentally
  supported location in the cytosol/cytoplasm, nucleus, mitochondrion,
  plastid (including chloroplast), peroxisome, or cell membrane. The location
  must be annotated for the entry's reference sequence. The `early_tm_90`
  column flags a transmembrane helix starting at residues 1–90.

A protein with only inferred Signal annotations stays outside both classes.
This keeps the negative set from treating a lack of experimental signal
evidence as evidence that the signal is absent.

The collector uses Requests. Install the project dependencies from the
repository root with `.venv/bin/python -m pip install -e .`.

Run a two-page preview from the repository root using a new output directory.

```sh
.venv/bin/python src/signal_peptide_prediction/collect_uniprot.py \
  --max-pages 2 --output-dir data/processed/my_uniprot_preview
```

For the full collection, omit `--max-pages` and choose a new output directory.
To walk through saved entry JSON files offline, use
`--input-json-dir data/raw/pilot` in place of `--max-pages 2`.

Both classes have a TSV joined by accession to a FASTA of full sequences.
`filter_counts.tsv` shows the input, kept, and dropped counts at each step.
Both class branches start from the last shared count. `run.json` records the
query, release, and collection size for live API runs.

## Planned methods

- von Heijne position-specific scoring method
- Support Vector Machine
- optional additional predictive approaches

## Repository structure

-  — raw and processed biological data
-  — exploratory analysis
-  — reusable Python source code
-  — executable workflows
-  — automated tests
-  — figures and tables
-  — scientific project report

## Status

Work in progress.
