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
