# Mock pipeline

## Description

The `mock pipeline` is a minimal pipeline used for testing the core components of the CAMEL pathogen characterization 
pipelines. The mock pipeline supports the following input types: `illumina`, `ont`, and `hybrid`.

## Included assays

- Downsampling
- Read trimming
- *de novo* assembly
- Contamination check (Kraken2+, ConFindr)
- Advanced QC checks

*Note:* All assays are compatible with the supported input types.

## Test dataset

The `mock pipeline` testdata set consists of at ~10 kb section of the E. coli genome, contaning the X AMR gene.
