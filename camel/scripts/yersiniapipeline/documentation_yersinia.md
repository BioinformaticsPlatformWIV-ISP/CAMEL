# Overview
The *Yersinia* pipeline performs complete characterization of *Yersinia enterocolitica* and *Yersinia pseudotuberculosis* isolates.

version: **0.1**

## 1. Coverage check
The workflow starts by checking the coverage of the input FASTQ datasets. Coverage is estimated by dividing the total number of bases by the average size of the the `T00469` *Y. enterocolitica* and the `T00195` *Y. pseudotuberculosis*. The total number of bases in the FASTQ file is determined using the `size` function of `seqtk 1.4`.

Datasets with an estimated coverage >=100x are downsampled to ~100x using the `subsample` funcion of `seqtk 1.4`.

## 2. Read trimming
TODO