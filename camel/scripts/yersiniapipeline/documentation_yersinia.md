# Overview (TODO rewrite after new requirements)
The *Yersinia* pipeline performs complete characterization of *Yersinia enterocolitica* and *Yersinia pseudotuberculosis* isolates.

version: **1.0**

## 1. Coverage check
The workflow starts by checking the coverage of the input FASTQ datasets. Coverage is estimated by dividing the total number of bases by the size of the `NC_GCA_02575835.1` *Y. enterocolitica* reference genome. The total number of bases in the FASTQ file is determined using the `size` function of `seqtk 1.4`.

Datasets with an estimated coverage >=100x are downsampled to ~100x using the `subsample` funcion of `seqtk 1.4`.

## 2. Read trimming
Afterwards, reads are trimmed using `trimmomatic 0.39` with the following options:
```
-phred33
ILLUMINACLIP:NexteraPE-PE.fa:2:30:10
LEADING:10 
TRAILING:10 
SLIDINGWINDOW:4:20 
MINLEN:40
```

Quality reports are generated before and after trimming using `fastqc 0.11.7`.

## 3. Assembly
Processed reads are assembled using `SPAdes 3.15.5` with the following option:
```
--cov-cutoff 'off'
--isolate
```

`QUAST 5.2.0` is then used to check the quality of the resulting assembly with the following options:
```
-r {ref_genome_fasta}
--features {ref_genome_gff3}
--pe1 {forward_reads}
--pe2 {reverse_reads}
```

The completeness of the assembly is checked using `BUSCO 5.5.0` with the following options:
```
--mode genome
--offline
--lineage_dataset bacteria_odb10
```

## 4. Advanced QC

### Kraken 2

The trimmed paired-end reads are checked for contamination using `kraken2 2.1.1` against an in-house database with microbial genomes. The date of the last database update is included in the output report.

### ConFindr

The samples are screened for inter- and intra-species contamination using `ConFindr 0.8.1` with the ribosomal MLST database.

### Quality checks

An overview of the quality checks is provided below. Warnings are included for quality checks that fail but do not stop the pipeline execution. 

| **metric**                             | **warning threshold**  | **fail threshold**   | **description**                                                                                                                                                                                                        |
|----------------------------------------|------------------------|----------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Kraken: contaminants                   | 1.00%                  | 5.00%                | Percentage of reads assigned to species other than *Y. enterocolitica* or *Y. pseudotuberculosis*                                                                                                                      |
| Typing loci detected (%)               | 90%                    | 95%                  | Percentage of cgMLST loci detected (or MLST loci when cgMLST is disabled)                                                                                                                                              |
| Coverage against assembled contigs     | 20x                    | 10x                  | Coverage of the reads mapped to the assembly (determined by QUAST)                                                                                                                                                     |
| Reads mapping to the assembled contigs | 95%                    | 90%                  | Percentage of reads mapping back to the assembly (determined by QUAST)                                                                                                                                                 |
| Total assembly length deviation        | 10%                    | 20%                  | Percent deviation from the expected genome size (determined from the reference genome)                                                                                                                                 |
| ConFindr: number of contaminating SNPs | 10                     | 20                   | Number of SNPs flagged as contaminant by ConFindr                                                                                                                                                                      | 
| Percentage of complete BUSCO genes     | 90%                    | 95%                  | Percentage of complete BUSCO genes identified                                                                                                                                                                          |
| FastQC: Average quality score          | 30                     | 25                   | Checks if the average read quality is above the given threshold.                                                                                                                                                       |
| FastQC: GC-content deviation           | 2.00%                  | 4.00%                | checks if the detected GC content is close enough to the expected GC content for this organism.                                                                                                                        |
| FastQC: Max. N-fraction                | 0.0050                 | 0.0100               | checks if the maximal N fraction at any read position is below the given threshold.                                                                                                                                    |
| FastQC: Per-base sequence content      | 3.00%                  | 6.00%                | checks if the difference between A-T and C-G is below the given threshold at every position. The first 20 and last 5 bases of the reads are skipped, as the peaks there can be caused by the library kit or trimming artifacts. |
| FastQC: Q-score drop                   | 200                    | 150                  | checks whether the average position in the reads where the mean Q-score drops below 30 is above the given threshold.                                                                                                   |
| FastQC: Sequence length distribution   | 66.67%                 | 40.00%               | checks if the median read length of the trimmed reads is below a threshold compared to the mode length of the raw input reads (251).                                                                                   |

**Note:** FastQC metrics are evaluated separately for the forward and reverse reads.

## 5. Gene detection
Gene detection is performed as described in [Bogaerts *et al.*](https://pubmed.ncbi.nlm.nih.gov/30894839/) using an 
updated version of blast (`blast 2.14.0`).
Alternative detection using `kma 1.4.12a` or `srst2 0.2.0` is available by changing the `--detection-method` parameter.

The following databases are available: 

| **name**      | **origin**                                                                         |
|---------------|------------------------------------------------------------------------------------|
| AMRFinderPlus | Antimicrobial resistance genes  from the AMRFinderPlus tool maintained by the NCBI |
| ResFinder     | Antimicrobial resistance genes from the ResFinder tool maintained by DTU           |
| VFDB core     | Databases from the VirulenceFactor Core database                                   | 

## 6. Genomic context investigation
The `MOB-recon` function of `MOB-suite 3.1.4` is used to reconstruct putative plasmids. The contigs assigned to putative
plasmids are cross-checked against the gene detection results for the virulence genes and AMR genes.

## 7. Sequence typing
Sequence typing is performed as described in [Bogaerts *et al.*](https://pubmed.ncbi.nlm.nih.gov/30894839/) with an 
updated version of blast (`blast 2.14.0`). 
Alternative detection using `kma 1.4.12a` or `srst2 0.2.0` is available by changing the `--detection-method` parameter.

The following typing schemes are available:

| **name**        | **origin** |
|-----------------|------------|
| MLST (Achtman)  | Enterobase |
| MLST (McNally)  | Enterobase |
| cgMLST          | BigSDB     |
| cgMLST-YE       | BigSDB     |
| cgMLST-YP       | BigSDB     |
| cgMLST-Yersinia | Enterobase |
| rMLST           | PubMLST    |



