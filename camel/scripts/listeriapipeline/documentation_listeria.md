# Overview
The *Listeria* pipeline performs complete characterization of *Listeria monocytogenes* isolates.

Version: **1.3**

# Components

## 1. Human read removal (optional) 

If enabled, human reads are removed using the NCBI Human Read Removal Tool (HRRT) 2.2.1.
The tool is executed with default options.

## 2. Coverage check
The workflow starts by checking the coverage of the input FASTQ datasets. 
Coverage is estimated by dividing the total number of bases by the size of the `NC_003210.1` *L. monocytogenes* 
EGD-e reference genome. The total number of bases in the FASTQ files is determined using the `size` function of 
`seqtk 1.4`.

Datasets with an estimated coverage >=100x are downsampled to ~100x using the `subsample` function of `seqtk 1.4`.

## 3. Read trimming

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

## 4. Assembly

Processed reads are assembled using `SPAdes 3.15.5` with the following options:
```
--cov-cutoff 10
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

## 5. Advanced QC

### Kraken 2

The trimmed paired-end reads are checked for contamination using `kraken2 2.1.1` against an in-house database with 
microbial genomes. The date of the last database update is included in the output report.

### ConFindr

The samples are screened for inter- and intra-species contamination using `ConFindr 0.8.1` with the ribosomal MLST 
database.

### Quality checks

An overview of the quality checks is provided below. Warnings are included for quality checks that fail but do not stop 
the pipeline execution. 

| **metric**                             | **warning threshold**  | **fail threshold**   | **description**                                                                                                                                                                                                                 |
|----------------------------------------|------------------------|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Kraken: contaminants                   | 1.00%                  | 5.00%                | Percentage of reads assigned to species other than *L. monocytogenes*                                                                                                                                                           |
| Typing loci detected (%)               | 90%                    | 95%                  | Percentage of cgMLST loci detected (or MLST loci when cgMLST is disabled)                                                                                                                                                       |
| Coverage against assembled contigs     | 20x                    | 10x                  | Coverage of the reads mapped to the assembly (determined by QUAST)                                                                                                                                                              |
| Reads mapping to the assembled contigs | 95%                    | 90%                  | Percentage of reads mapping back to the assembly (determined by QUAST)                                                                                                                                                          |
| Total assembly length deviation        | 10%                    | 20%                  | Percent deviation from the expected genome size (determined from the reference genome)                                                                                                                                          |
| ConFindr: number of contaminating SNPs | 10                     | 20                   | Number of SNPs flagged as contaminant by ConFindr                                                                                                                                                                               | 
| Percentage of complete BUSCO genes     | 90%                    | 95%                  | Percentage of complete BUSCO genes identified                                                                                                                                                                                   |
| FastQC: Average quality score          | 30                     | 25                   | Checks if the average read quality is above the given threshold.                                                                                                                                                                |
| FastQC: GC-content deviation           | 2.00%                  | 4.00%                | checks if the detected GC content is close enough to the expected GC content for this organism (38.00%).                                                                                                                        |
| FastQC: Max. N-fraction                | 0.0050                 | 0.0100               | checks if the maximal N fraction at any read position is below the given threshold.                                                                                                                                             |
| FastQC: Per-base sequence content      | 3.00%                  | 6.00%                | checks if the difference between A-T and C-G is below the given threshold at every position. The first 20 and last 5 bases of the reads are skipped, as the peaks there can be caused by the library kit or trimming artifacts. |
| FastQC: Q-score drop                   | 200                    | 150                  | checks whether the average position in the reads where the mean Q-score drops below 30 is above the given threshold.                                                                                                            |
| FastQC: Sequence length distribution   | 66.67%                 | 40.00%               | checks if the median read length of the trimmed reads is below a threshold compared to the mode length of the raw input reads (251).                                                                                            |

**Note:** FastQC metrics are evaluated separately for the forward and reverse reads.

## 6. Gene detection

Gene detection is performed as described in [Bogaerts *et al.*](https://pubmed.ncbi.nlm.nih.gov/30894839/) using an 
updated version of blast (`blast 2.14.0`).
Alternative detection using `kma 1.4.12a` or `srst2 0.2.0` is available by changing the `--detection-method` parameter.

The following databases are available: 

| **name**       | **origin**                                                               |
|----------------|--------------------------------------------------------------------------|
| NDARO          | Antimicrobial resistance genes from the NCBI NDARO database              |
| PlasmidFinder  | Plasmid replicons from the PlasmidFinder tool maintained by DTU          |
| ResFinder      | Antimicrobial resistance genes from the ResFinder tool maintained by DTU |
| VFDB core      | Databases from the VirulenceFactor Core database                         | 
| VirulenceFiner | Virulence genes from the VirulenceFinder tool maintained by DTU          |

## 7. Sequence typing

Sequence typing is performed as described in [Bogaerts *et al.*](https://pubmed.ncbi.nlm.nih.gov/30894839/) with an 
updated version of blast (`blast 2.14.0`). 
Alternative detection using `kma 1.4.12a` or `srst2 0.2.0` is available by changing the `--detection-method` parameter.

The following typing schemes are available:

| **name**             | **origin**                |
|----------------------|---------------------------|
| rMLST                | BIGSdb (Institut Pasteur) |
| Classic MLST         | BIGSdb (Institut Pasteur) |
| cgMLST               | BIGSdb (Institut Pasteur) |
| PCR-serogroup        | BIGSdb (Institut Pasteur) |
| Species confirmation | BIGSdb (Institut Pasteur) |
| AMR typing           | BIGSdb (Institut Pasteur) |
| Virulence typing     | BIGSdb (Institut Pasteur) |
