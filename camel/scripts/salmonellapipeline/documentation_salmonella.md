# Overview
The *Salmonella* pipeline performs complete characterization of Salmonella isolates.

Version: **0.3**

# Components

**Note:** If the input type is `fasta`, pre-processing steps 2 to 4 are skipped.

## 1. Human read removal (optional) 

If enabled, human reads are removed using the NCBI Human Read Removal Tool (HRRT) 2.2.1.
The tool is executed with default options.

## 2. Coverage check
The workflow starts by checking the coverage of the input FASTQ datasets. 
Coverage is estimated by dividing the total number of bases by the size of the `NC_003197.2` *S. enterica* 
reference genome. The total number of bases in the FASTQ files is determined using the `size` function of 
`seqtk 1.4`.

Datasets with an estimated coverage >=100x are downsampled to ~100x using the `subsample` function of `seqtk 1.4`.

## 3. Read trimming

Read trimming is performed using `fastp 0.23.4` (default) or `trimmomatic 0.39`.

For `fastp` the following options are used:
```
--compression 4
--detect_adapter_for_pe
--cut_front
--cut_front_window_size 1
--cut_front_mean_quality 10
--cut_tail
--cut_tail_window_size 1
--cut_tail_mean_quality 10
--cut_right
--cut_right_window_size 4
--cut_right_mean_quality 20
--length_required 40
```

For `trimmomatic` the following options are used:
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

The trimmed paired-end reads or contigs are checked for contamination using `kraken2 2.1.1` against an in-house database
with microbial genomes. The date of the last database update is included in the output report.

### ConFindr

The samples are screened for inter- and intra-species contamination using `ConFindr 0.8.1` with the ribosomal MLST 
database.

### Quality checks

An overview of the quality checks is provided below. Warnings are included for quality checks that fail but do not stop 
the pipeline execution. 

| **metric**                             | **warning threshold**  | **fail threshold**   | **description**                                                                                                                                                                                                                 |
|----------------------------------------|------------------------|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Kraken: contaminants                   | 1.00%                  | 5.00%                | Percentage of reads / contigs assigned to species other than *N. meningitidis*                                                                                                                                                  |
| Typing loci detected (%)               | 90%                    | 95%                  | Percentage of cgMLST loci detected (or MLST loci when cgMLST is disabled)                                                                                                                                                       |
| Coverage against assembled contigs     | 20x                    | 10x                  | Coverage of the reads mapped to the assembly (determined by QUAST)                                                                                                                                                              |
| Reads mapping to the assembled contigs | 95%                    | 90%                  | Percentage of reads mapping back to the assembly (determined by QUAST)                                                                                                                                                          |
| Total assembly length deviation        | 10%                    | 20%                  | Percent deviation from the expected genome size (determined from the reference genome)                                                                                                                                          |
| ConFindr: number of contaminating SNPs | 10                     | 20                   | Number of SNPs flagged as contaminant by ConFindr                                                                                                                                                                               | 
| Percentage of complete BUSCO genes     | 90%                    | 95%                  | Percentage of complete BUSCO genes identified                                                                                                                                                                                   |
| FastQC: Average quality score          | 30                     | 25                   | Checks if the average read quality is above the given threshold.                                                                                                                                                                |
| FastQC: GC-content deviation           | 2.00%                  | 4.00%                | checks if the detected GC content is close enough to the expected GC content for this organism (50.50%).                                                                                                                        |
| FastQC: Max. N-fraction                | 0.0050                 | 0.0100               | checks if the maximal N fraction at any read position is below the given threshold.                                                                                                                                             |
| FastQC: Per-base sequence content      | 3.00%                  | 6.00%                | checks if the difference between A-T and C-G is below the given threshold at every position. The first 20 and last 5 bases of the reads are skipped, as the peaks there can be caused by the library kit or trimming artifacts. |
| FastQC: Q-score drop                   | 200                    | 150                  | checks whether the average position in the reads where the mean Q-score drops below 30 is above the given threshold.                                                                                                            |
| FastQC: Sequence length distribution   | 66.67%                 | 40.00%               | checks if the median read length of the trimmed reads is below a threshold compared to the mode length of the raw input reads (251).                                                                                            |

**Note:** FastQC metrics are evaluated separately for the forward and reverse reads.
**Note:** *Escherichia* is not considered a contaminant for the Kraken 2 QC check.

The QC checks enabled for the supported input types are listed in the table below.

| **metric**                             | **illumina** | **fasta** |
|----------------------------------------|--------------|-----------|
| Kraken: contaminants                   | Yes          | Yes       | 
| Typing loci detected (%)               | Yes          | Yes       | 
| Coverage against assembled contigs     | Yes          | No        | 
| Reads mapping to the assembled contigs | Yes          | No        | 
| Total assembly length deviation        | Yes          | Yes       | 
| ConFindr: number of contaminating SNPs | Yes          | No        |  
| Percentage of complete BUSCO genes     | Yes          | Yes       | 
| FastQC: Average quality score          | Yes          | No        | 
| FastQC: GC-content deviation           | Yes          | No        | 
| FastQC: Max. N-fraction                | Yes          | No        | 
| FastQC: Per-base sequence content      | Yes          | No        | 
| FastQC: Q-score drop                   | Yes          | No        | 
| FastQC: Sequence length distribution   | Yes          | No        |


## 6. *Salmonella* serotyping

### SISTR

`SISTR (v1.1.1)` (*Salmonella* In Silico Typing Resource) is able to rapidly type and subtype *Salmonella* genome 
assemblies. See [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147101) for more details.

### Seqsero

`Seqsero (v1.2.1)` determines *Salmonella* serotypes based on raw sequencing reads or genome assemblies. See 
[here](https://journals.asm.org/doi/10.1128/jcm.00323-15) for more details.

Based on the input type, both modes (reads + assembly) or only one mode (assembly) is carried out.

### Mykrobe

`Mykrobe (v0.13.0)` has a genotyping scheme for *Salmonella typhi*. It efficiently identifies the genotypes, AMR and 
plasmid markers from WGS data or assemblies. See [here](https://github.com/Mykrobe-tools/mykrobe?tab=readme-ov-file) 
for more details.

## 7. Antimicrobial resistance (AMR) detection

### AbriTAMR

`AbriTAMR 1.0.14` utilises NCBI’s AMRFinderPlus to detect genes associated with AMR. A validated antibiogram is 
reported for *Salmonella*.

The database version is indicated in the output report and summary output file.


### ResFinder4

`ResFinder4 4.4.2` is used with the following options to detect genes and mutation associated with AMR:

```
--min_cov 0.6
--acquired
--threshold 0.9
--species "Escherichia coli"
```

The database version is indicated in the output report and summary output file.

## 8. Pathogenicity island determination

### SPIFinder
`SPIFinder 1.4.12a` identifies *Salmonella* Pathogenicity Islands in sequencing data or assemblies.

Based on the input type, both modes (reads + assembly) or only one mode (assembly) is carried out.

The database version is indicated in the output report and summary output file.

## 9. Virulence characterization

### Gene detection

Gene detection is performed as described in [Bogaerts *et al.*](https://pubmed.ncbi.nlm.nih.gov/30894839/) using an 
updated version of blast (`blast 2.14.0`). Alternative detection using `kma 1.4.12a` or `srst2 0.2.0` is available by 
changing the `--detection-method` parameter.

The following databases are available: 

| **name**  | **origin**                                       |
|-----------|--------------------------------------------------|
| VFDB core | Databases from the VirulenceFactor Core database |


## 10. Plasmid characterization

### MOB-suite & genomic context

The `MOB-recon` function of `MOB-suite 3.1.4` is used to reconstruct putative plasmids. The contigs assigned to putative
plasmids are cross-checked against the gene detection results for the virulence genes and AMR genes.

## 11. Sequence typing

Sequence typing is performed as described in [Bogaerts *et al.*](https://pubmed.ncbi.nlm.nih.gov/30894839/) with an 
updated version of blast (`blast 2.14.0`). 
Alternative detection using `kma 1.4.12a` or `srst2 0.2.0` is available by changing the `--detection-method` parameter.

The following typing schemes are available:

| **name**               | **origin**     |
|------------------------|----------------|
| Classic MLST (Warwick) | EnteroBase     |
| cgMLST                 | EnteroBase     |
| rMLST                  | PubMLST        |
