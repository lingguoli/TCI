# TCI - Tissue Contribution Index

A bioinformatics toolkit for converting various file to **Tissue Contribution Index**.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [fq2tci](#fq2tci)
  - [bam2tci](#bam2tci)
  - [bed2tci](#bed2tci)
  - [finaleDB2tci](#finaledb2tci)
- [Requirements](#requirements)


## Features

- Convert FASTQ, BAM, BED, and FinaleDB files to TCI format
- Comprehensive parameter validation
- Automatic output directory creation
- Detailed help documentation for all commands

## Installation

```bash
git clone https://github.com/lingguoli/TCI.git
cd TCI
chmod +x TCI
```

## Usage
```bash
./TCI <command> [options]
```

### fq2tci
Convert FASTQ files to TCI.csv

Options:
```
Option	Description
-i, --input1	Input FASTQ1 file (required)
-I, --input2	Input FASTQ2 file (required unless -s)
-a, --adapter_sequence	Adapter sequence (required if -s)
-o, --output_directory	Output directory (default: current directory)
-r, --reference	    Reference genome (default: ./Ref/hg38.fa)
-s, --single_end	Single end mode (default: false)
-h, --help	        Show help message
```

Example:

```bash
# paired end
./TCI fq2tci -i ./demo/test_1.fq.gz -I ./demo/test_2.fq.gz -o ./demo/fq2tci_PE
# single end
./TCI fq2tci -i ./demo/test_1.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o ./demo/fq2tci_SE -s
```

### bam2tci
Convert BAM files to TCI.csv

Options:
```
Option	Description
  -i, --bam                Input BAM file    (required, both .bam and .bam.bai must exist)
  -o, --output_directory   Output directory  (default: current directory)
  -r, --reference          Reference genome  (default: ./Ref/hg38.fa)
  -s, --single_end         Single end input  (default: false)
  -h, --help               Show help message

```

### bed2tci
Convert BED files to TCI.csv

Options:
```
Option Description
  -i, --bed                Input BED file    (required)
  -o, --output_direc       Output directory  (default: current directory)
  -r, --reference          Reference genome  (default: ./Ref/hg38.fa)
  -s, --single_end         Single end input  (default: false)
  -h, --help               Show help message
```
### finaleDB2tci
Convert FinaleDB to TCI.csv

```
Option Description
  -i, --input              Input FinaleDB .tsv.bgz (required)
  -o, --output_directory   Output directory 
  -h, --help               Show help message
```
Example:

```bash
./TCI finaleDB2tci -i ./demo/EE87874.hg38.frag.tsv.bgz -o ./demo/finaleDB2tci
```

## Requirements
```
fastp 0.23.4
minimap2 2.28
samblaster 0.1.26
samtools 1.21
bedtools 2.25.0
python3 3.12.0
pysam 0.22.1
pandas 2.2.2
numpy 1.26.4
matplotlib 3.9.2
pybedtools 0.10.0
```
