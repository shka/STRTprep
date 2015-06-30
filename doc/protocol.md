# Protocol of preprocessing and analysis

## Table of contents

1. [Overview of a typical protocol](#typical-protocol)
2. Protocol details
  - [Preparation of project folder](#preparation-of-project-folder)
  - [Preparation of indexes](#preparation-of-indexes)
  - [Design of experiments](#design-of-experiments)

## Typical protocol

- **Step 1: Installation.** Prepare [your project folder](#preparation-of-project-folder) with [indexes of reference genome and transcriptome](#preparation-of-indexes), like as
```bash
git clone https://github.com/shka/STRTprep.git STRTprep2.test
cd STRTprep2.test
bin/install.sh
. bin/setup.sh
bin/index_hg19_ercc92_ynbA_u13369.sh
bin/index_refGene.sh hg19 src/ebwt/hg19_ercc92_ynbA_u13369/ref
```

> You can pull the bug fixes and further more new features by `git pull` within the project folder.

> `[v3dev]` When you try new features in the version 3 beta, create your project folder by `git clone -b v3dev https://github.com/shka/STRTprep.git`. Although it might contain bugs, but it will update frequently.

- **Step 2: [Design of your experiments](#design-of-experiments).** Edit `conf.yaml` and `src/samples.csv`.

- **Step 3: Preprocessing and quality-check.** Run as follows, then check `out/byGene/samples.xls`; see also [a section in "Interpretation of the results"](result.md#quality-check).
```bash
rake -m -j `gnproc` qc > qc.log 2>&1
```
- **Step 4: Differential expression tests.** Run as follows, and check `out/byGene/diffexp.xls` and `out/byGene/plugin_*`, then analyze more further or go back to the step 2; see also [a section in "Interpretation of the results"](result.md#differential-expression-tests).
```bash
rake gene > gene.log 2>&1
```

> You need `. bin/setup.sh` again before the steps 3 and 4, when you logout/exit before the step.

## Protocol details

### Preparation of project folder

STRTprep consists one folder containing configurations and scripts. It creates output files also in the folder after the execution. Therefore, you need to prepare one STRTprep folder per project. Here we use `STRTprep2.test` for example.

```
git clone https://github.com/shka/STRTprep.git STRTprep2.test
cd STRTprep2.test
```

> Do not contain space, special characters, which you need escape, in the pathname of your project folder.

There is a script for installation of the additionally required softwares. You need to run it only once when you create the project folder.

```
bin/install.sh
```

### Preparation of indexes

You need to prepare at least two indexes to align the STRT reads with (i) reference genome and spike-in, and (ii) reference transcriptome. Moreover, you may need to prepare one more index to exclude phyX sequences in the raw output reads, when the sequencing facility added. The following is an example to build the indexes on hg19 and RefSeq.

```bash
. bin/setup.sh
bin/index_phyX.sh
bin/index_hg19_ercc92_ynbA_u13369.sh
bin/index_refGene.sh hg19 src/ebwt/hg19_ercc92_ynbA_u13369/ref
```

> These indexes are for bowtie version 1 and tophat, so you can build index as you want.

#### Index of reference genome and spike-in

Currently there are five scripts to prepare index of the reference genome and the spike-ins.

- `bin/index_hg19_ercc92_ynbA_u13369.sh` creates index of (i) human hg19 genome, (ii) ERCC92 spike-ins, (iii) ynbA spike-in, and (iv) U13369 human ribosomal DNA complete repeating unit, at `src/ebwt/hg19_ercc92_ynbA_u13369/ref`.
- `bin/index_mm9_ercc92_ynbA_bk000964.sh` creates index of (i) mouse mm9 genome, (ii) ERCCC92 spike-ins, (iii) ynbA spike-in, and (iv) BK000964 mouse ribosomal DNA complete repeating unit, at `src/ebwt/mm9_ercc92_ynbA_bk000964/ref`.
- `bin/index_canFam3_ercc92_ynbA.sh` creates index of (i) dog canFam3 genome, (ii) ERCCC92 spike-ins, and (iii) ynbA spike-in, at `src/ebwt/canFam3_ercc92_ynbA/ref`.
- `bin/index_susScr3_ercc92_ynbA.sh` creates index of (i) pig susScr3 genome, (ii) ERCCC92 spike-ins, and (iii) ynbA spike-in, at `src/ebwt/susScr3_ercc92_ynbA/ref`.

#### Index of reference transcriptome

Currently there are three scripts to prepare index of the reference transcriptome. It requires two options, (i) genome version, and (ii) genome+spike-in index; for example `bin/index_refGene.sh hg19 src/ebwt/hg19_ercc92_ynbA_13369`.

- `bin/index_ensGene.sh` creates index based on ENSEMBL at `src/ebwt/ver_ensGene/ref`.
- `bin/index_knownGene.sh` creates index based on UCSC known genes at `src/ebwt/ver_knownGene/ref`.
- `bin/index_refGene.sh` creates index based on RefSeq at `src/ebwt/ver_refGene/ref`.

While update of the genome sequence is rarely, the reference transcriptome updates very frequently, twice a year as usual.

#### Index of PhyX

You need to build the third index when you need to exclude PhyX control sequences.

- `bin/index_phyX.sh` creates index at `src/ebwt/phyX/ref`

### Design of experiments

There are two configuration files in your STRTprep project folder, (i) `conf.yaml`, and (ii) `src/samples.csv`.

#### [`conf.yaml`](https://github.com/shka/STRTprep/blob/master/conf.yaml); STRT library information

This is information about your libraries, written by [YAML](http://www.yaml.org/start.html) format. You can edit it by any text editor.

- Keys `UMI`, `BARCODE`, `GAP`, `CDNA` (lines 2-5): Lengths of UMI, barcode, gap and cDNA in each STRT read.
- `LAYOUT` (line 6): Barcode layout file with barcode+gap sequences and well names; although there are two predefined mapping for 48-plex STRT with 6 bp barcodes (and 3 bp gap), (i) `src/barcodes.old.txt`, and (ii) `src/barcodes.txt`, but you can specify your own layout.
- `PHYX` (line 7): PhyX index; you need to specify if you need to exclude phyX sequences from the raw sequences.
- `GENOMESPIKERIBO` (line 8): Genome index
- `TRANSCRIPT` (line 9): Transcriptome index
- `FLUCTUATION` (line 10): Fluctuation p-value threshold for plugins
- `DIFFEXP` (line 11): Differential expression q-value threshold for plugins
- `TEST1` (line 13): `TEST1` itself is name of the first library; you can give any name as you want.
- Line 16 & 17 after `FASTQS`: Raw sequences (fastq format) for the first library; you must give at least one, and you can give more sequences.
- Line 18-22, 23-27, and 28-32: Configurations of the second to fourth libraries; you must give at least one, and you can give more libraries.

#### [`src/samples.csv`](https://github.com/shka/STRTprep/blob/master/src/samples.csv); sample information and study design

This is sample information and the study design, written by [CSV](https://tools.ietf.org/html/rfc4180) format. You can edit it by various spreadsheet editor, for example Microsoft Excel.

> Microsoft Office uses ";" (semicolon) instead of "," (colon) as column separater, when you your computer with European environment. In case of OSX, please change the primary "Language" of your OSX at System Preferences to, for example, US.

- Column `LIBRARY`: Library name, given in `conf.yaml`
- `WELL`: Well name, based on the barcode layout, specified in `conf.yaml`
- `NAME`: Sample name; "NA" in case of empty or ignore; " " (space) characters at the end of the name values are prohibited.
- `CLASS.n` (n=0, 1, ...; optional): Sample class number for differential expression test; 1 or 2 for two class comparison; 1, 2, 3, … for multiclass comparison; "NA" in case of empty or ignore
- `BLOCK.n` (n=0, 1, ...; optional): Permutation block number for the differential expression test; 1, 2, … ; identical number is assigned to a pair for the paired comparison, to a same gender, and/or to a library for canceling of library bias, for example; "NA" in case of empty or ignore
- `FORCE_APPROVAL` (optional): If `TRUE` as for inter-library control samples etc., we can ignore the automatic outlier check for the samples. The other samples to be checked must be `FALSE`, when you use this column.

> `CLASS.n`, and `BLOCK.n` for design of the differential expression tests are used by SAMstrt [[Katayama et al. 2013](http://www.ncbi.nlm.nih.gov/pubmed/?term=23995393)]. See also SAMseq [[Li and Tibshirani 2013](http://www.ncbi.nlm.nih.gov/pubmed?term=22127579)] and the manual in the [official web page](http://statweb.stanford.edu/~tibs/SAM/) to understand the statistical background and the usage.

The `src/samples.csv` file can contain any other columns about the sample/layout information.
