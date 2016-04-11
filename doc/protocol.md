# Protocol of preprocessing and analysis

## Table of contents

1. [Overview of a typical protocol](#typical-protocol)
2. Standard protocol details
  - [Preparation of project folder](#preparation-of-project-folder)
  - [Preparation of indexes](#preparation-of-indexes)
  - [Design of experiments](#design-of-experiments)
3. Special protocols
  - [Overexpression study](#overexpression-study)
  - [Rescue of overlapping gene](#rescue-of-overlapping-gene)
  - [STRT by Fluidigm C1](#strt-by-fluidigm-c1)

## Typical protocol

- **Step 1: Installation.** Prepare [your project folder](#preparation-of-project-folder) with [indexes of reference genome and transcriptome](#preparation-of-indexes), like as
```bash
git clone -b v3dev https://github.com/shka/STRTprep.git STRTprep3.test
cd STRTprep3.test
bin/install.sh
. bin/setup.sh
bin/index_hg19_ercc92_ynbA_u13369.sh
bin/index_refGene.sh hg19 src/ebwt/hg19_ercc92_ynbA_u13369/ref
```

> You can pull the bug fixes and further more new features by `git pull` within the project folder.

- **Step 2: [Design of your experiments](#design-of-experiments).** Edit `src/conf.yaml` and `src/samples.csv`.

> At first, you should specify only the three required columns for `src/samples.csv` and go through until step 4. Then you repeat from the step 2 again.

- **Step 3: Preprocessing and quality-check.** Run as follows, then check `out/byGene/samples.xls`; see also "[Quality check](result.md#quality-check)".
```bash
rake -m -j `gnproc` qc > qc.log 2>&1
```

> Since many intermediate tasks in the step 3 are skipped, it would finish quickly in the second round or later.

> You need `. bin/setup.sh` again before the steps 3, 4 and 5, when you restert after logout/exit.

- **Step 4: Differential expression tests.** Run as follows, and check `out/byGene/diffexp.xls` and reports by plugins, then go back to the step 2 if you add or change the test designs; see also "[Differential expression tests](result.md#differential-expression-tests)".
```bash
rake gene > gene.log 2>&1
```

- **Step 5: High resolution analysis.** Run as follows, and check `out/byTFE/diffexp.xls` and reports by plugins, then enjoy further downstream analysis!
```bash
rake > tfe.log 2>&1
```

## Protocol details

### Preparation of project folder

STRTprep consists one folder containing configurations and scripts. It creates output files also in the folder after the execution. Therefore, you need to prepare one STRTprep folder per project. Here we use `STRTprep3.test` for example.

```
git clone -b v3dev https://github.com/shka/STRTprep.git STRTprep3.test
cd STRTprep3.test
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

#### Index of reference genome and spike-in

This is bowtie version 1 index containing reference genomes, spike-ins, and any other sequences of your interest. There are several scripts to prepare index for major species.

Script | Genome | Spike-ins | Ribosomal DNA unit | Created index location
-------|--------|-----------|--------------------|---------
`bin/index_hg19_ercc92_ynbA_u13369.sh` | hg19 (human) | ERCC92 & ynbA | U13369 | `src/ebwt/hg19_ercc92_ynbA_u13369/ref`
`bin/index_hg38_ercc92_ynbA_u13369.sh` | hg38 (human) | ERCC92 & ynbA | U13369 | `src/ebwt/hg38_ercc92_ynbA_u13369/ref`
`bin/index_mm9_ercc92_ynbA_bk000964.sh` | mm9 (mouse) | ERCC92 & ynbA | BK000964 | `src/ebwt/mm9_ercc92_ynbA_u13369/ref`
`bin/index_canFam3_ercc92_ynbA.sh` | canFam3 (dog) | ERCC92 & ynbA | | `src/ebwt/canFam3_ercc92_ynbA/ref`
`bin/index_susScr3_ercc92_ynbA.sh` | susScr3 (pig) | ERCC92 & ynbA | | `src/ebwt/susScr3_ercc92_ynbA/ref`
`bin/index_danRer7_ercc92_ynbA.sh` | danRer7 (zebrafish) | ERCC92 & ynbA | | `src/ebwt/danRer7_ercc92_ynbA/ref`

`bin/index_hg19_ercc92_ynbA_u13369.sh` and `bin/index_hg38_ercc92_ynbA_u13369.sh` accept two options as follows.

- **1st:** Path of index location
- **2nd:** Path of fasta format file, which contains extra sequences, for example expression vectors.

#### Index of reference transcriptome

This is index for tophat. Currently there are three scripts to prepare index of the reference transcriptome. It requires two options, (i) genome version, and (ii) location of genome+spike-in index.

Script | Transcriptome | Created index location
-------|---------------|-----------------------
`bin/index_ensGene.sh` | ENSEMBL | `src/ebwt/ver_ensGene/ref`
`bin/index_knownGene.sh` | UCSC known genes | `src/ebwt/ver_knownGene/ref`
`bin/index_refGene.sh` | NCBI RefSeq | `src/ebwt/ver_refGene/ref`

> Although update of the genome sequence is rarely, the reference transcriptome updates very frequently, twice a year as usual. Therefore, you should create new transcriptome index for each project, while you can copy genome index used in the previous project.

#### Index of PhyX

You need to build the third index when you need to exclude PhyX control sequences.

- `bin/index_phyX.sh` creates index at `src/ebwt/phyX/ref`

### Design of experiments

There are two important configuration files in your STRTprep project folder, (i) `src/conf.yaml`, and (ii) `src/samples.csv`. And sometimes you may need to edit/change a barcode layout file.

#### `src/conf.yaml`; STRT library information and study design

This is information about your libraries and your study, written by [YAML](http://www.yaml.org/start.html) format, so you can edit it by any text editor. You can use [`src/template-conf.yaml`](https://github.com/shka/STRTprep/blob/v3dev/src/template-conf.yaml) for a template. It consists three sections as below.

```yaml
---
PREPROCESS:
  # parameters for preprocessing
LIBRARIES:
  # description of your libraries and the raw sequences
PLUGINS:
  # parameters for plugins
```

##### `PREPROCESS` section

This section contains preprocessing parameters. All key-value pairs below are essential.

Key | Type | Value
----|------|------
`UMI` | Integer | Length of UMI
`BARCODE` | Integer | Length of barcode
`GAP` | Integer | Length of gap
`CDNA` | Integer | Length of cDNA part
`GENOMESPIKERIBO` | String | Path of genome index
`TRANSCRIPT` | String | Path of transcriptome index
`CUSTOMTSS` | String | (Optional) Path of bed-format file, which defines transcription start regions for hypothetical genes or expression vectors; see also [Overexpression study](#overexpression-study)
`GENEMASKING` | Strings | (Optional) Genes to be masked, mainly for rescue of completely overlapping genes in gene-based analysis; see also [Rescue of overlapping gene](#rescue-of-overlapping-gene)
`LAYOUT` | String | (Optional) Path of [barcode layout file](#barcode-layout-file) as default; see also [`LIBRARIES` section](#libraries-section) and [Barcode layout file](#barcode-layout-file)
`MULTIPLEXTYPE` | String | (Optional) Barcode reads and UMI+gap+cDNA reads must be sequenced separately if you specify `C1` as `MULTIPLEXTYPE`. Otherwise (ex. no specification), the barcode+UMI+gap+cDNA must be contained in each read. See also [`LIBRARIES` section](#libraries-section) and [STRT by Fluidigm C1](#strt-by-fluidigm-c1).
`PHYX` | String | (Optional) Path of PhyX index
`QUALITYBASE` | Integer | (Optional) Either 33 (as default) or 64; see also [FASTQ format in Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format)

> You can give empty `PHYX` value when the exclusion is unnecessary.

```yaml
# Example of PREPROCESS section
PREPROCESS:
  UMI: 6
  BARCODE: 6
  GAP: 3
  CDNA: 44
  LAYOUT: src/barcodes-May2015.txt
  GENOMESPIKERIBO: src/ebwt/hg19_ercc92_ynbA_u13369/ref
  TRANSCRIPT: src/ebwt/hg19_refGene/ref
```

##### `LIBRARIES` section

This section defines your library names, and location of the raw sequence files. In the example, `TEST1` and `TEST2` at the second level are the library names themselves, and the third levels are specification of each library as below.

Key | Type | Value
----|------|------
`FASTQS` | Array (of string) | Locations of the raw sequences. It must contain barcode+UMI+gap+cDNA (as default), or at least UMI+gap+cDNA (if `MULTIPLEXTYPE == C1` at [`PREPROCESS` section](#preprocess-section)).
`FASTQS2` | Array (of string) | (Optional) Locations of the raw barcode sequences; see also [`PREPROCESS` section](#preprocess-section) and [STRT by Fluidigm C1](#strt-by-fluidigm-c1).
`LAYOUT` | String | (Optional) Path of [barcode layout file](#barcode-layout-file) for each library; see also [`PREPROCESS` section](#preprocess-section) and [Barcode layout file](#barcode-layout-file)

```yaml
# Example of LIBRARIES section
LIBRARIES:
  TEST1:
    FASTQS:
    - src/test1-1.fq.gz
    - src/test1-2.fq.gz
  TEST2:
    FASTQS:
    - src/test2-1.fq.gz
    - src/test2-2.fq.gz
```

##### `PLUGINS` section

Plugin for STRTprep is a script that adds a proximal downstream analysis. And the `PLUGINS` section if for choice of the pulgins, and to give parameters for the plugins. You can find all registered plugins and the parameters at "[Plugin framework](plugin.md)".

#### `src/samples.csv`; sample information and study design

This is sample information and the study design, written by [CSV](https://tools.ietf.org/html/rfc4180) format. You can edit it by various text or spreadsheet editor, for example Microsoft Excel. You can use [`src/template-samples.csv`](https://github.com/shka/STRTprep/blob/v3dev/src/template-samples.csv) as a template.

> Microsoft Office uses ";" (semicolon) instead of "," (colon) as column separater, when you use your computer with European environment. If so, you can split it appropriately after loading into Excel by "[Convert Text to Columns Wizard](https://support.office.com/en-gb/article/Split-names-by-using-the-Convert-Text-to-Columns-Wizard-39f7b055-6b39-4cb5-9512-13cc19b3a807)".

Important columns are as below, and you can add the other columns for description of sample properties, which would be referred by plugins.

Column | Type | Value
-------|------|------
`LIBRARY` | Word | Library name, given by `src/conf.yaml`
`WELL` | Word | Well name, based on the barcode layout, specified in `src/conf.yaml`
`NAME` | Word | Sample name; `NA` in case of empty well or ignoring
`CLASS.TFE` | Word | (Required only for the step 5) Class name for transcript assembly as TFE definition; "NA" in case of empty or ignore
`CLASS.n` (n=0, 1, ...) | Integer | (Optional) Sample class number for differential expression test `n`; 1 or 2 for two class comparison; 1, 2, 3, ... for multiclass comparison; "NA" in case of empty well or ignoring
`BLOCK.n` (n=0, 1, ...) | Integer | (Optional) Permutation block number (1, 2, ...) for the differential expression test `n`; identical number is assigned (i) to a pair for the paired comparison, (ii) to a same gender, or (iii) to a library for canceling of library bias, for example; "NA" in case of empty or ignore
`FORCE_APPROVAL` | Boolean | (Optional) If `TRUE` as for inter-library control samples etc., we can ignore the automatic outlier check for the samples. The other samples to be checked must be `FALSE`, when you use this column.

> Especially in use of Excel, there could be " " (space) characters at the end of the values unexpectedly, but these are prohibited.

> `CLASS.n`, and `BLOCK.n` for design of the differential expression tests are used by SAMstrt [[Katayama et al. 2013](http://www.ncbi.nlm.nih.gov/pubmed/?term=23995393)]. See also SAMseq [[Li and Tibshirani 2013](http://www.ncbi.nlm.nih.gov/pubmed?term=22127579)] and the manual in the [official web page](http://statweb.stanford.edu/~tibs/SAM/) to understand the statistical background and the usage.

#### Barcode layout file

This is tab-delimited text file linking between well names and barcode+gap sequence. Although you can specify your own layout and sequence, there are three predefined mapping for 48-plex STRT with 6 bp barcodes and 3 bp gap, (i) [`src/barcodes.old.txt`](https://github.com/shka/STRTprep/blob/v3dev/src/barcodes.old.txt),  (ii) [`src/barcodes.txt`](https://github.com/shka/STRTprep/blob/v3dev/src/barcodes.txt), and (iii) [`src/barcodes-May2015.txt`](https://github.com/shka/STRTprep/blob/v3dev/src/barcodes-May2015.txt).

## Special protocols

### Overexpression study

To estimate level of your over-expressed genes in the samples, you need two following files.

- **Fasta format file:** This file must contain DNA sequences around actual start site of the expression vector. To distinguish between the exogenous and endogenous expressions, the start site sequence should be unique. An example is at [`src/PiggyBacCMV.fa`](https://github.com/shka/STRTprep/blob/v3dev/src/PiggyBacCMV.fa).
- **Bed format file:** This file defines transcription start region within the vector sequence as quantification unit like as genes. Defined regions will be analyzed in gene-based analysis. An example is at [`src/PiggyBacCMV.bed`](https://github.com/shka/STRTprep/blob/v3dev/src/PiggyBacCMV.bed).

With these two files, you need to create a custom genome index.

```bash
. bin/setup_uppmax.sh
bin/index_hg19_ercc92_ynbA_u13369.sh src/ebwt/hg19_ercc92_ynbA_u13369_PiggyBacCMV src/PiggyBacCMV.fa
bin/index_refGene.sh hg19 src/ebwt/hg19_ercc92_ynbA_u13369_PiggyBacCMV/ref
```

Also, do not forget to specify `CUSTOMTSS` at `PREPROCESS` configuration.

```yaml
PREPROCESS:
  UMI: 6
  BARCODE: 6
  GAP: 3
  CDNA: 44
  LAYOUT: src/barcodes-Sep2015.txt
  PHYX:
  GENOMESPIKERIBO: src/ebwt/hg19_ercc92_ynbA_u13369_PiggyBacCMV/ref
  TRANSCRIPT: src/ebwt/hg19_refGene/ref
  CUSTOMTSS: src/PiggyBacCMV.bed
```

### Rescue of overlapping gene

Gene-based analysis (not TFE-based analysis) ignores genes completely overlapping, since it is not sure which the genes are origin of aligned reads at the overlapping regions. However, it's sometimes inconvenient, for example

- Readthrough gene definition hides important primary gene; ex. INS-IGF2, which is readthrough from INS (insulin) to IGF2 on hg19 (NCBI definition 21-Nov-2015).
- Subspecie-specific gene definition hides each other; ex. Hbb-bt, which is a haplotype of hemoglobin beta in C57BL/-type strain, overlaps Hbb-b2, which is a haplotype of hemoglobin beta in BALB/c and 129Sv on mm9 (NCBI definition 13-Nov-2015).

To rescue either overlapping genes, you can specify `GENEMASKING` at `PREPROCESS` configuration. This is an example to rescue Hbb-b2 and Prnp.

```yaml
PREPROCESS:
  UMI: 4
  BARCODE: 6
  GAP: 3
  CDNA: 46
  LAYOUT: src/barcodes.old.txt
  PHYX:
  GENOMESPIKERIBO: src/ebwt/mm9_ercc92_ynbA_bk000964/ref
  TRANSCRIPT: src/ebwt/mm9_refGene/ref
  GENEMASKING:
  - Prn
  - Hbb-bt
```

> Confirmation of readthrough level is recommended, when the rescued genes are differentially regulated.

### STRT by Fluidigm C1

You can use STRTprep not ony for single-end STRT (= barcode, UMI, gap for template-switching and cDNA at the same end of sequencing target, like in [Islam et al. 2011](http://genome.cshlp.org/content/21/7/1160.long) and [Krjutškov et al. 2016](http://humrep.oxfordjournals.org/content/31/4/844.long)) but also pair-end STRT (UMI, gap and cDNA at either end, and barcode at the other end, like by [Fluidigm C1](https://www.fluidigm.com/c1openapp/scripthub/script/2015-06/strt2fc1-protocol-1434125971861-2)). For the latter case, you must specify

- `C1` at `MULTIPLEXTYPE` of [`PREPROCESS` section](#preprocess-section)
- `FASTQS2` for each library at [`LIBRARY` section](#library-section)

This is an example.

```yaml
PREPROCESS:
  UMI: 6
  BARCODE: 8
  GAP: 3
  CDNA: 42
  GENOMESPIKERIBO: src/ebwt/hg19_ercc92_ynbA_u13369/ref
  TRANSCRIPT: src/ebwt/hg19_refGene/ref
  QUALITYBASE: 64
  MULTIPLEXTYPE: C1
LIBRARIES:
  TEST1:
    FASTQS:
      - src/Run00354_L2_1_160226_read1_indexC1-1.fq.gz
    FASTQS2:
      - src/Run00354_L2_1_160226_read2_indexC1-1.fq.gz
    LAYOUT: src/barcodes.test1.txt
```
