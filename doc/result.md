# Interpretation of the results

## Table of contents

1. [Quality check](#quality-check)
2. [Differential expression tests](#differential-expression-tests)
3. High resolution analysis
4. [For Genome Browsers](#for-genome-browsers) to understand the results

## Quality check

### Typical protocol

The most important file in the quality check reports is `out/byGene/samples.xls`, which is a list of qualified samples and the statistics on quality. You must check whether there are enough number of qualified samples for your studies. If not, you should check the statistics also about the rejected samples at `out/byGene/samples_all.csv`; it may suggest that you can rescue some rejected samples by the `FORCE_APPROVAL` option in `src/samples.csv`.

### Detail of quality check report files

#### `out/byGene/samples_all.csv`

This table is a quality check report by samples. The following is the table legend. Several columns were identical with `src/samples.csv`; the followings are extra columns as the quality check reports.

Columns | Value
--------|------
`QUALIFIED_READS` | Read count just before the redundant-read exclusion
`TOTAL_READS` | Read count just after the redundant-read exclusion
`REDUNDANCY` | = `QUALIFIED_READ` / `TOTAL_READS`
`MAPPED_READS` | Read count, which are aligned at only one location
`MAPPED_RATE` | = `MAPPED_READS` / `TOTAL_READS`
`SPIKEIN_READS` | Read count, which are aligned on spike-in sense strand
`MAPPED/SPIKEIN` | ~ Relative abundance of endogenous polyA+ transcripts
`SPIKEIN_5END_READS` | Read count, which are aligned within 5'-end 50 nt region of spike-in sense strand
`SPIKEIN_5END_RATE` | = `SPIKEIN_5END_READS` / `SPIKEIN_READS`
`CODING_READS` | Read count, which are aligned within any exon or the proximal upstream (500bp) of coding gene sense strand
`CODING_5END_READS` | Read count, which are aligned within 5'-UTR or the proximal upstream of coding gene sense strand
`CODING_5END_RATE` | = `CODING_5END_READS` / `CODING_READS`
`SPIKEIN_READS.OUTLIER` | `TRUE` if it is outlier of `SPIKEIN_READS` in the library
`MAPPED/SPIKEIN.OUTLIER` | `TRUE` if it is outlier of `MAPPED/SPIKEIN` in the library
`SPIKEIN_5END_RATE.OUTLIER` | `TRUE` if it is outlier of `SPIKEIN_5END_RATE` in the library
`CODING_5END_RATE.OUTLIER` | `TRUE` if it is outlier of `CODING_5END_RATE` in the library

> You can find the distribution and the outliers by figures `out/byGene/fig_outliers_*.pdf`.

#### `out/byGene/samples.xls`

It contains only the qualified samples in `out/byGene/samples_all.csv` - all samples must be that (i) `NAME` is not `NA`, and (ii) all `*.OUTLIER` values are `FALSE`. Content of `out/byGene/samples.csv` is identical with this, but the format is different.

## Differential expression tests

### Typical protocol

Two figures, (i) `out/byGene/plugin_correlation_samples_global.pdf` and (ii) `out/byGene/plugin_pca_global.pdf`, are the most important in the first step of interpretation of the results. These are unsupervised clustering of the qualified samples based on expression levels of significantly fluctuated genes. Therefore, these figures would suggest major variation in your samples, library bias, sample batch bias, outlier samples in the expression profile, and so on.

For example, in case of comparison between cases and controls, the first branch of dendrogram in the sample correlation map (i), and the first principal component in the PCA map (ii), should be separate the case samples and the control samples. If not, for instance the first branch separated your samples by the libraries, and the second branch separated between the case and the controls, you should give different `BLOCK` number by libraries to avoid the library bias in the pairwise comparison.

If you find the appropriateness of your test design, you can extract more detailed results from `out/byGene/diffexp.xls`. For example,

1. You can find the differentially regulated genes from this table. Reasonable thresholding is by both of `qvalue` and `fluctuation`, but sometimes `pvalue` and `fluctuation` may work better when you prepared only a few biological replicates.
2. You can find the normalized expression levels in your genes of interest, and draw the figure by any of your favorite tools, for example R or Excel.
3. You can apply the other statistics using the other third-party tools with the raw read counts or the normalized expression levels.

Also the plugins in STRTprep would help your analysis. You may need to edit `src/conf.yaml` again and run more several times to get informative reports/figures/tables from the plugins.

### Detail of differential expression test report files

#### `out/byGene/diffexp.xls`

This is a big table containing expression levels and statistical results of the differential expression in the qualified samples. Content of `out/byGene/diffexp.csv` is identical with this, but the format is different.

Column | Value
-------|------
`Gene` | Gene name
`fluctuationScore.global` | Statistical score in degree of variation of the normalized expression levels between the qualified wells towards expected technical noise levels
`fluctuation.global` | Adjusted p-value of the degree of variation between the qualified wells
`diffexpScore.n` (n=0, 1, ...) | Statistical score of differential expression in the test `n`; `n` is  test ID given by `src/samples.csv`. Positive score in two class comparison is up-regulation of the class 2 samples than the class 1, and negative is down-regulation. In case of multiclass comparison, the score is all positive. Larger absolute score is more significant.
`pvalue.n` (n=0, 1, ...) | P-value of differential expression in the test `n`; corrected by Benjamini & Hochberg method
`qvalue.n` (n=0, 1, ...) | Q-value (a.k.a. FDR) of differential expression in the test `n`; corrected by Storey and Tibshirani method
`fluctuationScore.n` (n=0, 1, ...) | Statistical score in degree of variation of the normalized expression levels between targets of the test `n` towards expected technical noise levels
`fluctuation.n` (n=0, 1, ...) | Adjusted p-value of the degree of variation between the targets of the test `n`
`N`&#124;*`library.well`*&#124;*`name`* | Spike-in based normalized expression level
`R`&#124;*`library.well`*&#124;*`name`* | Raw aligned read count within 5'-UTR of the gene, or the proximal (500 nt) upstream; such transcripts starting these regions would be intact mRNAs as template for translation into proteins. You can find the regions [by Genome Browsers](Results#for-genome-browsers) with file `out/web/regions_byGene.bed.gz`.

> `diffexpScore.n`, `pvalue.n` and `qvalue.n` are calculated by SAMstrt called from STRTprep [[Katayama et al. 2013](http://www.ncbi.nlm.nih.gov/pubmed/?term=23995393)]. Briefly, in case of two class comparison, the score is average of Wilcoxon statistics with multiple Poisson resampling. In contrast, the score in multiple class comparison is average of Kruskal-Wallis statistic. See also SAMseq [[Li and Tibshirani 2013](http://www.ncbi.nlm.nih.gov/pubmed?term=22127579)] and the manual in the [official web page](http://statweb.stanford.edu/~tibs/SAM/) to understand the statistical background and the outputs.

> `fluctuationScore.*` and `fluctuation.*` are calculated by STRTprep itself [Krjutškov and Katayama et al., submitted]. Currently it was calculated by libraries for `global` test samples, or by blocks for the other tests, then merged.

> [AutoFilter](https://support.office.com/en-ca/article/Quick-start-Filter-data-by-using-an-AutoFilter-08647e19-11d1-42f6-b376-27b932e186e0?ui=en-US&rs=en-CA&ad=CA) by Microsoft Excel is very useful to select the significantly fluctuated genes. After the filtering, copy the gene names, then paste them into your favorite analysis tools!

## For Genome Browsers

Files in a folder `out/web` are appendix for understanding of the results on Genome Browsers; for example, you can [upload to UCSC Genome Browser](http://genome-euro.ucsc.edu/FAQ/FAQcustom.html#custom1).

### `out/web/regions_byGene.bed.gz`

The bands in this track are 5'-UTR of protein coding genes, or the proximal (500 nt) upstream. Aligned reads within the bands were used for the `byGene` analysis.

### `out/web/regions_byTFE.bed.gz`

The bands in this track are putative first exon according to assembling of STRT reads. Thick position in each band represents most frequent site of 5'-end of STRT reads, which would be major transcription start site. Aligned reads within the bands were used for the `byTFE` analysis.

### `out/web/*_transcripts.gff.gz`

These tracks are the putative transcripts by assembling of STRT reads.
