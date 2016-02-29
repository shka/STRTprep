# Plugin framework and extension of STRTprep

Plugin for STRTprep is a script that adds a proximal downstream analysis. Output files by the plugins have common prefix name `out/byGene/plugin_` or `out/byTFE/plugin_` followed by the plugin name and the test ID.

You can specify parameters for the plugins from `src/conf.yaml`. The second level is the plugin name, the third level is test ID, which you specify by `src/samples.csv`, and the fourth level is actual plugin parameters. In the following example, `transcripts` plugins was applied to tests `global` and `1`, but not to a test `0`.

```yaml
# Example of "transcripts" plugin parameters.
PLUGINS:
  transcripts:
    global:
      FLUCTUATIONP: 0.05
      COLOR: LIBRARY
      CLASSES:
      - LIBRARY
    1:
      FLUCTUATIONP: 0.05
      DIFFEXPQ: 0.05
      COLOR: GENDER
```

A list below is available plugins and the help documents.

- [`correlation_samples`](#plugin-correlation_samples)
- [`heatmap_diffexp`](#plugin-heatmap_diffexp)
- [`pca`](#plugin-pca)
- [`simple_diffexp`](#plugin-simple_diffexp)
- [`simple_gsea`](#plugin-simple_gsea)
- [`transcripts`](#plugin-transcripts)

Moreover,  you can add your plugins into STRTprep; see "[Extension of STRTprep](#extension-of-strtprep)".

## Anchor (&) and alias (\*)

```yaml
# Example without anchor & alias
PLUGINS:
  correation_samples:
    0:
      ANNOTATIONS:
      - DIAGNOSIS
      - TREATMENT
      FLUCTUATIONP: 0.05
    1:
      ANNOTATIONS:
      - DIAGNOSIS
      - TREATMENT
      FLUCTUATIONP: 0.05
    2:
      ANNOTATIONS:
      - DIAGNOSIS
      - TREATMENT
      FLUCTUATIONP: 0.05
  heatmap_diffexp:
    0:
      LABELS:
      - Control
      - Case
      ANNOTATIONS:
      - DIAGNOSIS
      - TREATMENT
      - GENDER
      FLUCTUATIONP: 0.05
      DIFFEXPQ: 0.05
    1:
      LABELS:
      - Control
      - Acute
      - Chronic
      ANNOTATIONS:
      - DIAGNOSIS
      - TREATMENT
      - GENDER
      FLUCTUATIONP: 0.05
      DIFFEXPQ: 0.05
```

```yaml
# Example with anchor & alias
PLUGINS:
  default0: &default_thresholds
    FLUCTUATIONP: 0.05
    DIFFEXPQ: 0.05
  default1: &default_annotationsA
    ANNOTATION:
    - DIAGNOSIS
    - TREATMENT
  default2: &default_annotationsB
    ANNOTATION:
    - DIAGNOSIS
    - TREATMENT
    - GENDER
  correation_samples:
    0:
      <<: [*default_thresholds, *default_annotationsA] # DIFFEXPQ is ignored
    1:
      <<: [*default_thresholds, *default_annotationsA] # DIFFEXPQ is ignored
    2:
      <<: [*default_thresholds, *default_annotationsA] # DIFFEXPQ is ignored
  heatmap_diffexp:
    0:
      LABELS:
      - Control
      - Case
      <<: [*default_thresholds, *default_annotationsB]
    1:
      LABELS:
      - Control
      - Acute
      - Chronic
      <<: [*default_thresholds, *default_annotationsB]
```

## Plugin `correlation_samples`

This plugin draws a heatmap of correlation coefficients between the samples, and a dendrogram of the correlation matrix, to elucidate variations between the samples. The heatmap and dendrogram are labeled by the sample names, and annotated by specified sample properties.

Parameter key | Type | Value
--------------|------|------
`FLUCTUATIONP` | Real, 0~1 | Threshold of fluctuation p-value
`ANNOTATIONS` | Words | Column name(s) of `src/samples.csv` to be annotated

```yaml
# Example of "correlations_samples" plugin parameters
PLUGINS:
  global:
    ANNOTATIONS:
    - LIBRARY
    - DIAGNOSIS
    - TREATMENT
    FLUCTUATIONP: 0.05
  0:
    ANNOTATIONS:
    - DIAGNOSIS
    - TREATMENT
    FLUCTUATIONP: 0.05
```

## Plugin `heatmap_diffexp`

This plugin draws (i) a expression heatmap of differentially regulated and/or fluctuated genes in the qualified samples, (ii) a dendrogram between the samples, and (iii) a dendrogram between the genes. The heatmap and dendrograms are labeled by names of the samples and the genes. Moreover, the samples are annotated by the specified properties.

Parameter key | Type | Value
--------------|------|------
`FLUCTUATIONP` | Real, 0~1 | (Optional) Threshold of fluctuation p-value
`DIFFEXPQ` | Real, 0~1 | (Optional; ignored in test `global`) Threshold of differential expression q-value
`DIFFEXPP` | Real, 0~1 | (Optional; ignored in test `global`) Threshold of differential expression p-value
`ANNOTATIONS` | Words | Column name(s) of `src/samples.csv` to be annotated
`LABELS` | Words | (Optional; ignored in test `global`) Class labels

```yaml
# Example of "heatmap_diffexp" plugin parameters
PLUGINS:
  heatmap_diffexp:
    global:
      ANNOTATIONS:
      - LIBRARY
      - DIAGNOSIS
      - TREATMENT
      FLUCTUATIONP: 0.05
    0:
      LABELS:
      - Control
      - Case
      ANNOTATIONS:
      - DIAGNOSIS
      - TREATMENT
      FLUCTUATIONP: 0.05
      DIFFEXPQ: 0.05
```

## Plugin `pca`

This plugin draws a PCA plot in the qualified samples with the fluctuated gene expression.

Parameter key | Type | Value
--------------|------|------
`FLUCTUATIONP` | Real, 0~1 | Threshold of fluctuation p-value
`COLOR` | Word | Column name of `src/samples.csv` to be distinguished by point colors
`POINT` | Word | (Optional) Column name of `src/samples.csv` to be distinguished by point characters
`COMPONENTS` | Integer | Number of principal components to be drawn

```yaml
# Example of "pca" plugin parameters
PLUGINS:
  pca:
    global:
      FLUCTUATIONP: 0.05
      COLOR: LIBRARY
      POINT: TREATMENT
      COMPONENTS: 5
    0:
      FLUCTUATIONP: 0.05
      COLOR: LIBRARY
      POINT: TREATMENT
      COMPONENTS: 3
```

## Plugin `simple_diffexp`

This plugin creates a subset of diffexp.csv, which contains only target samples and the significant features (ex. genes). It will provide all features if no specification of the significance thresholds.

Parameter key | Type | Value
--------------|------|------
`FLUCTUATIONP` | Real, 0~1 | (Optional) Threshold of fluctuation p-value
`DIFFEXPQ` | Real, 0~1 | (Optional; ignored in test `global`) Threshold of differential expression q-value
`DIFFEXPP` | Real, 0~1 | (Optional; ignored in test `global`) Threshold of differential expression p-value

```yaml
# Example of "simple_diffexp" plugin parameters
PLUGINS:
  simple_diffexp:
    global:
      FLUCTUATIONP: 0.05
    0:
      FLUCTUATIONP: 0.05
      DIFFEXPQ: 0.05
```

## Plugin `simple_gsea`

This plugin tests enrichment of  [MSigDB](http://software.broadinstitute.org/gsea/msigdb/index.jsp) genes in the differentially expressed genes; there are 13,311 gene sets in total at the version 5.1 in Jan 27, 2016. Significance of the enrichment is, simply, based on Pearson's Chi-squared test with Yates' continuity correction. This plugin is only for `byGene` quantitation, and not for global comparison. Several parameters are compatible with `heatmap_diffexp`; it is expected that you will apply this plugin after `heatmap_diffexp`.

* Output file `out/byGene/plugin_simple_gsea_n.csv` is the test result.
  * Column `ENRICHMENTP` is the significance with Benjamini & Hochberg correction.
  * `RESULT` == `OR` means over-representation of member genes in a gene set within the differentially expressed genes. In detail, (i) `ENRICHMENTP` is less than the threshold (see the parameter below), (ii) `EXP1` > `OBS1`, and (iii) `EXP1` > 5.
  * `EXP1` is observed number of significantly regulated member gene, while 'OBS1' is the expected number.
  * Sum of `EXP1`, `EXP2`, `EXP3` and `EXP4` is all genes detected by target samples in the comparison.
* Folder `plugin_dimple_gsea_n.figures` contains heatmaps, which illustrate expression profile of member genes in the over-represented gene sets (`RESULT` == `OR`)

Parameter key | Type | Value
--------------|------|------
`FLUCTUATIONP` | Real, 0~1 | (Optional) Threshold of fluctuation p-value
`DIFFEXPQ` | Real, 0~1 | (Optional) Threshold of differential expression q-value
`DIFFEXPP` | Real, 0~1 | (Optional) Threshold of differential expression p-value
`ANNOTATIONS` | Words | Column name(s) of `src/samples.csv` to be annotated in the heatmaps
`LABELS` | Words | (Optional) Class labels
`MSIGDB` | Word | Location of MSigDB xml file; you can download from the [website](http://software.broadinstitute.org/gsea/downloads.jsp)
`ENRICHMENTP` | Real, 0~1 | Threshold of the corrected significance of chi-squared test
`HOMOLOGENE$FILE` | Word | (Optional) Location of [NCBI HomoloGene](http://www.ncbi.nlm.nih.gov/homologene) table file; you can download from [FTP site](ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current). You can omit when you analyze human samples.
`HOMOLOGENE$TAXID` | Integer | (Optional) Taxonomy ID, which is maintained at [NCBI Taxonomy](http://www.ncbi.nlm.nih.gov/taxonomy/), of your sample. You can omit when you analyze human samples.

```yaml
# Example of "simple_gsea" plugin parameters
PLUGINS:
  simple_gsea:
    0:
      ANNOTATIONS:
      - DIAGNOSIS
      - TREATMENT
      FLUCTUATIONP: 0.05
      DIFFEXPQ: 0.05
      MSIGDB: src/msigdb_v5.1.xml
      HOMOLOGENE:
        FILE: src/homologene-68.data
        TAXID: 10090
      ENRICHMENTP: 0.05
```

## Plugin `transcripts`

This plugin draws boxplots in the relative poly(A)+ transcript contents.

Parameter key | Type | Value
--------------|------|------
`COLOR` | Word | Column name of `src/samples.csv` to be distinguished by point colors
`CLASSES` | Words | (Only for test `global`) Column name(s) of `src/samples.csv` to be distinguished as CLASSES
`LABELS` | Words | (Ignored in test `global`) Class labels

```yaml
# Example of "transcripts" plugin parameters
PLUGINS:
  transcripts:
    global:
      COLOR: DIAGNOSIS
      CLASSES:
      - LIBRARY
      - TREATMENT
    0:
      COLOR: DIAGNOSIS
      LABELS:
      - Control
      - Case
```

## Extension of STRTprep
(on going)
