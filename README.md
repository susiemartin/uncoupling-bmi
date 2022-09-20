# Martin et al. Disease consequences of higher adiposity uncoupled from its adverse metabolic effects using mendelian randomisation. eLife, 11:e72452, 2022.

Scripts used to conduct original analyses and updated results based on larger GWAS releases since publication.

## R scripts used to conduct analysis and generate results tables

Scripts are split into three parts:

### 1-MR_analysis.R

Conducts Mendelian Randomisation for each disease and sub-disease using IVW, weighted median, MR-Egger and penalised weighted median methods. Generates **Supplementary File 1g(i-iv)**. Requires the following arguments and will need to be run per exposure trait (BMI, BFP, FA or UFA):
```
Rscript 1-MR_analysis.R test_trait snp_info summ_stats_file
```
where
`test_trait` is the name of the exposure (BMI, BFP, FA or UFA),

`snp_info` is the location and name of the exposure summary statistics file, and

`summ_stats_file` is the location and name of the file containing individual BetaYG and seBetaYG for the disease outcomes.

**N.B.** The `snp_info` file should contain the summary statistics of the exposure genetic variants (as given in **Supplementary File 1d**) as well as appropriate proxy variants, and the `summ_stats_file` file should contain the exposure genetic variants (or best available proxies) extracted from the disease outcome GWAS.

### 2-meta_analysis.R

Conducts random-effects meta-analysis of published GWAS and FinnGen IVW MR results together for appropriate matching outcomes. Generates **Supplementary File 1e, 1f and 1h**. Requires the output from 1-MR_analysis.R as well as a file of diseases to be meta-analysed together (including meta-analysed trait name) titled `disease_meta_analysis.txt`, of the form:
```
Trait     Meta_name
CAD_GWAS  CAD_meta
CAD_FG    CAD_meta
PCOS_GWAS PCOS_meta
...       ...
```

### 3-correction_FDR_BH.R

Calculates Benjamini-Hochberg false discovery rates of the 37 disease IVW MR results. Requires output from 2-meta_analysis.R only.

## Updated results tables and plots based on larger GWAS releases since publication

The pdfs containing figures and supplementary material relating to updates based on larger FinnGen and published GWAS releases to be found in the `updated_results_*` folders. The write-ups corresponding to these releases can be found in the Comments of the publication (accessible here: https://elifesciences.org/articles/72452). The FinnGen releases used for each update are given below:

`updated_results_sept2022` (to be uploaded soon) - FinnGen release R7.

`updated_results_april2022` - FinnGen release R6.
