# Martin et al. Disease consequences of higher adiposity uncoupled from its adverse metabolic effects using mendelian randomisation. eLife, 11:e72452, 2022.

# R scripts used to conduct analysis and generate results tables

Scripts are split into three parts:
* 1-MR_analysis.R - conducts Mendelian Randomisation for each disease and sub-disease using IVW, weighted median, MR-Egger and penalised weighted median methods. Generates Supplementary File 1g(i-iv).
* 2-meta_analysis.R - conducts random-effects meta-analysis of published GWAS and FinnGen IVW MR results together for appropriate matching outcomes. Generates Supplementary File 1e, 1f and 1h.
* 3-correction_FDR_BH.R - calculates Benjamini-Hochberg false discovery rates of the 37 disease IVW MR results.
