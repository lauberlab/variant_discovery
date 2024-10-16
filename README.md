## Summary

This reposity contains code accompanying the article "*TMEM259* alleles modulate respiratory syncytial virus infection and ER-stress-triggered apoptosis".

**Abstract**

Respiratory syncytial virus (RSV) is a main cause of infant morbidity and mortality. Susceptibility factors for severe RSV bronchiolitis in previously healthy children are unclear.
We analyzed ...

Published in: ...

## Variant calling pipeline
This folder contains a Snakefile which was used to call SNPs and Indels of our cohorts. The pipeline employs the Genome Analysis Toolkit 4 (GATK4) Best Practices Workflows for Germline short variant discovery outlined by the Broad Institute. 

Full list of tools used in this pipeline
- GATK4
- BWA
- Picard Tools
- Samtools

## Meta-analysis
This folder contains the R script that was used to conduct the meta-analysis of IRIS1+2+3 vs. LoewenKIDS. It used the R package [meta](https://github.com/guido-s/meta/).

## Colocalization analysis
This folder contains the R script that was used to conduct the co-localization analysis. It uses the R package [coloc](https://github.com/chr1swallace/coloc). 
