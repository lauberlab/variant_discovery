# Summary

This reposity contains code accompanying the manuscript "*TMEM259* alleles modulate respiratory syncytial virus infection and ER-stress-triggered apoptosis".

**Abstract**
Respiratory syncytial virus (RSV) is a main cause of infant morbidity and mortality. Susceptibility factors for severe RSV bronchiolitis in previously healthy children are unclear.
We analyzed genetic variants in 5,142 genes involved in virus sensing, interferon (IFN) signaling and effector functions in a population of n=101 previously healthy infants with severe RSV bronchiolitis. Comparing allele frequencies of the patient cohort with the exome aggregation consortium dataset (ExAC) our analysis revealed 94 non-synonymous coding single nucleotide polymorphisms (SNPs) mapping to 79 potential risk genes. Follow up of genes expressed in human lung cells showed that TMEM259 (also known as membralin), a protein involved in regulating endoplasmic reticulum (ER) stress and neuronal cell survival, restricts RSV infection. Analysis of engineered human airway cells endogenously expressing the *TMEM259* major allele or the rs77868901 minor variant showed allele-dependent modulation of RSV infection and ER-stress induced apoptosis, but not IFN activity nor induction of IFN stimulated genes. TMEM259 abundance and polymorphisms modulate RSV infection likely by modulating the ER stress response and apoptosis. This data may aid future risk stratification and development of prevention and treatment strategies for RSV.

## Variant calling pipeline
This folder contains a Snakefile which was used to call SNPs and Indels of our cohorts. The pipeline employs the Genome Analysis Toolkit 4 (GATK4) Best Practices Workflows for Germline short variant discovery outlined by the Broad Institute. 

Full list of tools used in this pipeline
- GATK4
- BWA
- Picard Tools
- Samtools

# Meta-analysis
This folder contains the R script that was used to conduct the meta-analysis of IRIS1+2+3 vs. LoewenKIDS. It used the R package [meta](https://github.com/guido-s/meta/).

# Colocalization analysis
This folder contains the R script that was used to conduct the co-localization analysis. It uses the R package [coloc](https://github.com/chr1swallace/coloc). 
