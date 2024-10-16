
library(coloc)

# Load your datasets
plink_gwas_data <- read.csv("filtered_gwas_plink.csv", header = TRUE)
eqtl_data <- read.csv("filtered_eqtl.csv", header = TRUE)


plink_gwas_data$P <- as.numeric(plink_gwas_data$P)
plink_gwas_data$OR <- as.numeric(plink_gwas_data$OR)
plink_gwas_data$F_A <- as.numeric(plink_gwas_data$F_A)
plink_gwas_data$F_U <- as.numeric(plink_gwas_data$F_U)


plink_gwas_data$BETA <- log(plink_gwas_data$OR) 
#
plink_gwas_data$VARBETA <- 1 / (plink_gwas_data$P * (1 - plink_gwas_data$P) * 517)


plink_gwas_data$MAF <- (plink_gwas_data$F_A + plink_gwas_data$F_U) / 2


gwas_plink <- list(
  snp = plink_gwas_data$SNP,  # snp IDs
  pvalues = plink_gwas_data$P,  # p-value
  beta = plink_gwas_data$BETA,  # Beta (log(OR))
  varbeta = plink_gwas_data$VARBETA,  # var
  N = 517,  
  MAF = plink_gwas_data$MAF,  # 
  type = "cc"  # = case-control (binary trait)
)

# columns to numeric
eqtl_data$pvalue <- as.numeric(gsub(",", ".", eqtl_data$pvalue))
eqtl_data$beta <- as.numeric(gsub(",", ".", eqtl_data$beta))
eqtl_data$se <- as.numeric(gsub(",", ".", eqtl_data$se))
eqtl_data$maf <- as.numeric(gsub(",", ".", eqtl_data$maf))  
eqtl_data$ma_samples <- as.numeric(eqtl_data$ma_samples)  


eqtl <- list(
  snp = eqtl_data$rsid,  
  pvalues = eqtl_data$pvalue,  
  beta = eqtl_data$beta,  
  varbeta = (eqtl_data$se)^2,  # variance
  N = eqtl_data$ma_samples,  
  MAF = eqtl_data$maf,  
  type = "quant"  
)


results_plink <- coloc.abf(gwas_plink, eqtl)

# Print the results for PLINK
summary_results <- results_plink$summary
write.table(as.data.frame(summary_results), file = "coloc_summary_results_plink.txt", sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(results_plink$results)) {
  snp_results <- results_plink$results
  write.table(snp_results, file = "coloc_snp_results_plink.txt", sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("No SNP-level results found.")
}

print("Results:")
print(results_plink)
