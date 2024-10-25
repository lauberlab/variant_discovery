#install.packages("meta")
library(meta)

data <- data.frame(
  comparison = c("IRIS1 vs LoewenKIDS", "IRIS2 vs LoewenKIDS", "IRIS3 vs LoewenKIDS"),
  OR = c(     2.0364, 3.7437, 0.8365),  # odds ratios of Fisher-exact test
  p_value = c(0.2711, 0.0621, 0.7362),  #  p-values of Fisher-exact test
  n_cases = c(   101,     95,     32),  # Sample sizes of cases = IRIS
  n_controls = c(367,    367,    367)   # Sample sizes of LoewenKIDS
)

# Calculates Z-scores and standard errors from p-values
data$Z  <- qnorm(1 - data$p_value / 2)
data$SE <- abs( log(data$OR) / data$Z )

# run meta-analysis
meta_result <- metagen(
  TE = log(data$OR),  
  seTE = data$SE,   
  studlab = data$comparison,
  sm = "OR",          
  method.tau = "REML",
  data = data
)

# results
summary(meta_result)

png("forest_plot_single.png", width = 800, height = 600)
forest(meta_result)
dev.off()

# Save the funnel plot as a PNG
png("funnel_plot_single.png", width = 800, height = 600)
funnel(meta_result)
dev.off()
