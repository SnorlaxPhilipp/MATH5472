library(SPARK)

counts <- read.table("./Layer2_BC_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), 
                                               "[", 1)), y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

# SPARK
spark <- CreateSPARKObject(counts = t(counts), location = info[, 1:2], 
                           percentage = 0.1, min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 4, verbose = F, fit.maxiter = 100)
spark <- spark.test(spark, check_positive = T, verbose = F)

save(spark, file = "./Layer2_BC_spark.rds")


# for (iper in 1:10) {
#   set.seed(iper)
#   ran_index <- sample(1:nrow(info))
#   perm_info <- info[ran_index, ]
#   rownames(perm_info) <- rownames(info)
#   
#   spark <- CreateSPARKObject(counts = t(counts), location = perm_info[, 
#                                                                       1:2], percentage = 0.1, min_total_counts = 10)
#   
#   spark@lib_size <- apply(spark@counts, 2, sum)
#   spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
#                     num_core = 1, verbose = F)
#   spark <- spark.test(spark, check_positive = T, verbose = F)
#   save(spark, file = paste0("./Layer_BC_perm", iper, "_spark.rds"))
#   rm(spark)
# }




library(trendsceek)
library(spatstat)

counts <- read.table("./Layer2_BC_count_matrix-1.tsv", check.names = F)
counts <- t(counts)
counts <- genefilter_exprmat(counts, 5, 3)

vargenes_stats = calc_varstats(counts, counts, quant.cutoff = 0.9, method = "rlm")

topvar.genes = rownames(vargenes_stats[["real.stats"]])[1:500]

counts_sub = counts[topvar.genes, ]
dim(counts_sub)

cn <- colnames(counts_sub)
pp <- matrix(as.numeric(do.call(rbind, strsplit(cn, split = "x"))), ncol = 2)

pp <- ppp(x = pp[,1], y = pp[,2], range(pp[,1]), range(pp[,2]))
pp = set_marks(pp, counts_sub, log.fcn = log10)

topvar.genes = rownames(vargenes_stats[["real.stats"]])[1:500]
pp2plot = pp_select(pp, topvar.genes)

# Trendsceek
source("trendsceek.R")
set.seed(1)
trendstat_list = trendsceek_test_2(pp2plot,nrand = 50)
save(trendstat_list, time_comp, file = paste0("./BC_tsc.rds"))

# Plotting
library(ggplot2)

load("./Layer2_BC_spark.rds")
spatialDE <- read.csv("Layer2_BC_spe.csv")
alphas <- seq(0, 0.1, 0.001)
power_spark <- c()
power_spe <- c()
for (alpha in alphas){
  power_spark <- c(power_spark, sum(spark@res_mtest$adjusted_pvalue < alpha))
  power_spe <- c(power_spe, sum(spatialDE$qval < alpha))
}

powers <- data.frame(FDR = alphas, spark = power_spark, spatialDE = power_spe)

library(reshape2)
powers <- melt(powers, id.vars = "FDR",
               variable.name = "method",
               value.name = "power")
ggplot(powers, aes(x = FDR, y = power, color = method)) +
  geom_line()

# ggsave("powers_bc.png")

p_spe <- spatialDE$pval
p_spark <- spark@res_mtest$adjusted_pvalue
p_spe <- -log10(p_spe[p_spe > 0])
p_spark <- -log10(p_spark[p_spark > 0])
theoretical_spe <- -log10(seq(0,1,length.out = length(p_spe)+1))[-1]
theoretical_spark <- -log10(seq(0,1,length.out = length(p_spark)+1))[-1]
qqplot(x = theoretical_spe, y = p_spe,
       xlab = "Expected", ylab = "Observed", col = "red", pch = 20)
points(sort(theoretical_spark), sort(p_spark), col = "blue", pch = 20)
legend("topleft",c("SPARK","SpatialDE"), col = c("red","blue"),pch = 20)

