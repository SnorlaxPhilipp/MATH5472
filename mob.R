library(SPARK)

counts <- read.table("./Rep11_MOB_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)), 
                         y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

# SPARK
spark <- CreateSPARKObject(counts = t(counts), location = info[, 1:2], 
                           percentage = 0.1, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 1, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)

save(spark, file = "./MOB_spark.rds")


# Code for Trendsceek
# Due to storage limit the running cannot complete
library(trendsceek)
counts <- t(read.table("./Rep11_MOB_count_matrix-1.tsv", check.names = F))
tot <- apply(counts, 2, sum)
counts_filt = genefilter_exprmat(counts, 5, 3)
dim(counts_filt)

quantile.cutoff = 0.9
method = "rlm"
vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = 0.9,
                               method = "rlm")
n.topvar = 500
topvar.genes = rownames(vargenes_stats[["real.stats"]])[1:n.topvar]
counts_sub = counts_filt[topvar.genes, ]

cn <- colnames(counts_sub)
xy <- matrix(as.numeric(do.call(rbind, strsplit(cn, split = "x"))), ncol = 2)
pp <- ppp(x = xy[,1], y = xy[,2], range(xy[,1]), range(xy[,2]))
pp = set_marks(pp, counts_sub, log.fcn = log10)
topvar.genes = rownames(vargenes_stats[["real.stats"]])[1:n.topvar]
pp2plot = pp_select(pp, topvar.genes)

source("trendsceek.R")
i = 1
set.seed(i)
trendstat_list = trendsceek_test_2(pp,nrand = 50)
save(trendstat_list, time_comp, file = "./MOB_tsc.rds")

# Plotting
library(ggplot2)

load("./MOB_spark.rds")
spatialDE <- read.csv("MOB_spe.csv")
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

# ggsave("powers_mob.png")

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
