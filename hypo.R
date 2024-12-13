library(SPARK)

counts <- read.csv("./MERFISH_Animal18_Bregma0.11_countdata.csv", 
                   row.names = 1, check.names = F)
info <- read.csv("./MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)

set.seed(1)
ran_index <- sample(1:nrow(info))
perm_info <- info[ran_index, ]
rownames(perm_info) <- rownames(info)

spark <- CreateSPARKObject(counts = t(counts), location = perm_info[, 
                                                                    1:2], percentage = 0, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 1, verbose = T)
spark <- spark.test(spark, check_positive = T, verbose = T)
save(spark, file = "hypo_spark.rds")

library(spatstat)
library(trendsceek)
counts <- t(read.csv("./MERFISH_Animal18_Bregma0.11_countdata.csv", 
                     row.names = 1))
info <- read.csv("./MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)
counts = genefilter_exprmat(counts, 5, 3)


xy <- info[, 1:2]
pp <- ppp(x = xy[,1], y = xy[,2], range(xy[,1]), range(xy[,2]))
pp = set_marks(pp, counts, log.fcn = log10)

source("trendsceek.R")
set.seed(1)
trendstat_list = trendsceek_test_2(pp,nrand = 50)
save(trendstat_list, time_comp, file = "./hypo_tsc.rds")


# Plotting
library(ggplot2)
load("hypo_spark.rds")
spatialDE <- read.csv("hypo_spe.csv")

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