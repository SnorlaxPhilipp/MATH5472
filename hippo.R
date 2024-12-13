library(SPARK)
counts <- t(read.csv("./seqFISH_field43_countdata.csv", row.names = 1, 
                     check.names = F))
info <- read.csv("./seqFISH_field43_info.csv", row.names = 1)
spark <- CreateSPARKObject(counts = counts, location = info[, 1:2], 
                           percentage = 0.1, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 10, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)
save(spark, file = "./hippo_spark.rds")

winsorize <- function(x, n.win = 1) {
  n.vals = length(x)
  win.sorted.ind = order(x)
  x[win.sorted.ind[(n.vals - n.win + 1):n.vals]] = x[win.sorted.ind[n.vals - 
                                                                      n.win]]
  return(x)
}

library(trendsceek)
library(spatstat)
counts <- t(read.csv("./seqFISH_field43_countdata.csv", row.names = 1))
info <- read.csv("./seqFISH_field43_info.csv", row.names = 1)
win.counts <- t(apply(counts, 1, winsorize, 4))

counts_filt = genefilter_exprmat(win.counts, 5, 3)
xy <- info[, 1:2]
pp <- ppp(x = xy[,1], y = xy[,2], range(xy[,1]), range(xy[,2]))
pp = set_marks(pp, counts_filt, log.fcn = log10)

source("trendsceek.R")
set.seed(1)
trendstat_list = trendsceek_test_2(pp, nrand = 500, ncores = 1)
save(trendstat_list, file = "hippo_tsc.rds")

# Plotting
load("./hippo_tsc.rds")
load("./hippo_spark.rds")
spe <- read.csv("hippo_spe.csv")

alphas <- seq(0, 0.1, 0.001)
power_spark <- c()
power_spe <- c()
power_tsc <- c()
for (alpha in alphas){
  power_spark <- c(power_spark, sum(spark@res_mtest$adjusted_pvalue < alpha))
  power_spe <- c(power_spe, sum(spe$qval < alpha))
  power_tsc <- c(power_tsc, 
                 max(sum(trendstat_list$supstats_wide$p.bh_Emark < alpha),
                     sum(trendstat_list$supstats_wide$p.bh_Vmark < alpha),
                     sum(trendstat_list$supstats_wide$p.bh_markcorr < alpha),
                     sum(trendstat_list$supstats_wide$p.bh_markvario < alpha)))
}

powers <- data.frame(FDR = alphas, spark = power_spark, spatialDE = power_spe, trendsceek = power_tsc)

library(reshape2)
powers <- melt(powers, id.vars = "FDR",
               variable.name = "method",
               value.name = "power")

library(ggplot2)
ggplot(powers, aes(x = FDR, y = power, color = method)) +
  geom_line()

# ggsave("powers_hippo.png")

p_spe <- spe$pval
p_spark <- spark@res_mtest$adjusted_pvalue
p_Emark <- trendstat_list$supstats_wide$p.bh_Emark
p_Vmark <- trendstat_list$supstats_wide$p.bh_Vmark
p_markcorr <- trendstat_list$supstats_wide$p.bh_markcorr
p_markvario <- trendstat_list$supstats_wide$p.bh_markvario
p_spe <- -log10(p_spe[p_spe > 0])
p_spark <- -log10(p_spark[p_spark > 0])
p_Emark <- -log(p_Emark[p_Emark > 0])
p_Vmark <- -log10(p_Vmark[p_Vmark > 0])
p_markcorr <- -log10(p_markcorr[p_markcorr > 0])
p_markvario <- -log10(p_markvario[p_markvario > 0])

theoretical <- -log10(seq(0,1,length.out = length(p_spe)+1))[-1]
qqplot(x = theoretical, y = p_spe,
       xlab = "Expected", ylab = "Observed", col = "red", pch = 20)
points(sort(theoretical), sort(p_spark), col = "blue", pch = 20)
points(sort(theoretical), sort(p_Emark), col = "green", pch = 20)
points(sort(theoretical), sort(p_Vmark), col = "pink", pch = 20)
points(sort(theoretical), sort(p_markcorr), col = "yellow", pch = 20)
points(sort(theoretical), sort(p_markvario), col = "purple", pch = 20)

legend("topleft", c("SpatialDE","SPARK","Trendsceek.E",
                    expression(paste("Trendsceek.",rho)),
                    expression(paste("Trendsceek.",gamma)),"Trendsceek.V"), col = c("red","blue","green","yellow","purple","pink"),pch = 20)
