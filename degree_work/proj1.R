setwd("/home/sergey/GAGE/Staphylococ/alig_bwt")
library("entropy")
pl1_norm <- read.csv('pl1.txt', sep="\n")
pl1 <- read.csv('pl11.txt', sep="\n")
pl2 <- read.csv('pl2.txt', sep="\n")

p1 = pl1[,1]
p1_norm = pl1_norm[,1]
p2 = pl2[,1]
#p1 <- as.matrix(sapply(pl1, as.numeric))
#p2 <- as.matrix(sapply(pl2, as.numeric))

pv1 <- as.vector(sapply(p1, as.numeric))
pv1_norm <- as.vector(sapply(p1_norm, as.numeric))
pv2 <- as.vector(as.numeric(levels(p2))[p2])

entropy::chi2.plugin(pv1[8:49], pv2[8:49])
entropy::chi2.plugin(pv1_norm[8:49], pv2[8:49])

entropy::chi2.plugin(pv2[8:49], pv1[8:49])
entropy::chi2.plugin(pv2[8:49], pv1_norm[8:49])

chisq.test(p1_norm[8:49], p2[8:49])

cor(pv1_norm[8:49], pv2[8:49], method="spearman")
cor(pv1_norm[8:49], pv2[8:49], method="pearson")

cor(pv2[8:49], pv1_norm[8:49], method="spearman")
cor(pv2[8:49], pv1_norm[8:49], method="pearson")
res <- 1 - ecdf(40)

P = ecdf(pv2)    # P is a function giving the empirical CDF of X
P(40)  

pchisq(pv2[8:49], df=40)
#chi2inv(0.95, 40)
