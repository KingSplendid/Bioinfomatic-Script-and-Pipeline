install.packages("qqman")
library(qqman)

gwasResults <- read.table("C:\\Users\\Administrator\\Desktop\\hwe_gwas.txt",header=TRUE)
gwasResults

manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 100), cex = 0.5, 
          cex.axis = 0.5, col = c("blue", "orange", "grey"), suggestiveline = F, genomewideline = 9, 
          chrlabs = c(1:22,"X","Y"))
