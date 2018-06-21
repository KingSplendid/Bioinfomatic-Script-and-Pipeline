install.packages("devtools")

install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))

devtools::install_github("PheWAS/PheWAS")

library(PheWAS)

genotypes <- read.table("PheWAS_Genotype.txt",header=T, sep="\t")

phenotypes <- read.csv("PCV_Phenotype.csv")

results=phewas(predictors = genotypes,outcomes = csv.phenotypes[,c("id","ph","pa")],min.records=1,significance.threshold=c("bonferroni"))

write.table(results,file = "result.txt", row.names = F, quote = F, sep="\t")

