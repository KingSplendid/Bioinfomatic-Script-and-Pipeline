library(PheWAS)

genotypes <- read.table("PheWAS_Genotype.txt",header=T, sep="\t")

phenotypes <- read.csv("PCV_Phenotype.csv")

results=phewas(predictors = genotypes,outcomes = csv.phenotypes[,c("id","ph","pa","pna","pp","nodule","nodule_type11","nodule_type10","mb","drusen","sfct","mcvd","dls","irf","srf","rpei","peds","pedf","pedsa","peda","pedn","pedp","pedh","mped","bs","pachychoroid","polyps","cvh","gld","pcnv","tpcv")],min.records=1,significance.threshold=c("bonferroni"))

write.table(results,file = "result.txt", row.names = F, quote = F, sep="\t")

