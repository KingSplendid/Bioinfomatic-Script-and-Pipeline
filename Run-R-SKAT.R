source("http://www.bioconductor.org/biocLite.R")
biocLite("SKAT")
library(SKAT)
Generate_SSD_SetID("skat.bed", "skat.bim", "skat.fam", "skat.SetID", "skat.SSD", "skat.Info", Is.FlipGenotype=TRUE)
z <- Open_SSD("skat.SSD", "skat.Info")
yb1 <- rep(1,601)
yb2 <- rep(0,548)
y.b <- c(yb1,yb2)
obj<-SKAT_Null_Model(y.b ~ 1, out_type="D")
for (i in 1:4) { 
Z <- Get_Genotypes_SSD(SSD_INFO = z , Set_Index = i, is_ID = FALSE)
#SKAT(Z, obj, method="SKAT")$p.value
p <- SKAT(Z, obj, method="SKATO")$p.value
print(p)
#SKATBinary(Z, obj, kernel = "linear.weighted", method="SKAT", method.bin="Hybrid", weights.beta=c(1,25), weights = NULL, r.corr=0, impute.method = "bestguess", is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, N.Resampling=2 *10^6, seednum=100, epsilon=10^-6, SetID=NULL)
}