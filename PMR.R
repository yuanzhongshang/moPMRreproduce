library(PMR)
library(LaplacesDemon)

#load the scaled genenotype matrix in eQTL data (e.g. cis-SNPs of LIN9 gene from GEUVADIS data)
zx<-read.table("Zx.txt")
zx<-as.matrix(zx)
n1 = dim(zx)[1]
q = dim(zx)[2]

#load the scaled genenotype matrix in GWAS data (e.g. the same cis-SNPs from GERA data)
zy<-read.table("Zy.txt")
zy<-as.matrix(zy)
n2<-dim(zy)[1]

#set trait correlation matrix(e.g. the correlation among TC, HDL, LDL, TC in NFBC data)
squaresigma<- matrix(c(1,0.126,0.881,0.415,0.126,1,-0.138,-0.401,0.881,-0.138,1,0.330,0.415,-0.401,0.330,1),4,4)
k<-dim(squaresigma)[1] ##the number of phenotypes

#set PVE_zx to be 10%
squaresigma_beta<-0.1/q
squaresigma_x<-0.9

#set the common pleiotropy effect for each trait to be 0,0.001,0.005,0.01, respectively
gamma=matrix(c(0,0.001,0.005,0.01),k,1)

#set the causal effect to be 0
alpha=matrix(sqrt(0),k,1)

#get the simulated gene expression data
beta <- matrix(rnorm(q, 0, sd = sqrt(squaresigma_beta)),q,1)
epison_x <- matrix(rnorm(n1, 0, sd = sqrt(squaresigma_x)), n1, 1)
x <- zx %*% beta + epison_x
x<-as.matrix(x)

#get the simulated phenotype
y_mean=(gamma %*% matrix(1,1,q)) %*% t(zy)
epison_y<-rmatrixnorm(matrix(0,k,n2), squaresigma, diag(n2))
y<-y_mean+epison_y

#run the 4 PMR models using PMR_individual function
yin=as.vector(x)
zin=as.matrix(y)
x1in=zx
x2in=zy
z1in=zin[1,]
z2in=zin[2,]
z3in=zin[3,]
z4in=zin[4,]
result1<-PMR_individual(yin, z1in, x1in, x2in, method = "PMR_individual_Egger",  max_iterin = 1000, epsin = 1e-05,  Heritability_geneexpression_threshold = 1e-04)
result2<-PMR_individual(yin, z2in, x1in, x2in, method = "PMR_individual_Egger",  max_iterin = 1000, epsin = 1e-05,  Heritability_geneexpression_threshold = 1e-04)
result3<-PMR_individual(yin, z3in, x1in, x2in, method = "PMR_individual_Egger",  max_iterin = 1000, epsin = 1e-05,  Heritability_geneexpression_threshold = 1e-04)
result4<-PMR_individual(yin, z4in, x1in, x2in, method = "PMR_individual_Egger",  max_iterin = 1000, epsin = 1e-05,  Heritability_geneexpression_threshold = 1e-04)

#get the estimate of the causal effect for each trait
alpha1<-result1$causal_effect
alpha2<-result2$causal_effect
alpha3<-result3$causal_effect
alpha4<-result4$causal_effect

#get the estimate of the pleiotropy effect for each trait
gamma1<-result1$pleiotropy_effect
gamma2<-result2$pleiotropy_effect
gamma3<-result3$pleiotropy_effect
gamma4<-result4$pleiotropy_effect

#get the pvalue for the causal test for each trait
pvalue_alpha1 = result1$causal_pvalue
pvalue_alpha2 = result2$causal_pvalue
pvalue_alpha3 = result3$causal_pvalue
pvalue_alpha4 = result4$causal_pvalue

#get the pvalue for the pleiotropy effect for each trait
pvalue_gamma1 = result1$pleiotropy_pvalue
pvalue_gamma2 = result2$pleiotropy_pvalue
pvalue_gamma3 = result3$pleiotropy_pvalue
pvalue_gamma4 = result4$pleiotropy_pvalue

#################################
###################################