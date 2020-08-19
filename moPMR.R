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

#run moPMR 
xin=as.vector(x)
Yin=as.matrix(t(y))
Zxin=zx
Zyin=zy

result<-moPMR_individual(xin, Yin, Zxin, Zyin, max_iterin = 1000, epsin = 1e-05,  Heritability_geneexpression_threshold = 1e-04)

#get the estimate of the causal effect
alpha<-result$causal_effect

#get the estimate of the pleiotropy effect
gamma<-result$pleiotropy_effect

#get the pvalue for the causal test
pvalue_alpha = result$causal_pvalue

#get the pvalue for the pleiotropy effect
pvalue_gamma = result$pleiotropy_pvalue
#############################################
###################################################
