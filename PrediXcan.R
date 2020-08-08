library(glmnet)
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

Y1=y[1,]
Y2=y[2,]
Y3=y[3,]
Y4=y[4,]

#Run PrediXcan
fit.elasnet.cv <- cv.glmnet(zx, x, type.measure="mse", alpha=0.5, family="gaussian")
elasnet_M <- predict(fit.elasnet.cv, s=fit.elasnet.cv$lambda.min, newx=zy)
hat=coefficients(summary(lm(Y1~elasnet_M)))
pvalue1=hat[2,4]

hat=coefficients(summary(lm(Y2~elasnet_M)))
pvalue2=hat[2,4]

hat=coefficients(summary(lm(Y3~elasnet_M)))
pvalue3=hat[2,4]

hat=coefficients(summary(lm(Y4~elasnet_M)))
pvalue4=hat[2,4]
#################################
###################################