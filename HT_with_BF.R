# R code to calculate Bayes Factor for hypothesis testing
# of H0: beta_4=0 vs H1: beta_4 <> 0 
# for the land example

# load R package that contains the function 'dmvt'
# to calculate the pdf of MVSt 
library(mvtnorm)

# read data from file
dir <- "C:/Users/bsorenson/pubh7444-lab3/"
land.data <- read.table(file=paste0(dir,"land_data.txt"),header=T,sep="")
ls(land.data)

Y = land.data$Y

# Frist, specify likelihood and prior parameters of model M_0 and M_1
X.all = cbind(rep(1,times=nrow(land.data)),as.matrix(land.data[,2:4]))

M0.para <- list(X = X.all[,-4],
                mu = rep(0,3),
                V = 10^4 * solve(t(X.all[,-4]) %*% X.all[,-4]),
                a = 0.001,
                b = 0.001
                )


M1.para <- list(X = X.all,
                mu = rep(0,4),
                V = 10^4* solve(t(X.all) %*% X.all),
                a = 0.001,
                b = 0.001
                )    

# define the function to calculate the posterior marginal 
# of a given model m(y|M)
log.post.marg <- function(Y,M.para) {
    
    n <- length(Y)
    X <- M.para$X
    mu <- M.para$mu
    V <- M.para$V
    a <- M.para$a
    b <- M.para$b
    
    y.nu <- 2 * a
    y.mean <- X %*% mu
    y.cov <- b / a * (diag(n) + X %*% V %*% t(X))
    
    dmvt(Y,y.mean,y.cov,df=y.nu)
}

# calculate the Bayes factor of M0 over M1
BF = exp(log.post.marg(Y,M0.para)-log.post.marg(Y,M1.para))
BF
