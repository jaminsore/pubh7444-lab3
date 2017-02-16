# code to conduct posterior inference & prediction
# for the linear regression model for
# the land data using conjugate NIG prior

set.seed(810973206)

# first load "MASS" package which include the function
# "mvrnorm" to sample from multivariate normal dist.
library(MASS)

# read data from file
# dir <- "C:/temp/"
land.data <- read.table(file="./land_data.txt",header=T,sep="")
ls(land.data)

# define function to generate NITER samples of (beta,sigma^2) 
# from the joint posterior using NIG prior with parameters "prior.para"
# and a given dataset "data"

post.sampling2 <- function(data, prior.para, NITER) {
    
    Y <- data[,'Y']
    X <- as.matrix(data[,2:4])
    X <- cbind(rep(1,times=length(Y)),X)
    n <- length(Y)
    p <- dim(X)[2]
    
    tXX <- t(X) %*% X
    tXX.inv <- solve(tXX)

    # extract the parameter in the NIG prior
    mu <- prior.para$mu
    V <- prior.para$V
    a <- prior.para$a
    b <- prior.para$b
    
    # calculate the posterior parameters
    V.star <- solve(solve(V) + tXX)
    mu.star <- V.star %*% (V %*% mu + t(X) %*% Y)
    a.star <- a + n/2
    b.star <- b + ( t(mu) %*% solve(V) %*% mu + t(Y) %*% Y
              - t(mu.star) %*% solve(V.star) %*% mu.star )/2
    
    
    # perform posterior sampling and return results
    sigma2 <- rep(NA, times = NITER)
    beta <- matrix(NA, nrow = NITER, ncol = p)
    colnames(beta) <- c('beta1','beta2','beta3','beta4')
    
    for (i in 1:NITER) {
        sigma2[i] <- 1/rgamma(1, a.star, rate=b.star)
        beta[i,] <- mvrnorm(1, mu.star, V.star)
    }
    
    cbind(beta,sigma2)
}

# specify the prior parameters and collect posterior samples
X <- cbind(rep(1,times=nrow(land.data)),as.matrix(land.data[,2:4]))
prior.para = list(mu = rep(0,4),
                  V = 10000 * solve(t(X) %*% X),
                  a = 0.001,
                  b = 0.001
                  )

land.samples2 <- post.sampling2(land.data,prior.para,NITER=5000)


# define the function to compile summary statistics 
sumstats <- function(vector){
    stats <- cbind(mean(vector),
                   sd(vector),
                   t(quantile(vector,c(.025,.5,.975))))
    names(stats) <- c('mean','sd','2.5%','50%','97.5%')
    stats
}

# summaries of the NEW samples given the NIG prior
t(apply(land.samples2,2,sumstats))


## Now we are to obtain the 95% posterior predictive interval

X.tilde <- c(1,0.870213,0.9453101,0.9198721)

post.predict <- function(X.tilde,post.samples) {
    NITER <- nrow(post.samples)
    y.tilde <- rep(NA,NITER)
    for (i in 1:NITER) {
        beta <- post.samples[i,1:4]
        sigma2 <- post.samples[i,5]
        y.mean <- X.tilde %*% beta
        y.tilde[i] <- rnorm(1,y.mean,sigma2)
    }
    y.tilde
}

# collect NEW predictive samples of y and calculate interval
pred.samples2 <- post.predict(X.tilde,land.samples2)

exp(quantile(pred.samples2,c(.025,.975)))