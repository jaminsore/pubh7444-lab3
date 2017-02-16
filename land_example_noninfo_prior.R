# code to conduct posterior inference & prediction
# for the linear regression model for
# the land data using a noninformative prior

set.seed(810973206)

# first load "MASS" package which include the function
# "mvrnorm" to sample from multivariate normal dist.
library(MASS)

# read data from file
# dir <- "C:/temp/"
land.data <- read.table("./land_data.txt",header=T,sep="")
ls(land.data)

# define function to generate NITER samples of (beta,sigma^2) 
# from the joint posterior using noninformative prior
# and a given dataset "data"

post.sampling1 <- function(data, NITER) {
    
        # extract the responses Y 
        Y <- data[,'Y']
        
        # extract the three predictors X1,X2,X3 
        X <- as.matrix(data[,2:4])
        
        # combine with the column of 1 (intercept term) 
        # to form the design matrix
        X <- cbind(rep(1,times=length(Y)),X)
        
        n <- length(Y)
        p <- dim(X)[2]
        
        # calculate (X^TX)^-1, frequentist unbiased estimators
        tXX.inv <- solve(t(X) %*% X)
        beta.hat <- tXX.inv %*% t(X) %*% Y
        s2 <- t(Y - X%*%beta.hat)%*%(Y - X%*%beta.hat)/(n - p)
        
        # define vector and matrix to store posterior samples 
        sigma2 <- rep(NA, times = NITER)
        beta <- matrix(NA, nrow = NITER, ncol = p)
        colnames(beta) <- c('beta1','beta2','beta3','beta4')
    
        # collect NITER samples from the joint posterior
        # using the for loop function
        for (i in 1:NITER) {
            sigma2[i] <- 1/rgamma(1, (n-p)/2, rate=(n-p)*s2/2)
            beta.var <- (tXX.inv) * sigma2[i]
            beta[i,] <- mvrnorm(1, beta.hat, beta.var)
        }

        # combine the posterior samples of beta and sigma^2 
        # by column, and return the result
        cbind(beta,sigma2)
}

# use the defined function to collect NITER=5000 samples 
# of beta and sigma^2 for the land dataset
land.samples <- post.sampling1(land.data,NITER=5000)


# define the function to compile summary statistics
# including mean, sd, 2.5th, 50th, and 97.5th quantile
# for a vector of samples

sumstats <- function(vector){
    stats <- cbind(mean(vector),
                   sd(vector),
                   t(quantile(vector,c(.025,.5,.975))))
    names(stats) <- c('mean','sd','2.5%','50%','97.5%')
    stats
}

# use the "apply" function to obtain summary statistics 
# for each column of matrix of posterior samples
# check the "apply" function with the command "?apply"
t(apply(land.samples,2,sumstats))

## Now we are to obtain the 95% posterior predictive interval

# new observation of X.tilde
X.tilde <- c(1,0.870213,0.9453101,0.9198721)

# define function to collect sample from 
# posterior predictive distribution
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

# use the function to collect predictive samples of y
pred.samples <- post.predict(X.tilde,land.samples)

# get the 2.5th and 97.5th quantile from the samples of y_tilde
# a 95% predictive interval of selling price is 
# exponential of the quantiles
exp(quantile(pred.samples,c(.025,.975)))

