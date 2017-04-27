library(spam)
library(magic)
library(RSEIS)
library(timeDate)

#function to create the diagonal matrix for pre-processing
Diagonal <- function(N){
    D1 <- diag(x=sample(c(1,-1),N,replace=TRUE),nrow=N) #Diagonal Matrix
    return(D1)
}

#Function to generate the G matrix(mxN)
G_mat <- function(N,m,method){
    G <- matrix(nrow=m,ncol=N) #Transformation Matrix
    if (method=="Unstructured"){
        for (i in seq(1:N)){
            G[,i] <- rnorm(m)
        }
    }
    else if (method=="Circulant"){
        G <- circulant(rnorm(N))[1:m,]
    }
    else if (method=="Toeplitz"){
        G <- as.matrix(toeplitz.spam(rnorm(N),rnorm(N))[1:m,])
    }
    else if (method=="Hankel"){
        G <- mirror.matrix(as.matrix(toeplitz.spam(rnorm(N),rnorm(N))[1:m,]))
    }
        G = G/sqrt(m)
        return(G)
}

# Main JLT function
JLT <- function(N,m,data_size,method){
    X <- matrix(nrow=N,ncol=data_size) # Original dataset
    Y <- matrix(nrow=m,ncol=data_size) #Dataset with reduced dimensions
    
    # Generate the Dataset
    for (i in seq(1:data_size)){
        X[,i] <- runif(N)
    }
    
    #Dimension Reduction depending on value of method
    G <- G_mat(N,m,method)
    if (method == "Unstructured"){
        Y = G%*%X
    } else {
        Y1 = Diagonal(N) %*% X #Preprocessing
        Y = G%*%Y1
    }
    
    #Calculating Delta
    Delta <- as.numeric(vector(length=20))
    for (i in seq(1:20)){
        W <- sample(seq(1:data_size),2,replace=FALSE)
        Euc_dist_Y = sqrt(sum((Y[,W[1]]-Y[,W[2]])^2))
        Euc_dist_X = sqrt(sum((X[,W[1]]-X[,W[2]])^2))
        Delta[i]=  abs(Euc_dist_Y - Euc_dist_X)
    }
    return(mean(Delta))
}
# Calling the JLT function : Calculate Delta for different values of m
M <- data.frame(x=1:30,y=1:30)
m_start <- round(log(100))
for (i in seq(1:30)){
    M$x[i] <- m_start+i
    M$y[i] <- JLT(N = 100,m = m_start+i,data_size =1000,method="Hankel")
}
# Using linear regression to find optimal range for m
model1 <- lm(M$y ~ (M$x))
m1 <- ceiling((0.5-coef(model1)[[1]])/coef(model1)[[2]])
plot((M$x),M$y,xlab="Reduced Dimension(m)",ylab="Delta Values",main="Fast JLT using Hankel Matrix")
legend("topright",c("N = 100","Data size = 1000",paste("m*=",m1)),cex=0.7)
abline(model1,col="gray20")
abline((0.5),0,col='red')
summary(model1)


