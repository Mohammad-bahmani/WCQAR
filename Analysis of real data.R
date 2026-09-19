###########################################

####################========== Packages
#install.packages("ald")
#install.packages("MASS")
#install.packages("VGAM")
#install.packages("quantreg")
#install.packages("GIGrvg")
#install.packages("cqrReg")
#install.packages("ggplot2")
#install.packages("stats")
#install.packages("Kendall")
#install.packages("tseries") 

library(tseries)
library(MASS)
library(Kendall)
library(ald)
library(MASS)
library(GIGrvg)
library(stats)
library(readxl)
library(forecast)
library(ggplot2)
library(ggpubr)
library(forecast)
library(qicharts)

####################========== Main Function
CQRQA=function(K,y,p,n,n.burn,n.sample){
sample=n.burn+n.sample
pp  = p + 1
Hfi = matrix(0, ncol = pp , nrow = sample)
####-----------------
Y=matrix(0,ncol= p ,nrow=n)
for(JJ in pp:n){
BB= c()
for(jjj in 1:p){
BB[jjj] = y[JJ - jjj]
}
Y[JJ,]=c(  BB )
}
Y = cbind(1 , Y)
####-----------------
####-----------------
####------------------------ generate tua and theta1 and theta2 and intials
tau = (1:K)/(K+1)
theta1 = (1-2*tau)/(tau*(1-tau))
theta2 = 2 / (tau*(1-tau))
Sigma = rep(0.1 , K)
Fi = solve(t(Y[-c(1:p),])%*%Y[-c(1:p),])%*%t(Y[-c(1:p),])%*%(y[-c(1:p)]) 
Fi = matrix(c(Fi)  , ncol=1)
alpha = rep( 0.1 , K )
####------------------------ 
for(i in 1:sample){
  ####------------------------ Genarating Latent Varibles
  v= matrix(0 , ncol= K , nrow = n )
    for(Kk in pp:n){
    v1 = c()
       #---- 
       for(KK in 1:K){
       mu = y[Kk] - (Y[Kk,] %*% Fi)-alpha[KK]
       S_ = theta2[KK] * Sigma[KK]
       S1_ = ( theta1[KK]^(2) ) + ( 2 * theta2[KK] )
       v1[KK] = rgig(1 , lambda=0.5 , chi= ( ( mu^(2) ) / S_)  , psi= (S1_ / S_) )
       }
       #----
    v[Kk,]=v1
    }
####------------------ Estimate Sigma
  Sigma=c()
    for(M in 1:K){
    vv=v[ ,M]
    Delta1= c()
    Delta1[1:p]=0
      #---- 
      for(MM in pp:n){
      muu = y[MM]- (Y[MM,] %*% Fi) - alpha[M]
      SS = 2 * theta2[M]
      SS. = ( theta1[M]^(2) ) + ( 2 * theta2[M] )
      Delta1[MM]= ( (muu^(2)/SS)* ( vv[MM]^(-1) ) ) + ((SS./SS)*vv[MM]) - ( (theta1[M]*muu) / theta2[M] ) 
      }
      #---- 
   
    Sigma[M]= (2/ (3*(n-p)) ) * sum(Delta1)
    }
####------------------ Estimate Alpha
  alpha= c()
  for(hh in 1:K){
  vv.=v[ ,hh]
  Alph=c()
  Alph[1:p]=0
    
    #----
    for(Hh in pp:n){
    Alph[Hh]= ( y[Hh] - (Y[Hh , ] %*% Fi ) ) * vv.[Hh]
    }
    #---- 
  
  alpha[hh]= ( sum(Alph) - ((n-p)*theta1[hh]) )  / ( sum(vv.^(-1) ) )   
  }
####------------------ Estimate Fi
  d= c()
  w= c()
  d[1:p]= 0
  w[1:p]= 0
  for(ll in pp:n){
  vv..=v[ll, ]
  d.= c()
  w.= c()
    #----
    for(lll in 1:K){
    d.[lll]= (theta1[lll] + alpha[lll] * ( vv..[lll]^(-1) )) / ( theta2[lll] * Sigma[lll] )
    w.[lll]= ( vv..[lll]^(-1) ) / ( theta2[lll] * Sigma[lll] )
    }
    #----
  d[ll]= sum(d.)
  w[ll]= sum(w.)
  }

w = w[-c(1:p)]
d = d[-c(1:p)]
D = matrix(d , ncol = 1 , nrow = n-p )
W = diag( w )
Y_= matrix(y[-c(1:p)] , ncol= 1 , nrow = n-p)
YY = Y[-c(1:p) , ]
Y_Bar = Y_ - ginv(W) %*% D  
Fi =ginv( t(YY) %*% W %*% YY ) %*% t(YY) %*% W %*% Y_Bar

fi=matrix(Fi,nrow=1)
Hfi[i,]=fi
}
HFI=Hfi[(n.burn+1):sample,]
return(list(HFI=HFI))
}


EMAR <- function(K, y, p, max_iter = 1000, tol = 0.01) {
  
  n <- length(y)
  m <- n - p  
  
  X <- matrix(1, nrow = m, ncol = p + 1)
  for (i in 1:p) {
    X[, i + 1] <- y[(p - i + 1):(n - i)]
  }
  MATRIXE <- t(X)
  Y_target <- matrix(y[(p + 1):n], ncol = 1)
  
  tau <- (1:K) / (K + 1)
  theta1 <- (1 - 2 * tau) / (tau * (1 - tau))
  theta2 <- 2 / (tau * (1 - tau))
  
  ALPHA <- rep(0.1, K)
  Sigma <- rep(0.5, K)
  
  ee <- arima(y, order = c(p, 0, 0))
  FI <- c(ee$coef["intercept"], ee$coef[paste0("ar", 1:p)])
  FI[is.na(FI)] <- 0 
  
  for (iter in 1:max_iter) {
    HHFI <- FI
    
    residuals_vec <- Y_target - X %*% FI
    
    alpha <- numeric(K)
    for (j in 1:K) {
      diff_res <- residuals_vec - ALPHA[j]
      delta2 <- sqrt(theta1[j]^2 + 2 * theta2[j]) / abs(diff_res)
      sum_term <- sum(residuals_vec * delta2) - (m * theta1[j]) 
      alpha[j] <- sum_term / sum(delta2)
    }
    
    A <- mean(alpha)
    ALPHA <- alpha - A
    
    for (j in 1:K) {
      diff_res <- residuals_vec - ALPHA[j]
      delta2 <- sqrt(theta1[j]^2 + 2 * theta2[j]) / abs(diff_res)
      delta3 <- (theta2[j] * Sigma[j]) / (theta1[j]^2 + 2 * theta2[j]) +
        abs(diff_res) / sqrt(theta1[j]^2 + 2 * theta2[j])
      
      DELTAA <- ((diff_res^2 * delta2) / (2 * theta2[j])) +
        (((theta1[j]^2 + 2 * theta2[j]) * delta3) / (2 * theta2[j])) -
        ((theta1[j] * diff_res) / theta2[j])
      
      Sigma[j] <- (2 / (3 * m)) * sum(DELTAA)
    }
    
    W <- numeric(m)
    D <- numeric(m)
    
    for (i in 1:m) {
      diff_res <- residuals_vec[i]
      
      Wt_k <- numeric(K)
      DD1_k <- numeric(K)
      
      for (k in 1:K) {
        delta222 <- sqrt(theta1[k]^2 + 2 * theta2[k]) / abs(diff_res - ALPHA[k])
        Wt_k[k] <- delta222 / (theta2[k] * Sigma[k])
        DD1_k[k] <- (theta1[k] + ALPHA[k] * delta222) / (theta2[k] * Sigma[k])
      }
      
      W[i] <- sum(Wt_k)
      D[i] <- sum(DD1_k)
    }
    
    MATRIXW <- diag(W)
    MATRIXinvers <- ginv(MATRIXE %*% MATRIXW %*% t(MATRIXE))
    
    Ybar <- Y_target - (D / W)
    
    FI <- MATRIXinvers %*% MATRIXE %*% (MATRIXW %*% Ybar)
    FI <- as.vector(FI)
    
    if (max(abs(FI - HHFI)) <= tol) {
      break
    }
  }
  
  return(list(FI = FI, Iterations = iter))
}






####--------------- MAP WCQAR
MAP_WCQAR=function(y. , m , p , h ){
  y = y.[c(1:m)] ; pp = p + 1    ; hh = h + p
  fi_ =   CQRQA( K =9 , y = y , p = p ,  n=length(y)  , n.burn=2000 ,  n.sample=2000)
  f_  = c( mean(fi_$HFI[,1]) , mean(fi_$HFI[,2]) )

   
    e1=c() ; y_hat=c()
    e1[1:p]  = 0
    y_hat[1] = y.[m]
      for(i in pp:hh){
      y_hat[i] =  f_[1] +(y_hat[i-1] * f_[2])
      e1[i]    =  abs(y.[i+m-p] - y_hat[i])
      }
E=e1[-c(1:p)]
return(list(E=E, yhat = y_hat[-c(1)] ))
}

####--------------- MAP WCQAREM
MAP_WCQAREM <- function(y. , m , p , h ){
  y <- y.[c(1:m)] ; pp = p + 1    ; hh = h + p
  fi_ <-  EMAR(K =9, y = y, p = p, max_iter = 1000, tol = 0.01) 
  f_  <- fi_$FI

   
  e1=c() ; y_hat=c()
  e1[1:p]  = 0
  y_hat[1] = y.[m] 
  for(i in pp:hh){
    y_hat[i] =  f_[1] +(y_hat[i-1] * f_[2]) 
    e1[i]    =  abs(y.[i+m-p] - y_hat[i])
  }
  E=e1[-c(1:p)]
  return(list(E=E, yhat = y_hat[-c(1)] ))
}


#######################
####################### Start Pridiction
#######################

data <- read_excel("E:/My Articles/New Article 1/data 3.xlsx")
y <- as.numeric(data$Data) 
p <- 1




n_start = 38
h  = 2
h. = 6
p  = 1



Mat_MAP_OLS = matrix( 0 , ncol = h , nrow = h. )
y.hat1 <- c()

  for(g in 1:h.){
    nn=( g * h) + ( n_start - h)
    fi_ <- ar.ols(y[1:nn], FALSE, p)
    Predict = predict(fi_, n.ahead = h)
    Error = abs(as.numeric(Predict$pred) - y[(nn+1) : (nn+h)])
    y.hat1 <- c(y.hat1, as.numeric(Predict$pred))
    Mat_MAP_OLS[g,] = Error
  }
 


Mat_MAP_MLE = matrix( 0 , ncol = h , nrow = h. )
y.hat2 <- c()

  for(j in 1:h.){
    nn=( j * h) + ( n_start - h)
    fit_ar = arima(y[1:nn] , order=c(p , 0 , 0 )) #ar(y[1:nn], order.max = p, method = "mle")
    pred <- predict(fit_ar, n.ahead = h)
    Pred <- as.numeric(pred$pred)
    Error = abs(Pred - y[(nn+1) : (nn+h)])
    y.hat2 <- c(y.hat2, Pred)
    Mat_MAP_MLE[j,] = Error
  }
 


Mat_MAP_WCQAR = matrix( 0 , ncol = h , nrow = h. )
y.hat3 <- c()

  for(jj in 1:h.){
    nn=( jj * h) + ( n_start - h)
    Predict = MAP_WCQAR(y. = y , m = nn , p = p ,  h = h)
    Error = Predict$E
    y.hat3 <- c(y.hat3, Predict$yhat)
    Mat_MAP_WCQAR[jj,] = Error
print(jj)   
  }


Mat_MAP_WCQAREM = matrix( 0 , ncol = h , nrow = h. )
y.hat4 <- c()

  for(jj in 1:h.){
    nn=( jj * h) + ( n_start - h)
    Predict = MAP_WCQAREM(y. = y , m = nn , p = p ,  h = h)
    Error = Predict$E
    y.hat4 <- c(y.hat4, Predict$yhat)
    Mat_MAP_WCQAREM[jj,] = Error
  print(jj)   
  }



mean(Mat_MAP_OLS)
sd(Mat_MAP_OLS)

mean(Mat_MAP_MLE)
sd(Mat_MAP_MLE)

mean(Mat_MAP_WCQAR)
sd(Mat_MAP_WCQAR)

mean(Mat_MAP_WCQAREM)
sd(Mat_MAP_WCQAREM)


























###########################################
####################========== Packages

# Install packages if needed
# install.packages(c("ald", "MASS", "VGAM", "quantreg", "GIGrvg", 
#                    "cqrReg", "ggplot2", "stats", "Kendall", "tseries", 
#                    "readxl", "forecast", "ggpubr", "qicharts"))

library(tseries)
library(MASS)
library(Kendall)
library(ald)
library(GIGrvg)
library(readxl)
library(forecast)
library(ggplot2)
library(ggpubr)
library(qicharts)

####################========== Main Function

CQRQA <- function(K, y, p, n, n.burn, n.sample) {
  sample_size <- n.burn + n.sample
  pp <- p + 1
  Hfi <- matrix(0, ncol = pp, nrow = sample_size)
  
  # Generate Y matrix
  Y <- matrix(0, ncol = p, nrow = n)
  for (JJ in pp:n) {
    BB <- c()
    for (jjj in 1:p) {
      BB[jjj] <- y[JJ - jjj]
    }
    Y[JJ, ] <- BB
  }
  Y <- cbind(1, Y)
  
  # Generate tau, theta1, theta2 and initials
  tau <- (1:K) / (K + 1)
  theta1 <- (1 - 2 * tau) / (tau * (1 - tau))
  theta2 <- 2 / (tau * (1 - tau))
  Sigma <- rep(0.1, K)
  
  Fi <- solve(t(Y[-c(1:p), ]) %*% Y[-c(1:p), ]) %*% t(Y[-c(1:p), ]) %*% y[-c(1:p)]
  Fi <- matrix(c(Fi), ncol = 1)
  alpha <- rep(0.1, K)
  
  for (i in 1:sample_size) {
    
    # Generating Latent Variables
    v <- matrix(0, ncol = K, nrow = n)
    for (Kk in pp:n) {
      v1 <- c()
      for (KK in 1:K) {
        mu <- y[Kk] - (Y[Kk, ] %*% Fi) - alpha[KK]
        S_ <- theta2[KK] * Sigma[KK]
        S1_ <- (theta1[KK]^2) + (2 * theta2[KK])
        v1[KK] <- rgig(1, lambda = 0.5, chi = (mu^2 / S_), psi = (S1_ / S_))
      }
      v[Kk, ] <- v1
    }
    
    # Estimate Sigma
    Sigma <- c()
    for (M in 1:K) {
      vv <- v[, M]
      Delta1 <- c()
      Delta1[1:p] <- 0
      for (MM in pp:n) {
        muu <- y[MM] - (Y[MM, ] %*% Fi) - alpha[M]
        SS <- 2 * theta2[M]
        SS. <- (theta1[M]^2) + (2 * theta2[M])
        Delta1[MM] <- ((muu^2 / SS) * (vv[MM]^(-1))) + ((SS. / SS) * vv[MM]) - ((theta1[M] * muu) / theta2[M])
      }
      Sigma[M] <- (2 / (3 * (n - p))) * sum(Delta1)
    }
    
    # Estimate Alpha
    alpha <- c()
    for (hh in 1:K) {
      vv. <- v[, hh]
      Alph <- c()
      Alph[1:p] <- 0
      for (Hh in pp:n) {
        Alph[Hh] <- (y[Hh] - (Y[Hh, ] %*% Fi)) * vv.[Hh]
      }
      alpha[hh] <- (sum(Alph) - ((n - p) * theta1[hh])) / sum(vv.^(-1))
    }
    
    # Estimate Fi
    d <- c()
    w <- c()
    d[1:p] <- 0
    w[1:p] <- 0
    for (ll in pp:n) {
      vv.. <- v[ll, ]
      d. <- c()
      w. <- c()
      for (lll in 1:K) {
        d.[lll] <- (theta1[lll] + alpha[lll] * (vv..[lll]^(-1))) / (theta2[lll] * Sigma[lll])
        w.[lll] <- (vv..[lll]^(-1)) / (theta2[lll] * Sigma[lll])
      }
      d[ll] <- sum(d.)
      w[ll] <- sum(w.)
    }
    
    w <- w[-c(1:p)]
    d <- d[-c(1:p)]
    D <- matrix(d, ncol = 1, nrow = n - p)
    W <- diag(w)
    Y_ <- matrix(y[-c(1:p)], ncol = 1, nrow = n - p)
    YY <- Y[-c(1:p), ]
    Y_Bar <- Y_ - ginv(W) %*% D
    Fi <- ginv(t(YY) %*% W %*% YY) %*% t(YY) %*% W %*% Y_Bar
    
    fi <- matrix(Fi, nrow = 1)
    Hfi[i, ] <- fi
  }
  
  HFI <- Hfi[(n.burn + 1):sample_size, ]
  return(list(HFI = HFI))
}


EMAR <- function(K, y, p, max_iter = 1000, tol = 0.01) {
  n <- length(y)
  m <- n - p
  
  X <- matrix(1, nrow = m, ncol = p + 1)
  for (i in 1:p) {
    X[, i + 1] <- y[(p - i + 1):(n - i)]
  }
  
  MATRIXE <- t(X)
  Y_target <- matrix(y[(p + 1):n], ncol = 1)
  
  tau <- (1:K) / (K + 1)
  theta1 <- (1 - 2 * tau) / (tau * (1 - tau))
  theta2 <- 2 / (tau * (1 - tau))
  
  ALPHA <- rep(0.1, K)
  Sigma <- rep(0.5, K)
  
  ee <- arima(y, order = c(p, 0, 0))
  FI <- c(ee$coef["intercept"], ee$coef[paste0("ar", 1:p)])
  FI[is.na(FI)] <- 0
  
  for (iter in 1:max_iter) {
    HHFI <- FI
    residuals_vec <- Y_target - X %*% FI
    
    # Estimate alpha
    alpha <- numeric(K)
    for (j in 1:K) {
      diff_res <- residuals_vec - ALPHA[j]
      delta2 <- sqrt(theta1[j]^2 + 2 * theta2[j]) / abs(diff_res)
      sum_term <- sum(residuals_vec * delta2) - (m * theta1[j])
      alpha[j] <- sum_term / sum(delta2)
    }
    
    A <- mean(alpha)
    ALPHA <- alpha - A
    
    # Estimate Sigma
    for (j in 1:K) {
      diff_res <- residuals_vec - ALPHA[j]
      delta2 <- sqrt(theta1[j]^2 + 2 * theta2[j]) / abs(diff_res)
      delta3 <- (theta2[j] * Sigma[j]) / (theta1[j]^2 + 2 * theta2[j]) +
        abs(diff_res) / sqrt(theta1[j]^2 + 2 * theta2[j])
      
      DELTAA <- ((diff_res^2 * delta2) / (2 * theta2[j])) +
        (((theta1[j]^2 + 2 * theta2[j]) * delta3) / (2 * theta2[j])) -
        ((theta1[j] * diff_res) / theta2[j])
      
      Sigma[j] <- (2 / (3 * m)) * sum(DELTAA)
    }
    
    # Estimate Fi
    W <- numeric(m)
    D <- numeric(m)
    
    for (i in 1:m) {
      diff_res <- residuals_vec[i]
      Wt_k <- numeric(K)
      DD1_k <- numeric(K)
      
      for (k in 1:K) {
        delta222 <- sqrt(theta1[k]^2 + 2 * theta2[k]) / abs(diff_res - ALPHA[k])
        Wt_k[k] <- delta222 / (theta2[k] * Sigma[k])
        DD1_k[k] <- (theta1[k] + ALPHA[k] * delta222) / (theta2[k] * Sigma[k])
      }
      
      W[i] <- sum(Wt_k)
      D[i] <- sum(DD1_k)
    }
    
    MATRIXW <- diag(W)
    MATRIXinvers <- ginv(MATRIXE %*% MATRIXW %*% t(MATRIXE))
    Ybar <- Y_target - (D / W)
    
    FI <- MATRIXinvers %*% MATRIXE %*% (MATRIXW %*% Ybar)
    FI <- as.vector(FI)
    
    if (max(abs(FI - HHFI)) <= tol) {
      break
    }
  }
  
  return(list(FI = FI, Iterations = iter))
}


####--------------- MAP WCQAR
MAP_WCQAR <- function(y., m, p, h) {
  y <- y.[c(1:m)]
  pp <- p + 1
  hh <- h + p
  
  fi_ <- CQRQA(K = 9, y = y, p = p, n = length(y), n.burn = 2000, n.sample = 2000)
  f_ <- c(mean(fi_$HFI[, 1]), mean(fi_$HFI[, 2]))
  
  e1 <- c()
  y_hat <- c()
  e1[1:p] <- 0
  y_hat[1] <- y.[m]
  
  for (i in pp:hh) {
    y_hat[i] <- f_[1] + (y_hat[i - 1] * f_[2])
    e1[i] <- abs(y.[i + m - p] - y_hat[i])
  }
  
  E <- e1[-c(1:p)]
  return(list(E = E, yhat = y_hat[-c(1)]))
}


####--------------- MAP WCQAREM
MAP_WCQAREM <- function(y., m, p, h) {
  y <- y.[c(1:m)]
  pp <- p + 1
  hh <- h + p
  
  fi_ <- EMAR(K = 9, y = y, p = p, max_iter = 1000, tol = 0.01)
  f_ <- fi_$FI
  
  e1 <- c()
  y_hat <- c()
  e1[1:p] <- 0
  y_hat[1] <- y.[m]
  
  for (i in pp:hh) {
    y_hat[i] <- f_[1] + (y_hat[i - 1] * f_[2])
    e1[i] <- abs(y.[i + m - p] - y_hat[i])
  }
  
  E <- e1[-c(1:p)]
  return(list(E = E, yhat = y_hat[-c(1)]))
}


#######################
####################### Start Prediction
#######################

data <- read_excel("E:/My Articles/New Article 1/data 3.xlsx")
y <- as.numeric(data$Data)

p <- 1
n_start <- 38
h <- 4
h. <- 3


# 1. OLS Model
Mat_MAP_OLS <- matrix(0, ncol = h, nrow = h.)
y.hat1 <- c()

for (g in 1:h.) {
  nn <- (g * h) + (n_start - h)
  fi_ <- ar.ols(y[1:nn], FALSE, p)
  Predict <- predict(fi_, n.ahead = h)
  Error <- abs(as.numeric(Predict$pred) - y[(nn + 1):(nn + h)])
  y.hat1 <- c(y.hat1, as.numeric(Predict$pred))
  Mat_MAP_OLS[g, ] <- Error
}

# 2. MLE Model
Mat_MAP_MLE <- matrix(0, ncol = h, nrow = h.)
y.hat2 <- c()

for (j in 1:h.) {
  nn <- (j * h) + (n_start - h)
  fit_ar <- arima(y[1:nn], order = c(p, 0, 0))
  pred <- predict(fit_ar, n.ahead = h)
  Pred <- as.numeric(pred$pred)
  Error <- abs(Pred - y[(nn + 1):(nn + h)])
  y.hat2 <- c(y.hat2, Pred)
  Mat_MAP_MLE[j, ] <- Error
}

# 3. WCQAR Model
Mat_MAP_WCQAR <- matrix(0, ncol = h, nrow = h.)
y.hat3 <- c()

for (jj in 1:h.) {
  nn <- (jj * h) + (n_start - h)
  Predict <- MAP_WCQAR(y. = y, m = nn, p = p, h = h)
  Error <- Predict$E
  y.hat3 <- c(y.hat3, Predict$yhat)
  Mat_MAP_WCQAR[jj, ] <- Error
  print(jj)
}

# 4. WCQAREM Model
Mat_MAP_WCQAREM <- matrix(0, ncol = h, nrow = h.)
y.hat4 <- c()

for (jj in 1:h.) {
  nn <- (jj * h) + (n_start - h)
  Predict <- MAP_WCQAREM(y. = y, m = nn, p = p, h = h)
  Error <- Predict$E
  y.hat4 <- c(y.hat4, Predict$yhat)
  Mat_MAP_WCQAREM[jj, ] <- Error
  print(jj)
}


####################### Results Summary
cat("### OLS Model ###\n")
cat("Mean Error:", mean(Mat_MAP_OLS), "\n")
cat("SD Error:", sd(Mat_MAP_OLS), "\n\n")

cat("### MLE Model ###\n")
cat("Mean Error:", mean(Mat_MAP_MLE), "\n")
cat("SD Error:", sd(Mat_MAP_MLE), "\n\n")

cat("### WCQAR Model ###\n")
cat("Mean Error:", mean(Mat_MAP_WCQAR), "\n")
cat("SD Error:", sd(Mat_MAP_WCQAR), "\n\n")

cat("### WCQAREM Model ###\n")
cat("Mean Error:", mean(Mat_MAP_WCQAREM), "\n")
cat("SD Error:", sd(Mat_MAP_WCQAREM), "\n")







