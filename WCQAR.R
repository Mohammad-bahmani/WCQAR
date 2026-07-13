###########################################

####################========== Packages
#install.packages("MASS")
#install.packages("GIGrvg")
#install.packages("stats")



library(MASS)
library(GIGrvg)
library(stats)


####################========== Main Function
CQRQA=function(K,y,p,n,n.burn,n.sample){
sample=n.burn+n.sample
pp  = p + 1
Hfi = matrix(0, ncol = p , nrow = sample)
####-----------------
Y=matrix(0,ncol= p ,nrow=n)
for(JJ in pp:n){
BB= c()
for(jjj in 1:p){
BB[jjj] = y[JJ - jjj]
}
Y[JJ,]=c(  BB )
}
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
