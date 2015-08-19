
####################################################################################
###  Simulation-based estimation of the marginal posteriors for
###  the unknown parameters of a bivariate normal linear regression
###  model with:
###
###              regression (slope) coefficient = beta2;
###              error standard deviation = sig
###
###  using Gibbs sampling
###
###  Comparison made with the analytical marginals (available in this case)
###  
####################################################################################



####################################################################################
###########     Procedures to be called in program below      ######################
###########    Need to be at the top of the program in R      ######################
####################################################################################


## Student t marginal posterior probability density function (pdf) of beta2  ##


margbeta2 = function(beta2){
  
  t = (beta2 - beta2hat)/(sighat*sqrt(x22))
  
  intcon = (gamma((1+df)/2)/(gamma(df/2)*sqrt(df*pi)))*
    (sighat*sqrt(x22))^(-1)
  
  studt = intcon*(1 + (1/df)*t^2)^(-(1+df)/2)
  
  return(studt)
  
}


## Inverted gamma marginal posterior probability density function (pdf) of sig ##


margsig = function(sig){
  
  intcon = (2/(gamma(df/2)))*(df*(sighat^2)/2)^(df/2)
  
  invgamma = intcon*(sig^(-(1+df)))*exp(-(df*sighat^2)/(2*sig^2))
  
  return(invgamma)
  
}


## Normal conditional posterior probability density function (pdf) of beta2  ##


condbeta2 = function(beta2,varbeta2){
  
  z = (beta2 - beta2hat)/sqrt(varbeta2)
  
  cond_pdf = ((2*pi*varbeta2)^(-1/2))*exp(-0.5*(z^2))
  
  return(cond_pdf)
  
}


## Inverted gamma conditional posterior probability density function (pdf) of sig  ##


conds = function(sig,a){
  
  intcon = (2/(gamma(nobs/2)))*(a/2)^(nobs/2)
  
  condfunc = intcon*(sig^(-(1+nobs)))*exp(-a/(2*sig^2))
  
  return(condfunc)
  
}


###########################################################
### Setting the grids for the marginal density plots
###########################################################


lbeta2 = 0.1
upbeta2 = 1.0
beta2grid = 0.01
numbeta2 = round((upbeta2-lbeta2)/beta2grid)

lsig = 0.1
upsig = 1.0
siggrid = 0.01
numsig = round((upsig-lsig)/siggrid)


###########################################################
########       Generation of data = y               #######
###########################################################


true_beta = matrix(c(2,0.5,0.7),3,1)      ## True value of beta (vector) used in the data generating process (dgp)

true_sig = 0.3                            ## True value of sig used in the data generating process (dgp)



## Setting the random number seed

set.seed(123456, kind = NULL, normal.kind = NULL)

## Number of observations

nobs=20
##nobs =200


x1 = matrix(1,nobs,1)
x2 = matrix(c(rnorm(nobs,0,1)),nobs,1)
x3 = matrix(c(rnorm(nobs,0,1)),nobs,1)
u  = matrix(c(rnorm(nobs,0,1)),nobs,1)%*%true_sig

x = cbind(x1,x2,x3)

y = x%*%true_beta + u


###########################################################
### Specification of sample statistics
###########################################################

invxx = solve(t(x)%*%x)
x22 = invxx[2,2]
betahat = invxx%*%t(x)%*%y

df = nobs - 3
beta2hat = betahat[2]
sighat = sqrt((t(y-x%*%betahat)%*%(y-x%*%betahat))/df)


###########################################################
###
###  Specification of the exact marginal posterior
###  Student t pdf for beta2
###  based on the noninformative Jeffreys prior
###  with df = nobs - 3 degrees of freedon
###
###  (E.g. Judge et al, 2nd ed. Chp 8)
###
###########################################################


beta2 = seq(lbeta2,upbeta2,by=beta2grid)
exactbeta2 = margbeta2(beta2)           ## Calling procedure that specificies marginal posterior pdf of beta2



###########################################################
###
###  Specification of the exact marginal posterior
###  inverted gamma pdf for sigma
###  based on the noninformative Jeffreys prior
###  with df = nobs - 3  degrees of freedon
###
###  (E.g. Judge et al, 2nd ed. Chp 8)
###
###########################################################


sig = seq(lsig,upsig,by=siggrid)
exactsig = margsig(sig)                 ## Calling procedure that specificies marginal posterior of sigma



###############################################################
###
###  Generation of draws of beta2 and sigma via Gibbs sampling
###
###############################################################

repl = 100    ## Number of Gibbs replications
##repl = 1000
##repl = 5000



betav = matrix(0,repl,3)     ## Dimensioning matrix of beta draws ##  
sigv = matrix(0,repl,1)

##con_density_beta2 = matrix(0,1,(numbeta2+1))
##con_density_sigma = matrix(0,1,(numsig+1))

cbeta2 = matrix(0,1,(numbeta2+1))
csig = matrix(0,1,(numsig+1))

stnbeta = matrix(c(rnorm((repl*3),0,1)),repl,3)  ## Generating underlying normal random draws to be used in Gibbs algorithm


sigv[1] = true_sig    ## Specifying starting value for Gibbs chain


for(j in 2:repl){
  
  ## Generation of beta (vector) via its known multivariate normal conditional
  
  meanbeta = betahat
  varbeta = (sigv[j-1])^2*invxx
    
  varbeta2=(sigv[j-1])^2%*%x22    ## Specification of variance of conditional normal posterior of beta 2

  
  betagen = chol(varbeta)%*%stnbeta[j,] + meanbeta
  betav[j,] = betagen
  
  
  ## Generation of sigma via its known inverted gamma conditional
  
  a = t(y-x%*%betav[j,])%*%(y-x%*%betav[j,])
  
  zvec = matrix(c(rnorm(nobs,0,1)),nobs,1)
  
  chisq = t(zvec)%*%zvec
  
  sigv[j] = sqrt(a/chisq)
  
  
  ##  Update standard deviation of conditional density of beta2
  
  varbeta2 = (sigv[j])^2%*%x22
  
  ## Evaluate the conditional density of beta2
  
  con_density_beta2 = condbeta2(beta2,varbeta2)
  
  ## Cummulate the nromal conditionals
  
  cbeta2 = cbeta2 + con_density_beta2
  
  
  ## Evaluate the conditional density of sigma
  
  con_density_sigma = conds(sig,a)
  
  ## Cummulate the inverted gamm conditionals
  
  csig = csig + con_density_sigma
  
  ## Plot selected conditional posteriors of beta2
  
  if(j==4){
    
  condbeta24=con_density_beta2
  layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
  par(mar=c(5,5,5,1))
  plot(beta2,con_density_beta2,lty=1,type='l',
       xlab=expression(beta[2]),ylab=expression(paste("p(",beta[2],"|",sigma,",y)"))
  )
  title(bquote(paste("Conditional posterior for ",beta[2],";          ",sigma^{(j)}," = ",.(round(sigv[j-1],digits=4)))), line=4) 
  title(bquote(paste("E(",beta[2],"|",sigma^{(j)},",y) = ",.(round(beta2hat,digits=4)),"       ",
                     "var(",beta[2],"|",sigma^{(j)},",y) =",.(round((sigv[j-1])^2%*%x22,digits=4)))), line=2)
  
  layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
  par(mar=c(5,5,5,1))
  plot(sig,con_density_sigma,lty=1,type='l',
       xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",beta[2],",y)"))
  )
  title(bquote(paste("Conditional posterior for ",sigma)), line=2)
  
  }
  
  if(j==5){
    
    condbeta25=con_density_beta2
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(beta2,con_density_beta2,lty=1,type='l',
         xlab=expression(beta[2]),ylab=expression(paste("p(",beta[2],"|",sigma,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",beta[2],";          ",sigma^{(j)}," = ",.(round(sigv[j-1],digits=4)))), line=4) 
    title(bquote(paste("E(",beta[2],"|",sigma^{(j)},",y) = ",.(round(beta2hat,digits=4)),"       ",
                       "var(",beta[2],"|",sigma^{(j)},",y) =",.(round((sigv[j-1])^2%*%x22,digits=4)))), line=2)
    
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(sig,con_density_sigma,lty=1,type='l',
         xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",beta[2],",y)"))
    )
    title(bquote(paste("Conditional posterior for ",sigma)), line=2)
    
  }
  
  if(j==6){
    
    condbeta26=con_density_beta2
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(beta2,con_density_beta2,lty=1,type='l',
         xlab=expression(beta[2]),ylab=expression(paste("p(",beta[2],"|",sigma,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",beta[2],";          ",sigma^{(j)}," = ",.(round(sigv[j-1],digits=4)))), line=4) 
    title(bquote(paste("E(",beta[2],"|",sigma^{(j)},",y) = ",.(round(beta2hat,digits=4)),"       ",
                       "var(",beta[2],"|",sigma^{(j)},",y) =",.(round((sigv[j-1])^2%*%x22,digits=4)))), line=2)
    
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(sig,con_density_sigma,lty=1,type='l',
         xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",beta[2],",y)"))
    )
    title(bquote(paste("Conditional posterior for ",sigma)), line=2)
    
  }
  
  if(j==7){
    
    condbeta27=con_density_beta2
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(beta2,con_density_beta2,lty=1,type='l',
         xlab=expression(beta[2]),ylab=expression(paste("p(",beta[2],"|",sigma,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",beta[2],";          ",sigma^{(j)}," = ",.(round(sigv[j-1],digits=4)))), line=4) 
    title(bquote(paste("E(",beta[2],"|",sigma^{(j)},",y) = ",.(round(beta2hat,digits=4)),"       ",
                       "var(",beta[2],"|",sigma^{(j)},",y) =",.(round((sigv[j-1])^2%*%x22,digits=4)))), line=2)
    
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(sig,con_density_sigma,lty=1,type='l',
         xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",beta[2],",y)"))
    )
    title(bquote(paste("Conditional posterior for ",sigma)), line=2)
    
  }
  
}

############################################################
###
###  Estimation of the marginal density of beta2 
###  as an average of the stored conditionals
###
############################################################


ave_con_density_beta2 = cbeta2/repl

ave_con_density_sigma = csig/repl


## Remove the starting value of beta2 for demonstration of histogram-based estimate

beta2v = betav[2:nrow(betav),2]




## Histogram of draws of beta2

hist(beta2v,30,ylab="Percentage frequencies",xlab=expression(beta[2]),
     main=expression(paste("Histogram of draws of ",beta[2])))


## kernel density plot of the draws of beta2

d=density(beta2v)
plot(d,main=expression(paste("Kernel density estimate of marginal of ",beta[2])),
     xlab=expression(beta[2]),ylab=expression(paste("p(",beta[2],"|y)")),xlim=c(lbeta2,upbeta2))


###############################################################################################
###
###  Plot of the selected conditional pdfs of beta2 against the average of all repl conditionals
###
###############################################################################################


layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
par(mar=c(5,5,2,1))
plot(beta2,ave_con_density_beta2,lty=2,type='l',
     main=expression(paste("Conditional posteriors for ",beta[2]," and average of all 100 conditionals")),
     xlab=expression(beta[2]),ylab=expression(paste("p(",beta[2],"|",sigma,",y)")),
     ylim=c(min(condbeta24,condbeta25,condbeta26,condbeta27,ave_con_density_beta2),
            max(condbeta24,condbeta25,condbeta26,condbeta27,ave_con_density_beta2)))
lines(beta2,condbeta24,type='l',lty=3)
lines(beta2,condbeta25,type='l',lty=4)
lines(beta2,condbeta26,type='l',lty=5)
lines(beta2,condbeta27,type='l',lty=6)
legend(0.7, max(condbeta24,condbeta25,condbeta26,condbeta27,ave_con_density_beta2), c("Conditional 1","Conditional 2","Conditional 3","Conditional 4","Average of conditionals"),
       lty = c(3,4,5,6,2),  merge = TRUE,cex=0.8)


################################################################################################################
###
###  Plot of the estimated marginal pdf (average of all repl conditionals) of mu against the exact marginal pdf
###
################################################################################################################


#beta2

layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
par(mar=c(5,5,2,1))
plot(beta2,ave_con_density_beta2,lty=2,type='l',main=expression(paste("Average of all 100 conditionals & exact marginal posterior for ",beta[2])),xlab=expression(beta[2]),ylab=expression(paste("p(",beta[2],"|y)")),
     ylim=c(min(ave_con_density_beta2,exactbeta2),max(ave_con_density_beta2,exactbeta2)),
     xlim=c(lbeta2,upbeta2))
lines(beta2,exactbeta2,type='l',lty=1)
legend(0.65, max(ave_con_density_beta2,exactbeta2), c("Average of all conditionals", "Exact marginal"),
       lty = c(2, 1),  merge = TRUE,cex=0.8)

#sigma

layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
par(mar=c(5,5,2,1))
plot(sig,ave_con_density_sigma,lty=1,type='l',main=expression(paste("Average of all 100 conditionals & exact marginal posterior for ",sigma)),xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|y)")),
     ylim=c(min(ave_con_density_sigma,exactsig),max(ave_con_density_sigma,exactsig)))
lines(sig,exactsig,type='l',lty=2)
legend(0.65, max(ave_con_density_sigma,exactsig), c("Average of all conditionals", "Exact marginal"),
       lty = c(2, 1),  merge = TRUE,cex=0.8)


