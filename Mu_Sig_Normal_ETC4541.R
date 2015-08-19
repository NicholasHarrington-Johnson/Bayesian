
####################################################################################
###
###  Simulation-based estimation of the marginal posteriors for
###  the unknown parameters:
###
###              mean = mu
###              standard deviation = sig
###
###  of a normal random variable, using Gibbs sampling
###
###  Comparison made with the analytical marginals (available in this case)
###  
####################################################################################


####################################################################################
###########     Procedures to be called in program below      ######################
###########    Need to be at the top of the program in R      ######################
####################################################################################



## Student t marginal posterior probability density function (pdf) of mu  ##


margmu = function(mu){
  
  t = (mu - ybar)/(sighat/sqrt(nobs))
  
  intcon = (gamma((1+df)/2)/(gamma(df/2)*sqrt(df*pi)))*(sighat/sqrt(nobs))^(-1)
  
  ## Remember that gamma(1/2) = sqrt(pi)   
  
  studt = intcon*(1+(1/df)*t^2)^(-(1+df)/2)
  
  return(studt)
  
}

## Inverted gamma marginal posterior probability density function (pdf) of sig ##

margsig = function(sig){
  
  intcon = (2/(gamma(df/2)))*(df*(sighat^2)/2)^(df/2)
  
  invgamma = intcon*(sig^(-(1+df)))*exp(-(df*sighat^2)/(2*sig^2))
  
  return(invgamma)
  
}


## Normal conditional posterior probability density function (pdf) of mu  ##

condmu = function(mu,varmu){
  
  z=(mu - ybar)/sqrt(varmu)
  
  cond_pdf=((2*pi*varmu)^(-1/2))*exp(-0.5*(z^2))
  
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

lmu=0.4
upmu=4.0
mugrid=0.01
nummu=round((upmu-lmu)/mugrid)


lsig=0.4
upsig=1.6
siggrid=0.01
numsig=round((upsig-lsig)/siggrid)




###########################################################
########       Generation of data = y               #######
###########################################################


true_mu = 2.0   ## True value of mu used in the data generating process (dgp)

true_sig = 1.0  ## True value of sig used in the data generating process (dgp)


## Setting the random number seed

set.seed(123458, kind = NULL, normal.kind = NULL)

## Number of observations

nobs=20
##nobs =200

y = true_mu + rnorm(nobs,0,1)*true_sig  ##  Check interpretation of 0

### print (y)

###########################################################
### Specification of sample statistics
###########################################################

ybar=mean(y)
sighat=sd(y)
df=nobs-1


###########################################################
###
###  Specification of the exact marginal posterior
###  Student t pdf for mu
###  based on the noninformative Jeffreys prior
###  with df = nobs -1  degrees of freedom
###
###  (E.g. Judge et al, 2nd ed. page 86, Eqn (5.4.4))
###
###########################################################


mu = seq(lmu,upmu,by=mugrid)

exactmu = margmu(mu)           ## Calling the procedure that specificies the marginal posterior pdf of mu

###########################################################
###
###  Specification of the exact marginal posterior
###  inverted gamma pdf for sigma
###  based on the noninformative Jeffreys prior
###  with df = nobs -1  degrees of freedon
###
###  (E.g. Judge et al, 2nd ed. page 86, Eqn (5.4.4))
###
###########################################################


sig = seq(lsig,upsig,by=siggrid)
exactsig = margsig(sig)         ## Calling the procedure that specificies marginal posterior pdf of sigma


############################################################
###
###  Generation of draws of mu and sigma via Gibbs sampling
###
############################################################

repl = 100      ## Number of Gibbs replications
##repl = 1000
##repl = 3000


muv = matrix(0,repl,1)    ## Dimensioning vector of mu draws ##  
sigv = matrix(0,repl,1)

con_density_mu = matrix(0,1,(nummu+1))     ## Dimensioning vector of ordinates of conditional pdf for mu
con_density_sig = matrix(0,1,(numsig+1))

cmu = matrix(0,1,(nummu+1))     ##  Dimensioning vector of grid values for mu
csig = matrix(0,1,(numsig+1))

stnmu = rnorm(repl,0,1)  ## Generating underlying normal random draws to be used in Gibbs algorithm

sigv[1,1] = true_sig  ## Play with this to check robustness of results to starting value of Gibbs algorithm 

for(j in 2:repl){
  
  ## Generation of mu via its known normal conditional
  
  meanmu = ybar
  
  sdmu = sigv[j-1]/sqrt(nobs)
  
  muv[j] = sdmu*stnmu[j] + meanmu
  
    
  ## Generation of sigma via its known inverted gamma conditional
  
  a = t(y-muv[j])%*%(y-muv[j])
  
  zvec = rnorm(nobs,0,1)
  
  chisq = t(zvec)%*%zvec
  
  sigv[j] = sqrt(a/chisq)
  
  ##  Update standard deviation of conditional density of mu
  
  sdmu = sigv[j]/sqrt(nobs)
  
  ## Evaluate the conditional density of mu
  
  con_density_mu = condmu(mu,(sdmu^2))
  
  ## Cummulate the normal conditionals
  
  cmu = cmu + con_density_mu
  
  ## Evaluate the conditional density of sigma
  
  con_density_sigma = conds(sig,a)
    
  ## Cummulate the inverted gamma conditionals
  
  csig = csig + con_density_sig
  
  ## Plot selected conditional posteriors of mu
  
  if(j==4){
    
    condmu4=con_density_mu
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(mu,con_density_mu,lty=1,type='l',
         xlab=expression(mu),ylab=expression(paste("p(",mu,"|",sigma,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",mu,";          ",sigma^{(j)}," = ",.(round(sigv[j],digits=4)))), line=4) 
    title(bquote(paste("E(",mu,"|",sigma^{(j)},",y) = ybar = ",.(round(ybar,digits=4)),"       ",
                       "var(",mu,"|",sigma^{(j)},",y) = ",sigma^{2(j)},"/n = ",.(round(sigv[j]^2/nobs,digits=4)))), line=2)
    
    
  }
  
  if(j==5){
    
    condmu5=con_density_mu
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(mu,con_density_mu,lty=1,type='l',
         xlab=expression(mu),ylab=expression(paste("p(",mu,"|",sigma,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",mu,";          ",sigma^{(j)}," = ",.(round(sigv[j],digits=4)))), line=4) 
    title(bquote(paste("E(",mu,"|",sigma^{(j)},",y) = ybar = ",.(round(ybar,digits=4)),"       ",
                       "var(",mu,"|",sigma^{(j)},",y) = ",sigma^{2(j)},"/n = ",.(round(sigv[j]^2/nobs,digits=4)))), line=2)
    
  }
  
  if(j==6){
    
    condmu6=con_density_mu
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(mu,con_density_mu,lty=1,type='l',
         xlab=expression(mu),ylab=expression(paste("p(",mu,"|",sigma,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",mu,";          ",sigma^{(j)}," = ",.(round(sigv[j],digits=4)))), line=4) 
    title(bquote(paste("E(",mu,"|",sigma^{(j)},",y) = ybar = ",.(round(ybar,digits=4)),"       ",
                       "var(",mu,"|",sigma^{(j)},",y) = ",sigma^{2(j)},"/n = ",.(round(sigv[j]^2/nobs,digits=4)))), line=2)
    
  }
  
  if(j==7){
    
    condmu7=con_density_mu
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(mu,con_density_mu,lty=1,type='l',
         xlab=expression(mu),ylab=expression(paste("p(",mu,"|",sigma,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",mu,";          ",sigma^{(j)}," = ",.(round(sigv[j],digits=4)))), line=4) 
    title(bquote(paste("E(",mu,"|",sigma^{(j)},",y) = ybar = ",.(round(ybar,digits=4)),"       ",
                       "var(",mu,"|",sigma^{(j)},",y) = ",sigma^{2(j)},"/n = ",.(round(sigv[j]^2/nobs,digits=4)))), line=2)
  }
  
  ## Plot selected conditional posteriors for sigma 
  
  if(j==4){
    
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(sig,con_density_sigma,lty=1,type='l',
         xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",mu,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",sigma)), line=2) 
    
  }
  
  if(j==5){
    
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(sig,con_density_sigma,lty=1,type='l',
         xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",mu,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",sigma)), line=2) 
    
  }
  
  if(j==6){
    
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(sig,con_density_sigma,lty=1,type='l',
         xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",mu,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",sigma)), line=2) 
    
  }
  
  if(j==7){
    
    layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
    par(mar=c(5,5,5,1))
    plot(sig,con_density_sigma,lty=1,type='l',
         xlab=expression(sigma),ylab=expression(paste("p(",sigma,"|",mu,",y)"))
    )
    title(bquote(paste("Conditional posterior for ",sigma)), line=2) 
    
  }
  
  
}



############################################################
###
###  Estimation of the marginal posterior pdf of mu 
###  as an average of the stored conditionals
###
############################################################

## only repped through 2:repl not 1:repl therefore (repl-1)
ave_con_density_mu = cmu/(repl-1)

## ave_con_density_sigma = csig/repl


## Remove the starting value of mu for demonstration of histogram-based estimate


muv_trunc = muv[2:nrow(muv)]


## Histogram of draws of mu

hist(muv_trunc,30,ylab="Percentage frequencies",xlab=expression(mu),
     main=expression(paste("Histogram of draws of ",mu)))

## Kernel density plot of the draws of mu

d=density(muv_trunc)
plot(d,main=expression(paste("Kernel density estimate of marginal of ",mu)),
     xlab=expression(mu),ylab=expression(paste("p(",mu,"|y)")),xlim=c(lmu,upmu))


############################################################################################
###
###  Plot of the selected conditional pdfs of mu against the average of all repl conditionals
###
############################################################################################



layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
par(mar=c(5,5,2,1))
plot(mu,ave_con_density_mu,lty=2,type='l',
     main=expression(paste("Conditional posteriors for ",mu," and average of all 100 conditionals")),
     xlab=expression(mu),ylab=expression(paste("p(",mu,"|",sigma,",y)")),
     ylim=c(min(condmu4,condmu5,condmu6,condmu7,ave_con_density_mu),
            max(condmu4,condmu5,condmu6,condmu7,ave_con_density_mu)))
lines(mu,condmu4,type='l',lty=3)
lines(mu,condmu5,type='l',lty=4)
lines(mu,condmu6,type='l',lty=5)
lines(mu,condmu7,type='l',lty=6)
legend(2.6, max(condmu4,condmu5,condmu6,condmu7,ave_con_density_mu), c("Conditional 1","Conditional 2","Conditional 3","Conditional 4","Average of conditionals"),
       lty = c(3,4,5,6,2),  merge = TRUE,cex=0.8)



################################################################################################################
###
###  Plot of the estimated marginal pdf (average of all repl conditionals) of mu against the exact marginal pdf
###
################################################################################################################

## mu

layout(mat = matrix(c(1,1),1,1,byrow=TRUE), height = c(1), TRUE)
par(mar=c(5,5,2,1))
plot(mu,ave_con_density_mu,lty=2,type='l',main=expression(paste("Average of all 100 conditionals & exact marginal posterior for ",mu)),xlab=expression(mu),ylab=expression(paste("p(",mu,"|y)")),
     ylim=c(min(ave_con_density_mu,exactmu),max(ave_con_density_mu,exactmu)))
lines(mu,exactmu,type='l',lty=1)
legend(2.6, max(ave_con_density_mu,exactmu), c("Average of all conditionals", "Exact marginal"),
       lty = c(2, 1),  merge = TRUE,cex=0.8)




