# Fit double observer model to matched data
library(R2jags)
#plot function
#plots detection rates from a JAGS model fit object, specific to this context and huggins5
plotdobs<-function(fit, position=1, group=1, nspp=18, pt.legend=c("FRONT", "BACK"), 
                   jit=0.1, spplab=NULL){
  out = fit$BUGSoutput$sims.list
  mu <- matrix(out$beta0[,group],nrow=dim(out[[1]])[1], ncol=nspp) + 
    out$alpha[,1:nspp]
  post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
  plot(1,1,type="n", xlim=c(0, nspp+1), ylim=c(0,1), xaxt="n", xlab="", ylab="")
  axis(side=1, las=2, at=1:nspp, labels=spplab)
  grid()
  arrows(x0=1:nspp, x1=1:nspp, y0=post[1,], y1=post[3,], length=0, lwd=2)
  points(1:nspp, post[2,], pch=21, bg=1)
  mu <- matrix(out$beta0[,group]+(out$beta1[,1]*position), nrow=dim(out[[1]])[1], ncol=nspp) + 
    out$alpha[,(nspp+1):(2*nspp)]
  post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
  jit <- 0.1
  arrows(x0=1:nspp+jit, x1=1:nspp+jit, y0=post[1,], y1=post[3,], length=0, lwd=2, col="lightgray")
  points(1:nspp+jit, post[2,], pch=21, bg="lightgray")
  legend("bottomleft", legend=pt.legend, lwd=2, col=c(1, "lightgray"))
}
#specify model
## Model 1: group + position + species x position interaction
cat(file="Model_1.txt", "
    model{
    for(i in 1:4){
    beta0[i] ~ dnorm(0, 0.01) #group size, single is intercept
    }
    beta1 ~ dnorm(0, 0.01) #position effect
    
    for( i in 1:36){
    alpha[i] ~dnorm(0, tau[1]) #species (i = 1 - 18 front) and position (i = 19 - 36 rear) 
    }                          #interaction effect
    sd[1] ~ dunif(0, 10)
    tau[1] <- pow(sd[1], -2)
    
    for(s in 1:18){
    for(g in 1:4){
    logit(p1[g,s]) <- beta0[g] + alpha[s]
    logit(p2[g,s]) <- beta0[g] + beta1 + alpha[s+18]
    p[1,g,s] <- (1 - p1[g,s])*p2[g,s]  #01 capture history                           
    p[2,g,s] <- p1[g,s]*(1 - p2[g,s])  #10 capture history
    p[3,g,s] <- p1[g,s]*p2[g,s]  #11 capture history 
    p[4,g,s] <- 1 - p[1,g,s] - p[2,g,s] - p[3,g,s]  #00 capture history
    
    for(i in 1:3){
    mu[i,g,s] <- p[i,g,s]/(1 - p[4,g,s])
    }
    y[1:3,g,s] ~ dmulti( mu[1:3,g,s], M[g,s] )
    }}
    }")

load(file = "Make Matched Data/.RData")
#why is there a "WWLNA" crew and a "HMWNA" crew?
match <- match[match$crew %in% c("HMWWWL", "HMWTKZ", "WWLTKZ"),]
jags.data = list(y=as.array(table(match$ch, match$grpcode, match$sppn)), 
                 M=as.matrix(table(match$grpcode, match$sppn)) )
inits <- function(){list( beta0 = runif(4, -0.1, 0.1), beta1 = runif(1, -0.1, 0.1), 
                          alpha=runif(36, -0.1, 0.1), sd = runif(1, 0.1, 0.2) )}
parms = c("beta0", "beta1", "alpha", "sd")
fit <- jags(jags.data, inits, parms, "Model_1.txt", n.chains=4, n.iter = 4000, n.burnin = 2000)
plot(fit)
summary(fit$BUGSoutput$summary[,"Rhat"])
plotdobs(fit=fit, position=1, spplab=sort(unique(match$sppn)))

#Model 2: group + species + position:species + position:crew
cat(file="Model_2.txt", "
    model{
    for(i in 1:4){
    beta0[i] ~ dnorm(0, 1) #group size
    }
    
    for( i in 1:18){
    alpha[i] ~dnorm(0, tau[1])  #spp effect
    }
    sd[1] ~ dunif(0, 1)
    tau[1] <- pow(sd[1], -2)
    for( i in 1:36){
    gamma[i] ~ dnorm(0, tau[2])  #position:spp effect
    }
    sd[2] ~ dunif(0, 1)
    tau[2] <- pow(sd[2], -2)
    
    for(i in 1:(2*N.crew)){
    delta[i] ~ dnorm(0, tau[3])  #position:crew effect
    }
    sd[3] ~ dunif(0, 1)
    tau[3] <- pow(sd[3], -2)
    
    
    for(s in 1:18){
    for(c in 1:N.crew){
    logit(p1[1,s,c]) <- beta0[1] + alpha[s] + gamma[s] + delta[c]
    logit(p2[1,s,c]) <- beta0[1] + alpha[s] + gamma[s+18] + delta[c+N.crew]
    for(g in 2:4){
    logit(p1[g,s,c]) <- beta0[1] + beta0[g] + alpha[s] + gamma[s] + delta[c]
    logit(p2[g,s,c]) <- beta0[1] + beta0[g] + alpha[s] + gamma[s+18] + delta[c+N.crew]
    }}}
    for(g in 1:4){
    for(s in 1:18){
    for(c in 1:N.crew){
    p[1,g,s,c] <- (1 - p1[g,s,c])*p2[g,s,c]  #01 capture history                           
    p[2,g,s,c] <- p1[g,s,c]*(1 - p2[g,s,c])  #10 capture history
    p[3,g,s,c] <- p1[g,s,c]*p2[g,s,c]  #11 capture history 
    p[4,g,s,c] <- 1 - p[1,g,s,c] - p[2,g,s,c] - p[3,g,s,c]  #00 capture history
    
    for(i in 1:3){
    mu[i,g,s,c] <- p[i,g,s,c]/(1 - p[4,g,s,c])
    }
    y[1:3,g,s,c] ~ dmulti( mu[1:3,g,s,c], M[g,s,c] )
    }}}
    }")
jags.data = list(y=as.array(table(match$ch, match$grpcode, match$sppn, match$crew)), 
                 M=as.array(table(match$grpcode, match$sppn, match$crew)), 
                 N.crew=3)
inits <- function(){list( beta0 = runif(4, -0.1, 0.1),  
                          alpha=runif(18, -0.1, 0.1), gamma = runif(36, -0.1, 0.1), 
                          delta = runif(6, -0.01, 0.01), sd = runif(3, 0.1, 0.2) )}
parms = c("beta0", "alpha", "gamma", "delta", "sd")
fit2 <- jags(jags.data, inits, parms, "Model_2.txt", n.chains=4, n.iter = 4000, n.burnin = 2000)
plot(fit2)
summary(fit2$BUGSoutput$summary[,"Rhat"])
#plot function for model 2
plotdobs2<-function(fit, position=1, crew=0, nspp=18, pt.legend=c("FRONT", "BACK"), jit=0.1, 
                    spplab=NULL){
  out = fit$BUGSoutput$sims.list
  if(crew == 0) {
    mu <- matrix(out$beta0[,1],nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
    out$gamma[,1:nspp] 
  } else {
    mu <- matrix(out$beta0[,1],nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,1:nspp] + matrix(out$delta[,crew],nrow=dim(out[[1]])[1], ncol=nspp)
  }
  post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
  plot(1,1,type="n", xlim=c(0, nspp+1), ylim=c(0,1), xaxt="n", xlab="", ylab="")
  axis(side=1, las=2, at=1:nspp, labels=spplab)
  grid()
  arrows(x0=1:nspp, x1=1:nspp, y0=post[1,], y1=post[3,], length=0, lwd=2)
  points(1:nspp, post[2,], pch=21, bg=1)
  if(crew == 0) {
    mu <- matrix(out$beta0[,1], nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,(nspp+1):(2*nspp)] 
  } else{
    mu <- matrix(out$beta0[,1], nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,(nspp+1):(2*nspp)] + matrix(out$delta[,crew],nrow=dim(out[[1]])[1], ncol=nspp)
  }
  post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
  jit <- 0.1
  arrows(x0=1:nspp+jit, x1=1:nspp+jit, y0=post[1,], y1=post[3,], length=0, lwd=2, col="lightgray")
  points(1:nspp+jit, post[2,], pch=21, bg="lightgray")
  legend("bottomleft", legend=pt.legend, lwd=2, col=c(1, "lightgray"))
}
plotdobs2(fit=fit2, position=1, crew=1, spplab=sort(unique(match$sppn)),  pt.legend=c("HMW", "TKZ"))
plotdobs2(fit=fit2, position=1, crew=2, spplab=sort(unique(match$sppn)),  pt.legend=c("HMW", "WWL"))
plotdobs2(fit=fit2, position=1, crew=3, spplab=sort(unique(match$sppn)),  pt.legend=c("WWL", "TKZ"))
plotdobs2(fit=fit2, position=1, crew=0, spplab=sort(unique(match$sppn)))

#Model 3:  group + species + species:position + species:position:observer
cat(file="Model_3.txt", "
    model{
    for(i in 1:4){
    beta0[i] ~ dnorm(0, 1) #group size
    }
    
    for( i in 1:18){
    alpha[i] ~dnorm(0, tau[1])  #spp effect
    }
    sd[1] ~ dunif(0, 1)
    tau[1] <- pow(sd[1], -2)
    for(i in 1:2){for( j in 1:18){
    gamma[i,j] ~ dnorm(0, tau[2])  #position:spp effect
    }}
    sd[2] ~ dunif(0, 1)
    tau[2] <- pow(sd[2], -2)
    
    for(i in 1:4){for(j in 1:18){
    delta[i,j] ~ dnorm(0, tau[3])  #position:spp:observer effect
    }}
    sd[3] ~ dunif(0, 1)
    tau[3] <- pow(sd[3], -2)
    
    
    for(s in 1:18){
    for(c in 1:3){
    logit(p1[1,s,c]) <- beta0[1] + alpha[s] + gamma[1, s ] + delta[ front[c],s ]
    logit(p2[1,s,c]) <- beta0[1] + alpha[s] + gamma[2, s ] + delta[ rear[c],s ]
    for(g in 2:4){
    logit(p1[g,s,c]) <- beta0[1] + beta0[g] + alpha[s] + gamma[1,s] + delta[ front[c],s ]
    logit(p2[g,s,c]) <- beta0[1] + beta0[g] + alpha[s] + gamma[2,s] + delta[ rear[c],s ]
    }}}
    for(g in 1:4){
    for(s in 1:18){
    for(c in 1:3){
    p[1,g,s,c] <- (1 - p1[g,s,c])*p2[g,s,c]  #01 capture history                           
    p[2,g,s,c] <- p1[g,s,c]*(1 - p2[g,s,c])  #10 capture history
    p[3,g,s,c] <- p1[g,s,c]*p2[g,s,c]  #11 capture history 
    p[4,g,s,c] <- 1 - p[1,g,s,c] - p[2,g,s,c] - p[3,g,s,c]  #00 capture history
    
    for(i in 1:3){
    mu[i,g,s,c] <- p[i,g,s,c]/(1 - p[4,g,s,c])
    }
    y[1:3,g,s,c] ~ dmulti( mu[1:3,g,s,c], M[g,s,c] )
    }}}
    }")
jags.data = list(y=as.array(table(match$ch, match$grpcode, match$sppn, match$crew)), 
                 M=as.array(table(match$grpcode, match$sppn, match$crew)), 
                 front=c(1,1,2), rear=c(4,3,4))
inits <- function(){list( beta0 = runif(4, -0.1, 0.1),  
                          alpha=runif(18, -0.1, 0.1), 
                          gamma = matrix(runif(36, -0.1, 0.1), 2, 18),  
                          delta = matrix(runif(4*18, -0.01, 0.01), 4, 18), 
                          sd = runif(3, 0.1, 0.2) )}
parms = c("beta0", "alpha", "gamma", "delta", "sd")
fit3 <- jags(jags.data, inits, parms, "Model_3.txt", n.chains=4, n.iter = 4000, n.burnin = 2000)
plot(fit3)
summary(fit3$BUGSoutput$summary[,"Rhat"])
plotdobs3<-function(fit, position=1, crew=0, nspp=18, pt.legend=c("FRONT", "BACK"), jit=0.1, 
                    spplab=NULL, main=NULL, front=c(1,1,2), rear=c(4,3,4)){
  out = fit$BUGSoutput$sims.list
  if(crew == 0) {
    mu <- matrix(out$beta0[,1],nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,1,1:nspp] 
  } else {
    mu <- matrix(out$beta0[,1],nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,1,1:nspp] + matrix(out$delta[,front[crew], 1:nspp],nrow=dim(out[[1]])[1], ncol=nspp)
  }
  post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
  plot(1,1,type="n", xlim=c(0, nspp+1), ylim=c(0,1), xaxt="n", xlab="", ylab="", main=main)
  axis(side=1, las=2, at=1:nspp, labels=spplab)
  grid()
  arrows(x0=1:nspp, x1=1:nspp, y0=post[1,], y1=post[3,], length=0, lwd=2)
  points(1:nspp, post[2,], pch=21, bg=1)
  if(crew == 0) {
    mu <- matrix(out$beta0[,1], nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,2,1:nspp] 
  } else{
    mu <- matrix(out$beta0[,1], nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,2,1:nspp] + matrix(out$delta[,rear[crew], 1:nspp],nrow=dim(out[[1]])[1], ncol=nspp)
  }
  post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
  jit <- 0.1
  arrows(x0=1:nspp+jit, x1=1:nspp+jit, y0=post[1,], y1=post[3,], length=0, lwd=2, col="lightgray")
  points(1:nspp+jit, post[2,], pch=21, bg="lightgray")
  legend("bottomleft", legend=pt.legend, lwd=2, col=c(1, "lightgray"))
}

plotdobs3(fit=fit3, crew=1, spplab=sort(unique(match$sppn)), pt.legend=c("HMW", "TKZ"))
plotdobs3(fit=fit3, crew=2, spplab=sort(unique(match$sppn)), pt.legend=c("HMW", "WWL"))
plotdobs3(fit=fit3, crew=3, spplab=sort(unique(match$sppn)), pt.legend=c("WWL", "TKZ"))
plotdobs3(fit=fit3, crew=0, spplab=sort(unique(match$sppn)), pt.legend=c("FRONT", "BACK"))

#Might fit model without species:position effect, gamma above

#Model_4:  similar to Model_3 but with year:position:species:observer
cat(file="Model_4.txt", "
    model{
    for(i in 1:4){
    beta0[i] ~ dnorm(0, 1) #group size
    }
    
    for( i in 1:18){
    alpha[i] ~dnorm(0, tau[1])  #spp effect
    }

    for(i in 1:2){for( j in 1:18){
    gamma[i,j] ~ dnorm(0, tau[2])  #position:spp effect
    }}

    for(k in 1:2){for(i in 1:4){for(j in 1:18){
    delta[k,i,j] ~ dnorm(0, tau[3])  #year:position:spp:observer effect
    }}}
    
    for( i in 1:3){
    sd[i] ~ dgamma(2, 15)
    tau[i] <- pow(sd[i], -2)
    }
    
    for(s in 1:18){
    for(c in 1:3){
    for(k in 1:2){
    logit(p1[1,s,c,k]) <- beta0[1] + alpha[s] + gamma[1, s ] + 
                            delta[ k,front[c],s ]
    logit(p2[1,s,c,k]) <- beta0[1] + alpha[s] + gamma[2, s ] + 
                            delta[ k,rear[c],s ]
    for(g in 2:4){
    logit(p1[g,s,c,k]) <- beta0[1] + beta0[g] + alpha[s] + gamma[1,s] + 
                            delta[ k,front[c],s ]
    logit(p2[g,s,c,k]) <- beta0[1] + beta0[g] + alpha[s] + gamma[2,s] + 
                            delta[ k,rear[c],s ]
    }}}}

    for(g in 1:4){
    for(s in 1:18){
    for(c in 1:3){
    for(k in 1:2){
    p[1,g,s,c,k] <- (1 - p1[g,s,c,k])*p2[g,s,c,k]  #01 capture history                           
    p[2,g,s,c,k] <- p1[g,s,c,k]*(1 - p2[g,s,c,k])  #10 capture history
    p[3,g,s,c,k] <- p1[g,s,c,k]*p2[g,s,c,k]  #11 capture history 
    p[4,g,s,c,k] <- 1 - p[1,g,s,c,k] - p[2,g,s,c,k] - p[3,g,s,c,k]  #00 capture history
    
    for(i in 1:3){
    mu[i,g,s,c,k] <- p[i,g,s,c,k]/(1 - p[4,g,s,c,k])
    }
    y[1:3,g,s,c,k] ~ dmulti( mu[1:3,g,s,c,k], M[g,s,c,k] )
    }}}}
    }")
jags.data = list(y=as.array(table(match$ch, match$grpcode, match$sppn, match$crew, match$yr)), 
                 M=as.array(table(match$grpcode, match$sppn, match$crew, match$yr)), 
                 front=c(1,1,2), rear=c(4,3,4))
inits <- function(){list( beta0 = runif(4, -0.1, 0.1),  
                          alpha=runif(18, -0.1, 0.1), 
                          gamma = matrix(runif(36, -0.1, 0.1), 2, 18),  
                          delta = array(runif(2*4*18, -0.01, 0.01), dim = c(2,4,18)),
                          sd = runif(3, 0.1, 0.2) )}
parms = c("beta0", "alpha", "gamma", "delta", "sd")
fit4 <- jags(jags.data, inits, parms, "Model_4.txt", n.chains=4, n.iter = 14000, n.burnin = 10000)
plot(fit4)
summary(fit4$BUGSoutput$summary[,"Rhat"])
plotdobs4<-function(fit, position=1, crew=0, nspp=18, pt.legend=c("FRONT", "BACK"), jit=0.1, 
                    spplab=NULL, main=NULL, front=c(1,1,2), rear=c(4,3,4)){
  out = fit$BUGSoutput$sims.list
  if(crew == 0) {
    mu <- matrix(out$beta0[,1],nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,1,1:nspp] 
    post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
    plot(1,1,type="n", xlim=c(0, nspp+1), ylim=c(0,1), xaxt="n", xlab="", ylab="", main=main)
    axis(side=1, las=2, at=1:nspp, labels=spplab)
    grid()
    arrows(x0=1:nspp, x1=1:nspp, y0=post[1,], y1=post[3,], length=0, lwd=2)
    points(1:nspp, post[2,], pch=21, bg=1)
    mu <- matrix(out$beta0[,1], nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,2,1:nspp] 
    post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
    jit <- 0.1
    arrows(x0=1:nspp+jit, x1=1:nspp+jit, y0=post[1,], y1=post[3,], length=0, lwd=2, col="lightgray")
    points(1:nspp+jit, post[2,], pch=21, bg="lightgray")
    legend("bottomleft", legend=pt.legend, lwd=2, col=c(1, "lightgray"))
    
  } else {
    mu <- matrix(out$beta0[,1],nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,1,1:nspp] + matrix(out$delta[,1,front[crew], 1:nspp],nrow=dim(out[[1]])[1], ncol=nspp)
    post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
    plot(1,1,type="n", xlim=c(0, nspp+1), ylim=c(0,1), xaxt="n", xlab="", ylab="", main=main)
    axis(side=1, las=2, at=1:nspp, labels=spplab)
    grid()
    arrows(x0=1:nspp, x1=1:nspp, y0=post[1,], y1=post[3,], length=0, lwd=2)
    points(1:nspp, post[2,], pch=21, bg=1)
    mu <- matrix(out$beta0[,1],nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,1,1:nspp] + matrix(out$delta[,2,front[crew], 1:nspp],nrow=dim(out[[1]])[1], ncol=nspp)
    post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
    arrows(x0=1:nspp-jit, x1=1:nspp-jit, y0=post[1,], y1=post[3,], length=0, lwd=2)
    points(1:nspp-jit, post[2,], pch=21, bg=1)
    mu <- matrix(out$beta0[,1], nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,2,1:nspp] + matrix(out$delta[,1,rear[crew], 1:nspp],nrow=dim(out[[1]])[1], ncol=nspp)
    post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
    arrows(x0=1:nspp+jit, x1=1:nspp+jit, y0=post[1,], y1=post[3,], length=0, lwd=2, col="lightgray")
    points(1:nspp+jit, post[2,], pch=21, bg="lightgray")
    mu <- matrix(out$beta0[,1], nrow=dim(out[[1]])[1], ncol=nspp) + out$alpha[,1:nspp] + 
      out$gamma[,2,1:nspp] + matrix(out$delta[,2,rear[crew], 1:nspp],nrow=dim(out[[1]])[1], ncol=nspp)
    post = apply(plogis(mu), 2, quantile, probs=c(0.025, 0.5, 0.975))
    arrows(x0=1:nspp+2*jit, x1=1:nspp+2*jit, y0=post[1,], y1=post[3,], length=0, lwd=2, col="lightgray")
    points(1:nspp+2*jit, post[2,], pch=21, bg="lightgray")
    legend("bottomleft", legend=pt.legend, lwd=2, col=c(1, "lightgray"))
  }
}

plotdobs4(fit=fit4, crew=1, spplab=sort(unique(match$sppn)), pt.legend=c("HMW", "TKZ"))
plotdobs4(fit=fit4, crew=2, spplab=sort(unique(match$sppn)), pt.legend=c("HMW", "WWL"))
plotdobs4(fit=fit4, crew=3, spplab=sort(unique(match$sppn)), pt.legend=c("WWL", "TKZ"))
plotdobs4(fit=fit4, crew=0, spplab=sort(unique(match$sppn)), pt.legend=c("FRONT", "BACK"))
