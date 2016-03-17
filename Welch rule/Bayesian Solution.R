##### Bayesian solution to Welch ######

bayes.sep.t <- function(n1, n2, var1, var2) {
  s1 <- sqrt(var1)
  s2 <- sqrt(var2)
  v1 <- n1-1
  v2 <- n2-1
  cos.phi.sq.num <- (s2^2/n2)
  cos.phi.sq.denom <- s1^2/n1 + s2^2/n2
  cos.phi.sq <- cos.phi.sq.num/cos.phi.sq.denom
  sin.phi.sq <- 1-cos.phi.sq
  f1 <- (v2/(v2-2))*cos.phi.sq + (v1/(v1-2))*sin.phi.sq
  f2 <- (v2^2/((v2-2)^2*(v2-4)))*cos.phi.sq^2 + (v1^2/((v1-2)^2*(v1-4)))*sin.phi.sq^2
  b <- 4 + f1^2/f2
  a.sq <- ((b-2)/b)*f1
  result <- list(a.sq=a.sq, b=b, f1=f1, f2=f2)
  return(result)
}
 
bayes.sep.t(n1=20, n2=12, var1=12, var2=40)

df.student <- function(n1, n2, var1, var2) {
  n1+n2-2
}

df.welch <- function(n1, n2, var1, var2) {
  welch.num <- (var1/n1 + var2/n2)^2
  welch.denom <- var1^2/(n1^2*(n1-1)) + var2^2/(n2^2*(n2-1))
  welch.num/welch.denom
}

df.bayes <- function(n1, n2, var1, var2) {
  s1 <- sqrt(var1)
  s2 <- sqrt(var2)
  v1 <- n1-1
  v2 <- n2-1
  cos.phi.sq.num <- (s2^2/n2)
  cos.phi.sq.denom <- s1^2/n1 + s2^2/n2
  cos.phi.sq <- cos.phi.sq.num/cos.phi.sq.denom
  sin.phi.sq <- 1-cos.phi.sq
  f1 <- (v2/(v2-2))*cos.phi.sq + (v1/(v1-2))*sin.phi.sq
  f2 <- (v2^2/((v2-2)^2*(v2-4)))*cos.phi.sq^2 + (v1^2/((v1-2)^2*(v1-4)))*sin.phi.sq^2
  b <- 4 + f1^2/f2
  b
}

df.defaults <- list('n1'=50,'n2'=50,'var1'=2,'var2'=2)

partial.df.student <- function(var = 'a', params=df.ratio.defaults){
  params[[var]]=as.name('x')
  function(x)do.call(df.student, params)
}
partial.df.welch <- function(var = 'a', params=df.ratio.defaults){
  params[[var]]=as.name('x')
  function(x)do.call(df.welch, params)
}
partial.df.bayes <- function(var = 'a', params=df.ratio.defaults){
  params[[var]]=as.name('x')
  function(x)do.call(df.bayes, params)
}

# same Ns, different variances
ggplot(data.frame(x=seq(2,10,.1)), aes(x)) +
  #stat_function(fun=partial.df.student.vars(var='var2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2))) + 
  stat_function(fun=partial.df.welch.vars(var='var2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='red') + 
  stat_function(fun=partial.df.bayes.vars(var='var2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='blue') + 
  labs(y='degrees of freedom', x=bquote(frac(sigma[1]^2,sigma[2]^2))) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c(2,4,6,8,10)/2)

# different Ns, same variances
ggplot(data.frame(x=seq(20,100,.1)), aes(x)) + 
  stat_function(fun=partial.df.welch(var='n2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='red') +
  stat_function(fun=partial.df.bayes(var='n2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='blue') +
  labs(y='degrees of freedom', x=bquote(frac('n'[1],'n'[2]))) +
  scale_x_continuous(breaks=c(20,30,40,50,60,70,80,90,100), labels=c(20,30,40,50,60,70,80,90,100)/50)

# different NS, different vars
ggplot(data.frame(x=seq(2,10,.1)), aes(x)) + 
  stat_function(fun=partial.df.welch(var='var2', params=list('n1'=75,'n2'=50,'var1'=2,'var2'=2)), aes(colour='welch', linetype='n1=75, n2=50')) +
  stat_function(fun=partial.df.welch(var='var2', params=list('n1'=50,'n2'=75,'var1'=2,'var2'=2)), aes(colour='welch', linetype='n1=50, n2=75')) +
  stat_function(fun=partial.df.bayes(var='var2', params=list('n1'=75,'n2'=50,'var1'=2,'var2'=2)), aes(colour='bayes', linetype='n1=75, n2=50')) +
  stat_function(fun=partial.df.bayes(var='var2', params=list('n1'=50,'n2'=75,'var1'=2,'var2'=2)), aes(colour='bayes', linetype='n1=50, n2=75')) +
  labs(y='degrees of freedom', x=bquote(frac(sigma[1]^2,sigma[2]^2))) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c(2,4,6,8,10)/2)


se.welch <- function(n1, n2, var1, var2) {
  sqrt(var1/n1 + var2/n2)
}

se.bayes <- function(n1, n2, var1, var2) {
  v1 <- n1-1
  v2 <- n2-1
  cos.phi.sq.num <- (var2/n2)
  cos.phi.sq.denom <- var1/n1 + var2/n2
  cos.phi.sq <- cos.phi.sq.num/cos.phi.sq.denom
  sin.phi.sq <- 1-cos.phi.sq
  f1 <- (v2/(v2-2))*cos.phi.sq + (v1/(v1-2))*sin.phi.sq
  f2 <- (v2^2/((v2-2)^2*(v2-4)))*cos.phi.sq^2 + (v1^2/((v1-2)^2*(v1-4)))*sin.phi.sq^2
  b <- 4 + f1^2/f2
  a.sq <- ((b-2)/b)*f1
  sqrt(a.sq)*sqrt(var1/n1 + var2/n2)
}

partial.se.welch <- function(var = 'a', params=df.ratio.defaults){
  params[[var]]=as.name('x')
  function(x)do.call(se.welch, params)
}
partial.se.bayes <- function(var = 'a', params=df.ratio.defaults){
  params[[var]]=as.name('x')
  function(x)do.call(se.bayes, params)
}

# same Ns, different vars
ggplot(data.frame(x=seq(2,10,.1)), aes(x)) +
  stat_function(fun=partial.se.welch(var='var2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='red') + 
  stat_function(fun=partial.se.bayes(var='var2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='blue') + 
  labs(y='degrees of freedom', x=bquote(frac(sigma[1]^2,sigma[2]^2))) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c(2,4,6,8,10)/2)

# different Ns, same variances
ggplot(data.frame(x=seq(20,100,.1)), aes(x)) + 
  stat_function(fun=partial.se.welch(var='n2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='red') +
  stat_function(fun=partial.se.bayes(var='n2', params=list('n1'=50,'n2'=50,'var1'=2,'var2'=2)), color='blue') +
  labs(y='degrees of freedom', x=bquote(frac('n'[1],'n'[2]))) +
  scale_x_continuous(breaks=c(20,30,40,50,60,70,80,90,100), labels=c(20,30,40,50,60,70,80,90,100)/50)

# different NS, different vars
ggplot(data.frame(x=seq(2,10,.1)), aes(x)) + 
  stat_function(fun=partial.se.welch(var='var2', params=list('n1'=75,'n2'=50,'var1'=2,'var2'=2)), aes(colour='welch', linetype='n1=75, n2=50')) +
  stat_function(fun=partial.se.welch(var='var2', params=list('n1'=50,'n2'=75,'var1'=2,'var2'=2)), aes(colour='welch', linetype='n1=50, n2=75')) +
  stat_function(fun=partial.se.bayes(var='var2', params=list('n1'=75,'n2'=50,'var1'=2,'var2'=2)), aes(colour='bayes', linetype='n1=75, n2=50')) +
  stat_function(fun=partial.se.bayes(var='var2', params=list('n1'=50,'n2'=75,'var1'=2,'var2'=2)), aes(colour='bayes', linetype='n1=50, n2=75')) +
  labs(y='degrees of freedom', x=bquote(frac(sigma[1]^2,sigma[2]^2))) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c(2,4,6,8,10)/2)


a.sq <- function(n1, n2, var1, var2) {
  v1 <- n1-1
  v2 <- n2-1
  cos.phi.sq.num <- (var2/n2)
  cos.phi.sq.denom <- var1/n1 + var2/n2
  cos.phi.sq <- cos.phi.sq.num/cos.phi.sq.denom
  sin.phi.sq <- 1-cos.phi.sq
  f1 <- (v2/(v2-2))*cos.phi.sq + (v1/(v1-2))*sin.phi.sq
  f2 <- (v2^2/((v2-2)^2*(v2-4)))*cos.phi.sq^2 + (v1^2/((v1-2)^2*(v1-4)))*sin.phi.sq^2
  b <- 4 + f1^2/f2
  a.sq <- ((b-2)/b)*f1
  a.sq
}

a.sq.defaults <- list('n1'=50,'n2'=50,'var1'=2,'var2'=2)

partial.a.sq <- function(var = 'a', params=a.sq.defaults){
  params[[var]]=as.name('x')
  function(x)do.call(a.sq, params)
}

ggplot(data.frame(x=seq(2,100,.1)), aes(x)) + 
  stat_function(fun=partial.a.sq(var='n2', params=list('n1'=10,'n2'=10,'var1'=2,'var2'=2))) 


##### Rewrite function to compute Bayes value as well #####
t.contrast.bayes <- function(dv, groups, contrast) {
  means <- by(dv, groups, mean)
  vars <- by(dv, groups, var)
  Ns <- by(dv, groups, length)
  ihat <- contrast %*% means
  
  #classic
  df.classic <- sum(Ns)-length(Ns)
  mse <- sum(vars*(Ns-1))/df.classic
  se.classic <- sqrt(mse*(contrast^2 %*% (1/Ns)))
  t.classic <- ihat/se.classic
  p.classic <- 2*(1-pt(abs(t.classic), df.classic))
  
  #welch
  df.welch <- (contrast^2 %*% (vars/Ns))^2/(contrast^2 %*% (vars^2/(Ns^2*(Ns-1))))
  se.welch <- sqrt(contrast^2 %*% (vars/Ns))
  t.welch <- ihat/se.welch
  p.welch <- 2*(1-pt(abs(t.welch), df.welch))
  
  #bayes
  s1 <- sqrt(vars[1])
  s2 <- sqrt(vars[2])
  n1 <- Ns[1]
  n2 <- Ns[2]
  v1 <- n1-1
  v2 <- n2-1
  cos.phi.sq.num <- (s2^2/n2)
  cos.phi.sq.denom <- s1^2/n1 + s2^2/n2
  cos.phi.sq <- cos.phi.sq.num/cos.phi.sq.denom
  sin.phi.sq <- 1-cos.phi.sq
  f1 <- (v2/(v2-2))*cos.phi.sq + (v1/(v1-2))*sin.phi.sq
  f2 <- (v2^2/((v2-2)^2*(v2-4)))*cos.phi.sq^2 + (v1^2/((v1-2)^2*(v1-4)))*sin.phi.sq^2
  b <- 4 + f1^2/f2
  df.bayes <- b
  se.bayes <- sqrt(contrast^2 %*% (vars/Ns))
  t.bayes <- ihat/se.bayes
  p.bayes <- 2*(1-pt(abs(t.bayes), df.bayes))
  
  result <- list(ihat=ihat, est.vars=vars, se.classic=se.classic, t.classic=t.classic, df.classic=df.classic, p.classic=p.classic, se.welch=se.welch, t.welch=t.welch, df.welch=df.welch, p.welch=p.welch, se.bayes=se.bayes, t.bayes=t.bayes, df.bayes=df.bayes, p.bayes=p.bayes)
  return(result)
}

t.compare.bayes <- function(nsims, Ns, means, vars, contrast) {
  sims <- vector('list',nsims)
  group <- rep(1:length(Ns), Ns)   #vector of length N with group codes
  #dv <- vector('numeric',sum(Ns))
  
  sim.results <- lapply(sims, function(x){
    dv <- rnorm(n=sum(Ns), mean=rep(means,times=Ns), sd=sqrt(rep(vars,times=Ns)))
    
    sim.data <- data.frame(group, dv)
    
    fit <- t.contrast.bayes(dv,group,contrast)
    ihat <- fit$ihat
    est.vars <- fit$est.vars
    
    #save classic
    se.classic <- fit$se.classic
    t.classic <- fit$t.classic
    df.classic <- fit$df.classic
    p.classic <- fit$p.classic
    
    #save welch
    se.welch <- fit$se.welch
    t.welch <- fit$t.welch
    df.welch <- fit$df.welch
    p.welch <- fit$p.welch
    
    #save bayes
    se.bayes <- fit$se.bayes
    t.bayes <- fit$t.bayes
    df.bayes <- fit$df.bayes
    p.bayes <- fit$p.bayes
    
    current.sim <- list(sim.data, contrast, ihat, est.vars, se.classic, t.classic, df.classic, p.classic, se.welch, t.welch, df.welch, p.welch, se.bayes, t.bayes, df.bayes, p.bayes) # add matrix(c(dv,group), ncol=2, dimnames=list(c(),c('dv','group'))) to save the data
    names(current.sim) <- c('sim.data', 'contrast', 'ihat', 'est.vars','se.classic', 't.classic', 'df.classic', 'p.classic', 'se.welch', 't.welch', 'df.welch', 'p.welch', 'se.bayes', 't.bayes', 'df.bayes', 'p.bayes') # add data if saving the data
    return(current.sim)
  })
  
  names(sim.results)[1:length(sim.results)] <- paste('sim',1:length(sim.results),sep='')
  
  # store proportion of rejected null hypotheses
  classic.reject <- sum(lapply(sim.results, '[[', 'p.classic')<=.05)/nsims
  welch.reject <- sum(lapply(sim.results, '[[', 'p.welch')<=.05)/nsims
  bayes.reject <- sum(lapply(sim.results, '[[', 'p.bayes')<=.05)/nsims
  
  #store df ratio
  df.classic.vector <- unlist(lapply(sim.results, '[[', 'df.classic'))
  df.welch.vector <- unlist(lapply(sim.results, '[[', 'df.welch'))
  df.ratio <- df.welch.vector/df.classic.vector
  df.ratio.avg <- mean(df.ratio)
  
  # compute coverage rate
  true.ihat <- contrast %*% means
  
  obs.ihat <- data.frame(lapply(sim.results, '[[', 'ihat'))
  classic.df <- data.frame(lapply(sim.results, '[[', 'df.classic'))
  welch.df <- data.frame(lapply(sim.results, '[[', 'df.welch'))
  bayes.df <- data.frame(lapply(sim.results, '[[', 'df.bayes'))
  classic.ses <- data.frame(lapply(sim.results, '[[', 'se.classic'))
  welch.ses <- data.frame(lapply(sim.results, '[[', 'se.welch'))
  bayes.ses <- data.frame(lapply(sim.results, '[[', 'se.bayes'))
  
  t.classic <- apply(classic.df, 1, function(x){qt(.025, df=x)})        
  classic.lb <- obs.ihat-t.classic*classic.ses
  classic.ub <- obs.ihat+t.classic*classic.ses
  classic.coverage.logical <- (classic.ub-true.ihat)*(true.ihat-classic.lb)>0
  
  t.welch <- apply(welch.df, 1, function(x){qt(.025, df=x)})        
  welch.lb <- obs.ihat-t.welch*welch.ses
  welch.ub <- obs.ihat+t.welch*welch.ses
  welch.coverage.logical <- (welch.ub-true.ihat)*(true.ihat-welch.lb)>0   
  
  t.bayes <- apply(bayes.df, 1, function(x){qt(.025, df=x)})
  bayes.lb <- obs.ihat-t.bayes*bayes.ses
  bayes.ub <- obs.ihat+t.bayes*bayes.ses
  bayes.coverage.logical <- (bayes.ub-true.ihat)*(true.ihat-bayes.lb)>0
  
  classic.coverage <- sum(classic.coverage.logical)/nsims
  welch.coverage <- sum(welch.coverage.logical)/nsims
  bayes.coverage <- sum(bayes.coverage.logical)/nsims
  
  # return data 
  return(list(classic.reject=classic.reject, welch.reject=welch.reject, bayes.reject=bayes.reject, df.ratio.avg=df.ratio.avg, classic.coverage=classic.coverage, welch.coverage=welch.coverage, bayes.coverage=bayes.coverage, df.ratio=df.ratio, sim.results=sim.results))
}

set.seed(2184)
ve.ns20.ne.nullT <- t.compare.bayes(nsims=1000, Ns=c(20,20), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
set.seed(2184)
biggroupsmallvar <- t.compare.bayes(nsims=1000, Ns=c(30,10), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
set.seed(2184)
biggroupbigvar <- t.compare.bayes(nsims=1000, Ns=c(10,30), means=c(6,6), vars=c(10,2), contrast=c(-1,1))
set.seed(2184)
biggroupsmallvar <- t.compare.bayes(nsims=1000, Ns=c(30,10), means=c(6,8), vars=c(2,10), contrast=c(-1,1))
set.seed(2184)
biggroupbigvar <- t.compare.bayes(nsims=1000, Ns=c(10,30), means=c(6,8), vars=c(10,2), contrast=c(-1,1))
