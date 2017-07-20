## function to compute classic and Welch contrasts
# contrast is an n x m matrix with each of the n rows as a contrast and each of the m columns representing a group
t.contrast <- function(dv, groups, contrast) {
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
    
    result <- list(ihat=ihat, est.vars=vars, se.classic=se.classic, t.classic=t.classic, df.classic=df.classic, p.classic=p.classic, se.welch=se.welch, t.welch=t.welch, df.welch=df.welch, p.welch=p.welch)
    return(result)
}

## function to run simulations
# contrast is an n x m matrix with each of the n rows as a contrast and each of the m columns representing a group
t.compare <- function(nsims, Ns, means, vars, contrast) {
    sims <- vector('list',nsims)
    group <- rep(1:length(Ns), Ns)   #vector of length N with group codes
    #dv <- vector('numeric',sum(Ns))
    
    sim.results <- lapply(sims, function(x){
        dv <- rnorm(n=sum(Ns), mean=rep(means,times=Ns), sd=sqrt(rep(vars,times=Ns)))

        sim.data <- data.frame(group, dv)
        
        fit <- t.contrast(dv,group,contrast)
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
        
        current.sim <- list(sim.data, contrast, ihat, est.vars, se.classic, t.classic, df.classic, p.classic, se.welch, t.welch, df.welch, p.welch) # add matrix(c(dv,group), ncol=2, dimnames=list(c(),c('dv','group'))) to save the data
        names(current.sim) <- c('sim.data', 'contrast', 'ihat', 'est.vars','se.classic', 't.classic', 'df.classic', 'p.classic', 'se.welch', 't.welch', 'df.welch', 'p.welch') # add data if saving the data
        return(current.sim)
    })
    
    names(sim.results)[1:length(sim.results)] <- paste('sim',1:length(sim.results),sep='')
    
    # store proportion of rejected null hypotheses
    classic.reject <- sum(lapply(sim.results, '[[', 'p.classic')<=.05)/nsims
    welch.reject <- sum(lapply(sim.results, '[[', 'p.welch')<=.05)/nsims
    
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
    classic.ses <- data.frame(lapply(sim.results, '[[', 'se.classic'))
    welch.ses <- data.frame(lapply(sim.results, '[[', 'se.welch'))
    
    t.classic <- apply(classic.df, 1, function(x){qt(.025, df=x)})        
    classic.lb <- obs.ihat-t.classic*classic.ses
    classic.ub <- obs.ihat+t.classic*classic.ses
    classic.coverage.logical <- (classic.ub-true.ihat)*(true.ihat-classic.lb)>0
    
    t.welch <- apply(welch.df, 1, function(x){qt(.025, df=x)})        
    welch.lb <- obs.ihat-t.welch*welch.ses
    welch.ub <- obs.ihat+t.welch*welch.ses
    welch.coverage.logical <- (welch.ub-true.ihat)*(true.ihat-welch.lb)>0       
    
    classic.coverage <- sum(classic.coverage.logical)/nsims
    welch.coverage <- sum(welch.coverage.logical)/nsims
    
    # return data 
    return(list(classic.reject=classic.reject, welch.reject=welch.reject, df.ratio.avg=df.ratio.avg, classic.coverage=classic.coverage, welch.coverage=welch.coverage, df.ratio=df.ratio, sim.results=sim.results))
}



##### SIMULATIONS #####

#Compute mean differences to obtain specific effect sizes
cohen.diff <- function(d, N1, SD1, N2, SD2) #insert the mean, sample size, and standard deviation for each group
{
    poolSD <- sqrt(((N1-1)*SD1^2+(N2-1)*SD2^2)/(N1+N2-2)) #computes the pooled standard deviation
    diff <- d*poolSD #computes difference needed
    print(diff) #displays difference needed
}


##### Equal Ns, equal variances simulations #####

# true null
set.seed(2184)
ve.ns20.ne.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns50.ne.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns100.ne.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,2), contrast=c(-1,1))

# small effect
cohen.diff(.2, 20, sqrt(2), 20, sqrt(2))
set.seed(2184)
ve.ns20.ne.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))
ve.ns50.ne.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))
ve.ns100.ne.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(2))
set.seed(2184)
ve.ns20.ne.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))
ve.ns50.ne.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))
ve.ns100.ne.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(2))
set.seed(2184)
ve.ns20.ne.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))
ve.ns50.ne.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))
ve.ns100.ne.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))


## Save simulations
save(ve.ns20.ne.nullT,ve.ns50.ne.nullT,ve.ns100.ne.nullT,ve.ns20.ne.smalld,ve.ns50.ne.smalld,ve.ns100.ne.smalld,ve.ns20.ne.midd,ve.ns50.ne.midd,ve.ns100.ne.midd,ve.ns20.ne.bigd,ve.ns50.ne.bigd,ve.ns100.ne.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/veNeSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/veNeSeed2184.Rdata')
venesim1 <- list(ve.ns20.ne.nullT=ve.ns20.ne.nullT$sim.results$sim1$sim.data,
                 ve.ns50.ne.nullT=ve.ns50.ne.nullT$sim.results$sim1$sim.data,
                 ve.ns100.ne.nullT=ve.ns100.ne.nullT$sim.results$sim1$sim.data,
                 ve.ns20.ne.smalld=ve.ns20.ne.smalld$sim.results$sim1$sim.data,
                 ve.ns50.ne.smalld=ve.ns50.ne.smalld$sim.results$sim1$sim.data,
                 ve.ns100.ne.smalld=ve.ns100.ne.smalld$sim.results$sim1$sim.data,
                 ve.ns20.ne.midd=ve.ns20.ne.midd$sim.results$sim1$sim.data,
                 ve.ns50.ne.midd=ve.ns50.ne.midd$sim.results$sim1$sim.data,
                 ve.ns100.ne.midd=ve.ns100.ne.midd$sim.results$sim1$sim.data,
                 ve.ns20.ne.bigd=ve.ns20.ne.bigd$sim.results$sim1$sim.data,
                 ve.ns50.ne.bigd=ve.ns50.ne.bigd$sim.results$sim1$sim.data,
                 ve.ns100.ne.bigd=ve.ns100.ne.bigd$sim.results$sim1$sim.data)
save(venesim1, file='/users/joshwondra/R-projects/Welch rule/veNeSeed2184-sim1.Rdata')


##### Equal Ns, equal variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.ve.ne <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.ve.ne[,1,1] <- c(ve.ns20.ne.nullT$classic.reject, ve.ns50.ne.nullT$classic.reject, ve.ns100.ne.nullT$classic.reject)
reject.null.ve.ne[,1,2] <- c(ve.ns20.ne.nullT$welch.reject, ve.ns50.ne.nullT$welch.reject, ve.ns100.ne.nullT$welch.reject)

# small effect
reject.null.ve.ne[,2,1] <- c(ve.ns20.ne.smalld$classic.reject, ve.ns50.ne.smalld$classic.reject, ve.ns100.ne.smalld$classic.reject)
reject.null.ve.ne[,2,2] <- c(ve.ns20.ne.smalld$welch.reject, ve.ns50.ne.smalld$welch.reject, ve.ns100.ne.smalld$welch.reject)

# medium effect
reject.null.ve.ne[,3,1] <- c(ve.ns20.ne.midd$classic.reject, ve.ns50.ne.midd$classic.reject, ve.ns100.ne.midd$classic.reject)
reject.null.ve.ne[,3,2] <- c(ve.ns20.ne.midd$welch.reject, ve.ns50.ne.midd$welch.reject, ve.ns100.ne.midd$welch.reject)

# large effect
reject.null.ve.ne[,4,1] <- c(ve.ns20.ne.bigd$classic.reject, ve.ns50.ne.bigd$classic.reject, ve.ns100.ne.bigd$classic.reject)
reject.null.ve.ne[,4,2] <- c(ve.ns20.ne.bigd$welch.reject, ve.ns50.ne.bigd$welch.reject, ve.ns100.ne.bigd$welch.reject)


##### Equal Ns, equal variances coverage rate #####

## store observed coverage rates
obs.coverage.ve.ne <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.ve.ne[,1,1] <- c(ve.ns20.ne.nullT$classic.coverage, ve.ns50.ne.nullT$classic.coverage, ve.ns100.ne.nullT$classic.coverage)
obs.coverage.ve.ne[,1,2] <- c(ve.ns20.ne.nullT$welch.coverage, ve.ns50.ne.nullT$welch.coverage, ve.ns100.ne.nullT$welch.coverage)

# small effect
obs.coverage.ve.ne[,2,1] <- c(ve.ns20.ne.smalld$classic.coverage, ve.ns50.ne.smalld$classic.coverage, ve.ns100.ne.smalld$classic.coverage)
obs.coverage.ve.ne[,2,2] <- c(ve.ns20.ne.smalld$welch.coverage, ve.ns50.ne.smalld$welch.coverage, ve.ns100.ne.smalld$welch.coverage)

# medium effect
obs.coverage.ve.ne[,3,1] <- c(ve.ns20.ne.midd$classic.coverage, ve.ns50.ne.midd$classic.coverage, ve.ns100.ne.midd$classic.coverage)
obs.coverage.ve.ne[,3,2] <- c(ve.ns20.ne.midd$welch.coverage, ve.ns50.ne.midd$welch.coverage, ve.ns100.ne.midd$welch.coverage)

# big effect
obs.coverage.ve.ne[,4,1] <- c(ve.ns20.ne.bigd$classic.coverage, ve.ns50.ne.bigd$classic.coverage, ve.ns100.ne.bigd$classic.coverage)
obs.coverage.ve.ne[,4,2] <- c(ve.ns20.ne.bigd$welch.coverage, ve.ns50.ne.bigd$welch.coverage, ve.ns100.ne.bigd$welch.coverage)



##### Equal Ns, equal variances df proportion #####


## store df proportion
df.ratio.ve.ne <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.ve.ne[,1,1] <- c(ve.ns20.ne.nullT$df.ratio.avg, ve.ns50.ne.nullT$df.ratio.avg, ve.ns100.ne.nullT$df.ratio.avg)  # true null
df.ratio.ve.ne[,2,1] <- c(ve.ns20.ne.smalld$df.ratio.avg, ve.ns50.ne.smalld$df.ratio.avg, ve.ns100.ne.smalld$df.ratio.avg)  # small d
df.ratio.ve.ne[,3,1] <- c(ve.ns20.ne.midd$df.ratio.avg, ve.ns50.ne.midd$df.ratio.avg, ve.ns100.ne.midd$df.ratio.avg)  # mid d
df.ratio.ve.ne[,4,1] <- c(ve.ns20.ne.bigd$df.ratio.avg, ve.ns50.ne.bigd$df.ratio.avg, ve.ns100.ne.bigd$df.ratio.avg)  # big d




##### Save tables from equal Ns, equal variances simulations #####
save(df.ratio.ve.ne, obs.coverage.ve.ne, reject.null.ve.ne, file='/users/joshwondra/R-projects/Welch rule/veNeSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/veNeSeed2184Tables.Rdata')




##### Equal Ns, double variance simulations #####

# true null
set.seed(2184)
v2.ns20.ne.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns50.ne.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns100.ne.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,4), contrast=c(-1,1))

# small effect
cohen.diff(.2, 20, sqrt(2), 20, sqrt(4))
set.seed(2184)
v2.ns20.ne.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))
v2.ns50.ne.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))
v2.ns100.ne.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(4))
set.seed(2184)
v2.ns20.ne.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))
v2.ns50.ne.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))
v2.ns100.ne.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(4))
set.seed(2184)
v2.ns20.ne.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))
v2.ns50.ne.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))
v2.ns100.ne.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))


## Save simulations
save(v2.ns20.ne.nullT,v2.ns50.ne.nullT,v2.ns100.ne.nullT,v2.ns20.ne.smalld,v2.ns50.ne.smalld,v2.ns100.ne.smalld,v2.ns20.ne.midd,v2.ns50.ne.midd,v2.ns100.ne.midd, v2.ns20.ne.bigd,v2.ns50.ne.bigd,v2.ns100.ne.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2NeSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2NeSeed2184.Rdata')
v2nesim1 <- list(v2.ns20.ne.nullT=v2.ns20.ne.nullT$sim.results$sim1$sim.data,
                 v2.ns50.ne.nullT=v2.ns50.ne.nullT$sim.results$sim1$sim.data,
                 v2.ns100.ne.nullT=v2.ns100.ne.nullT$sim.results$sim1$sim.data,
                 v2.ns20.ne.smalld=v2.ns20.ne.smalld$sim.results$sim1$sim.data,
                 v2.ns50.ne.smalld=v2.ns50.ne.smalld$sim.results$sim1$sim.data,
                 v2.ns100.ne.smalld=v2.ns100.ne.smalld$sim.results$sim1$sim.data,
                 v2.ns20.ne.midd=v2.ns20.ne.midd$sim.results$sim1$sim.data,
                 v2.ns50.ne.midd=v2.ns50.ne.midd$sim.results$sim1$sim.data,
                 v2.ns100.ne.midd=v2.ns100.ne.midd$sim.results$sim1$sim.data, 
                 v2.ns20.ne.bigd=v2.ns20.ne.bigd$sim.results$sim1$sim.data,
                 v2.ns50.ne.bigd=v2.ns50.ne.bigd$sim.results$sim1$sim.data,
                 v2.ns100.ne.bigd=v2.ns100.ne.bigd$sim.results$sim1$sim.data)
save(v2nesim1, file='/users/joshwondra/R-projects/Welch rule/v2NeSeed2184-sim1.Rdata')


##### Equal Ns, double variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v2.ne <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v2.ne[,1,1] <- c(v2.ns20.ne.nullT$classic.reject, v2.ns50.ne.nullT$classic.reject, v2.ns100.ne.nullT$classic.reject)
reject.null.v2.ne[,1,2] <- c(v2.ns20.ne.nullT$welch.reject, v2.ns50.ne.nullT$welch.reject, v2.ns100.ne.nullT$welch.reject)

# small effect
reject.null.v2.ne[,2,1] <- c(v2.ns20.ne.smalld$classic.reject, v2.ns50.ne.smalld$classic.reject, v2.ns100.ne.smalld$classic.reject)
reject.null.v2.ne[,2,2] <- c(v2.ns20.ne.smalld$welch.reject, v2.ns50.ne.smalld$welch.reject, v2.ns100.ne.smalld$welch.reject)

# medium effect
reject.null.v2.ne[,3,1] <- c(v2.ns20.ne.midd$classic.reject, v2.ns50.ne.midd$classic.reject, v2.ns100.ne.midd$classic.reject)
reject.null.v2.ne[,3,2] <- c(v2.ns20.ne.midd$welch.reject, v2.ns50.ne.midd$welch.reject, v2.ns100.ne.midd$welch.reject)

# large effect
reject.null.v2.ne[,4,1] <- c(v2.ns20.ne.bigd$classic.reject, v2.ns50.ne.bigd$classic.reject, v2.ns100.ne.bigd$classic.reject)
reject.null.v2.ne[,4,2] <- c(v2.ns20.ne.bigd$welch.reject, v2.ns50.ne.bigd$welch.reject, v2.ns100.ne.bigd$welch.reject)


##### Equal Ns, double variances coverage rate #####

## store observed coverage rates
obs.coverage.v2.ne <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v2.ne[,1,1] <- c(v2.ns20.ne.nullT$classic.coverage, v2.ns50.ne.nullT$classic.coverage, v2.ns100.ne.nullT$classic.coverage)
obs.coverage.v2.ne[,1,2] <- c(v2.ns20.ne.nullT$welch.coverage, v2.ns50.ne.nullT$welch.coverage, v2.ns100.ne.nullT$welch.coverage)

# small effect
obs.coverage.v2.ne[,2,1] <- c(v2.ns20.ne.smalld$classic.coverage, v2.ns50.ne.smalld$classic.coverage, v2.ns100.ne.smalld$classic.coverage)
obs.coverage.v2.ne[,2,2] <- c(v2.ns20.ne.smalld$welch.coverage, v2.ns50.ne.smalld$welch.coverage, v2.ns100.ne.smalld$welch.coverage)

# medium effect
obs.coverage.v2.ne[,3,1] <- c(v2.ns20.ne.midd$classic.coverage, v2.ns50.ne.midd$classic.coverage, v2.ns100.ne.midd$classic.coverage)
obs.coverage.v2.ne[,3,2] <- c(v2.ns20.ne.midd$welch.coverage, v2.ns50.ne.midd$welch.coverage, v2.ns100.ne.midd$welch.coverage)

# big effect
obs.coverage.v2.ne[,4,1] <- c(v2.ns20.ne.bigd$classic.coverage, v2.ns50.ne.bigd$classic.coverage, v2.ns100.ne.bigd$classic.coverage)
obs.coverage.v2.ne[,4,2] <- c(v2.ns20.ne.bigd$welch.coverage, v2.ns50.ne.bigd$welch.coverage, v2.ns100.ne.bigd$welch.coverage)


##### Equal Ns, double variances df proportion #####

## store df proportion
df.ratio.v2.ne <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v2.ne[,1,1] <- c(v2.ns20.ne.nullT$df.ratio.avg, v2.ns50.ne.nullT$df.ratio.avg, v2.ns100.ne.nullT$df.ratio.avg)  # true null
df.ratio.v2.ne[,2,1] <- c(v2.ns20.ne.smalld$df.ratio.avg, v2.ns50.ne.smalld$df.ratio.avg, v2.ns100.ne.smalld$df.ratio.avg)  # small d
df.ratio.v2.ne[,3,1] <- c(v2.ns20.ne.midd$df.ratio.avg, v2.ns50.ne.midd$df.ratio.avg, v2.ns100.ne.midd$df.ratio.avg)  # mid d
df.ratio.v2.ne[,4,1] <- c(v2.ns20.ne.bigd$df.ratio.avg, v2.ns50.ne.bigd$df.ratio.avg, v2.ns100.ne.bigd$df.ratio.avg)  # big d


##### Save tables from equal Ns, double variances simulations #####
save(df.ratio.v2.ne, obs.coverage.v2.ne, reject.null.v2.ne, file='/users/joshwondra/R-projects/Welch rule/v2NeSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v2NeSeed2184Tables.Rdata')



##### Equal Ns, 5x different variances simulations #####

# true null
set.seed(2184)
v5.ns20.ne.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns50.ne.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns100.ne.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,10), contrast=c(-1,1))

# small effect
cohen.diff(.2, 20, sqrt(2), 20, sqrt(10))
set.seed(2184)
v5.ns20.ne.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))
v5.ns50.ne.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))
v5.ns100.ne.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(10))
set.seed(2184)
v5.ns20.ne.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))
v5.ns50.ne.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))
v5.ns100.ne.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(10))
set.seed(2184)
v5.ns20.ne.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))
v5.ns50.ne.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))
v5.ns100.ne.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))


## Save simulations
save(v5.ns20.ne.nullT,v5.ns50.ne.nullT,v5.ns100.ne.nullT,v5.ns20.ne.smalld,v5.ns50.ne.smalld,v5.ns100.ne.smalld,v5.ns20.ne.midd,v5.ns50.ne.midd,v5.ns100.ne.midd, v5.ns20.ne.bigd,v5.ns50.ne.bigd,v5.ns100.ne.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5NeSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5NeSeed2184.Rdata')
v5nesim1 <- list(v5.ns20.ne.nullT=v5.ns20.ne.nullT$sim.results$sim1$sim.data,
                 v5.ns50.ne.nullT=v5.ns50.ne.nullT$sim.results$sim1$sim.data,
                 v5.ns100.ne.nullT=v5.ns100.ne.nullT$sim.results$sim1$sim.data,
                 v5.ns20.ne.smalld=v5.ns20.ne.smalld$sim.results$sim1$sim.data,
                 v5.ns50.ne.smalld=v5.ns50.ne.smalld$sim.results$sim1$sim.data,
                 v5.ns100.ne.smalld=v5.ns100.ne.smalld$sim.results$sim1$sim.data,
                 v5.ns20.ne.midd=v5.ns20.ne.midd$sim.results$sim1$sim.data,
                 v5.ns50.ne.midd=v5.ns50.ne.midd$sim.results$sim1$sim.data,
                 v5.ns100.ne.midd=v5.ns100.ne.midd$sim.results$sim1$sim.data, 
                 v5.ns20.ne.bigd=v5.ns20.ne.bigd$sim.results$sim1$sim.data,
                 v5.ns50.ne.bigd=v5.ns50.ne.bigd$sim.results$sim1$sim.data,
                 v5.ns100.ne.bigd=v5.ns100.ne.bigd$sim.results$sim1$sim.data)
save(v5nesim1, file='/users/joshwondra/R-projects/Welch rule/v5NeSeed2184-sim1.Rdata')


##### Equal Ns, 5x variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v5.ne <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v5.ne[,1,1] <- c(v5.ns20.ne.nullT$classic.reject, v5.ns50.ne.nullT$classic.reject, v5.ns100.ne.nullT$classic.reject)
reject.null.v5.ne[,1,2] <- c(v5.ns20.ne.nullT$welch.reject, v5.ns50.ne.nullT$welch.reject, v5.ns100.ne.nullT$welch.reject)

# small effect
reject.null.v5.ne[,2,1] <- c(v5.ns20.ne.smalld$classic.reject, v5.ns50.ne.smalld$classic.reject, v5.ns100.ne.smalld$classic.reject)
reject.null.v5.ne[,2,2] <- c(v5.ns20.ne.smalld$welch.reject, v5.ns50.ne.smalld$welch.reject, v5.ns100.ne.smalld$welch.reject)

# medium effect
reject.null.v5.ne[,3,1] <- c(v5.ns20.ne.midd$classic.reject, v5.ns50.ne.midd$classic.reject, v5.ns100.ne.midd$classic.reject)
reject.null.v5.ne[,3,2] <- c(v5.ns20.ne.midd$welch.reject, v5.ns50.ne.midd$welch.reject, v5.ns100.ne.midd$welch.reject)

# large effect
reject.null.v5.ne[,4,1] <- c(v5.ns20.ne.bigd$classic.reject, v5.ns50.ne.bigd$classic.reject, v5.ns100.ne.bigd$classic.reject)
reject.null.v5.ne[,4,2] <- c(v5.ns20.ne.bigd$welch.reject, v5.ns50.ne.bigd$welch.reject, v5.ns100.ne.bigd$welch.reject)


##### Equal Ns, 5x variances coverage rate #####

## store observed coverage rates
obs.coverage.v5.ne <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v5.ne[,1,1] <- c(v5.ns20.ne.nullT$classic.coverage, v5.ns50.ne.nullT$classic.coverage, v5.ns100.ne.nullT$classic.coverage)
obs.coverage.v5.ne[,1,2] <- c(v5.ns20.ne.nullT$welch.coverage, v5.ns50.ne.nullT$welch.coverage, v5.ns100.ne.nullT$welch.coverage)

# small effect
obs.coverage.v5.ne[,2,1] <- c(v5.ns20.ne.smalld$classic.coverage, v5.ns50.ne.smalld$classic.coverage, v5.ns100.ne.smalld$classic.coverage)
obs.coverage.v5.ne[,2,2] <- c(v5.ns20.ne.smalld$welch.coverage, v5.ns50.ne.smalld$welch.coverage, v5.ns100.ne.smalld$welch.coverage)

# medium effect
obs.coverage.v5.ne[,3,1] <- c(v5.ns20.ne.midd$classic.coverage, v5.ns50.ne.midd$classic.coverage, v5.ns100.ne.midd$classic.coverage)
obs.coverage.v5.ne[,3,2] <- c(v5.ns20.ne.midd$welch.coverage, v5.ns50.ne.midd$welch.coverage, v5.ns100.ne.midd$welch.coverage)

# big effect
obs.coverage.v5.ne[,4,1] <- c(v5.ns20.ne.bigd$classic.coverage, v5.ns50.ne.bigd$classic.coverage, v5.ns100.ne.bigd$classic.coverage)
obs.coverage.v5.ne[,4,2] <- c(v5.ns20.ne.bigd$welch.coverage, v5.ns50.ne.bigd$welch.coverage, v5.ns100.ne.bigd$welch.coverage)


##### Equal Ns, 5x variances df proportion #####

## store df proportion
df.ratio.v5.ne <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v5.ne[,1,1] <- c(v5.ns20.ne.nullT$df.ratio.avg, v5.ns50.ne.nullT$df.ratio.avg, v5.ns100.ne.nullT$df.ratio.avg)  # true null
df.ratio.v5.ne[,2,1] <- c(v5.ns20.ne.smalld$df.ratio.avg, v5.ns50.ne.smalld$df.ratio.avg, v5.ns100.ne.smalld$df.ratio.avg)  # small d
df.ratio.v5.ne[,3,1] <- c(v5.ns20.ne.midd$df.ratio.avg, v5.ns50.ne.midd$df.ratio.avg, v5.ns100.ne.midd$df.ratio.avg)  # mid d
df.ratio.v5.ne[,4,1] <- c(v5.ns20.ne.bigd$df.ratio.avg, v5.ns50.ne.bigd$df.ratio.avg, v5.ns100.ne.bigd$df.ratio.avg)  # big d



##### Save tables from equal Ns, 5x variances simulations #####
save(df.ratio.v5.ne, obs.coverage.v5.ne, reject.null.v5.ne, file='/users/joshwondra/R-projects/Welch rule/v5NeSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v5NeSeed2184Tables.Rdata')






##### Different Ns, equal variances simulations #####

##### 1.5x Ns, equal variances #####

# true null
set.seed(2184)
ve.ns20.1.5n.nullT <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns50.1.5n.nullT <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns100.1.5n.nullT <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6), vars=c(2,2), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(2),30,sqrt(2))
cohen.diff(.2,50,sqrt(2),75,sqrt(2))
cohen.diff(.2,100,sqrt(2),150,sqrt(2))
set.seed(2184)
ve.ns20.1.5n.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))
ve.ns50.1.5n.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))
ve.ns100.1.5n.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),30,sqrt(2))
cohen.diff(.5,50,sqrt(2),75,sqrt(2))
cohen.diff(.5,100,sqrt(2),150,sqrt(2))
set.seed(2184)
ve.ns20.1.5n.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))
ve.ns50.1.5n.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))
ve.ns100.1.5n.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),30,sqrt(2))
cohen.diff(.8,50,sqrt(2),75,sqrt(2))
cohen.diff(.8,100,sqrt(2),150,sqrt(2))
set.seed(2184)
ve.ns20.1.5n.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))
ve.ns50.1.5n.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))
ve.ns100.1.5n.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))

## Save simulations
save(ve.ns20.1.5n.nullT,ve.ns50.1.5n.nullT,ve.ns100.1.5n.nullT,ve.ns20.1.5n.smalld,ve.ns50.1.5n.smalld,ve.ns100.1.5n.smalld,ve.ns20.1.5n.midd,ve.ns50.1.5n.midd,ve.ns100.1.5n.midd,ve.ns20.1.5n.bigd,ve.ns50.1.5n.bigd,ve.ns100.1.5n.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/veN15Seed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/veN15Seed2184.Rdata')
ven15sim1 <- list(ve.ns20.1.5n.nullT=ve.ns20.1.5n.nullT$sim.results$sim1$sim.data,
               ve.ns50.1.5n.nullT=ve.ns50.1.5n.nullT$sim.results$sim1$sim.data,
               ve.ns100.1.5n.nullT=ve.ns100.1.5n.nullT$sim.results$sim1$sim.data,
               ve.ns20.1.5n.smalld=ve.ns20.1.5n.smalld$sim.results$sim1$sim.data,
               ve.ns50.1.5n.smalld=ve.ns50.1.5n.smalld$sim.results$sim1$sim.data,
               ve.ns100.1.5n.smalld=ve.ns100.1.5n.smalld$sim.results$sim1$sim.data,
               ve.ns20.1.5n.midd=ve.ns20.1.5n.midd$sim.results$sim1$sim.data,
               ve.ns50.1.5n.midd=ve.ns50.1.5n.midd$sim.results$sim1$sim.data,
               ve.ns100.1.5n.midd=ve.ns100.1.5n.midd$sim.results$sim1$sim.data,
               ve.ns20.1.5n.bigd=ve.ns20.1.5n.bigd$sim.results$sim1$sim.data,
               ve.ns50.1.5n.bigd=ve.ns50.1.5n.bigd$sim.results$sim1$sim.data,
               ve.ns100.1.5n.bigd=ve.ns100.1.5n.bigd$sim.results$sim1$sim.data)
save(ven15sim1, file='/users/joshwondra/R-projects/Welch rule/veN15Seed2184-sim1.Rdata')


##### 1.5x Ns, equal variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.ve.1.5n <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.ve.1.5n[,1,1] <- c(ve.ns20.1.5n.nullT$classic.reject, ve.ns50.1.5n.nullT$classic.reject, ve.ns100.1.5n.nullT$classic.reject)
reject.null.ve.1.5n[,1,2] <- c(ve.ns20.1.5n.nullT$welch.reject, ve.ns50.1.5n.nullT$welch.reject, ve.ns100.1.5n.nullT$welch.reject)

# small effect
reject.null.ve.1.5n[,2,1] <- c(ve.ns20.1.5n.smalld$classic.reject, ve.ns50.1.5n.smalld$classic.reject, ve.ns100.1.5n.smalld$classic.reject)
reject.null.ve.1.5n[,2,2] <- c(ve.ns20.1.5n.smalld$welch.reject, ve.ns50.1.5n.smalld$welch.reject, ve.ns100.1.5n.smalld$welch.reject)

# medium effect
reject.null.ve.1.5n[,3,1] <- c(ve.ns20.1.5n.midd$classic.reject, ve.ns50.1.5n.midd$classic.reject, ve.ns100.1.5n.midd$classic.reject)
reject.null.ve.1.5n[,3,2] <- c(ve.ns20.1.5n.midd$welch.reject, ve.ns50.1.5n.midd$welch.reject, ve.ns100.1.5n.midd$welch.reject)

# large effect
reject.null.ve.1.5n[,4,1] <- c(ve.ns20.1.5n.bigd$classic.reject, ve.ns50.1.5n.bigd$classic.reject, ve.ns100.1.5n.bigd$classic.reject)
reject.null.ve.1.5n[,4,2] <- c(ve.ns20.1.5n.bigd$welch.reject, ve.ns50.1.5n.bigd$welch.reject, ve.ns100.1.5n.bigd$welch.reject)


##### 1.5x Ns, equal variances coverage rate #####

## store observed coverage rates
obs.coverage.ve.1.5n <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.ve.1.5n[,1,1] <- c(ve.ns20.1.5n.nullT$classic.coverage, ve.ns50.1.5n.nullT$classic.coverage, ve.ns100.1.5n.nullT$classic.coverage)
obs.coverage.ve.1.5n[,1,2] <- c(ve.ns20.1.5n.nullT$welch.coverage, ve.ns50.1.5n.nullT$welch.coverage, ve.ns100.1.5n.nullT$welch.coverage)

# small effect
obs.coverage.ve.1.5n[,2,1] <- c(ve.ns20.1.5n.smalld$classic.coverage, ve.ns50.1.5n.smalld$classic.coverage, ve.ns100.1.5n.smalld$classic.coverage)
obs.coverage.ve.1.5n[,2,2] <- c(ve.ns20.1.5n.smalld$welch.coverage, ve.ns50.1.5n.smalld$welch.coverage, ve.ns100.1.5n.smalld$welch.coverage)

# medium effect
obs.coverage.ve.1.5n[,3,1] <- c(ve.ns20.1.5n.midd$classic.coverage, ve.ns50.1.5n.midd$classic.coverage, ve.ns100.1.5n.midd$classic.coverage)
obs.coverage.ve.1.5n[,3,2] <- c(ve.ns20.1.5n.midd$welch.coverage, ve.ns50.1.5n.midd$welch.coverage, ve.ns100.1.5n.midd$welch.coverage)

# big effect
obs.coverage.ve.1.5n[,4,1] <- c(ve.ns20.1.5n.bigd$classic.coverage, ve.ns50.1.5n.bigd$classic.coverage, ve.ns100.1.5n.bigd$classic.coverage)
obs.coverage.ve.1.5n[,4,2] <- c(ve.ns20.1.5n.bigd$welch.coverage, ve.ns50.1.5n.bigd$welch.coverage, ve.ns100.1.5n.bigd$welch.coverage)



##### 1.5 Ns, equal variances df proportion #####

## store df proportion
df.ratio.ve.1.5n <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.ve.1.5n[,1,1] <- c(ve.ns20.1.5n.nullT$df.ratio.avg, ve.ns50.1.5n.nullT$df.ratio.avg, ve.ns100.1.5n.nullT$df.ratio.avg)  # true null
df.ratio.ve.1.5n[,2,1] <- c(ve.ns20.1.5n.smalld$df.ratio.avg, ve.ns50.1.5n.smalld$df.ratio.avg, ve.ns100.1.5n.smalld$df.ratio.avg)  # small d
df.ratio.ve.1.5n[,3,1] <- c(ve.ns20.1.5n.midd$df.ratio.avg, ve.ns50.1.5n.midd$df.ratio.avg, ve.ns100.1.5n.midd$df.ratio.avg)  # mid d
df.ratio.ve.1.5n[,4,1] <- c(ve.ns20.1.5n.bigd$df.ratio.avg, ve.ns50.1.5n.bigd$df.ratio.avg, ve.ns100.1.5n.bigd$df.ratio.avg)  # big d



##### Save tables from 1.5x Ns, equal variances simulations #####
save(df.ratio.ve.1.5n, obs.coverage.ve.1.5n, reject.null.ve.1.5n, file='/users/joshwondra/R-projects/Welch rule/ve1andhalfnSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/ve1andhalfnSeed2184Tables.Rdata')




##### Double Ns, equal variances #####

# true null
set.seed(2184)
ve.ns20.2n.nullT <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns50.2n.nullT <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns100.2n.nullT <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6), vars=c(2,2), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(2),40,sqrt(2))
cohen.diff(.2,50,sqrt(2),100,sqrt(2))
cohen.diff(.2,100,sqrt(2),200,sqrt(2))
set.seed(2184)
ve.ns20.2n.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))
ve.ns50.2n.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))
ve.ns100.2n.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.28), vars=c(2,2), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),40,sqrt(2))
cohen.diff(.5,50,sqrt(2),100,sqrt(2))
cohen.diff(.5,100,sqrt(2),200,sqrt(2))
set.seed(2184)
ve.ns20.2n.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))
ve.ns50.2n.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))
ve.ns100.2n.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.71), vars=c(2,2), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),40,sqrt(2))
cohen.diff(.8,50,sqrt(2),100,sqrt(2))
cohen.diff(.8,100,sqrt(2),200,sqrt(2))
set.seed(2184)
ve.ns20.2n.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))
ve.ns50.2n.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))
ve.ns100.2n.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.13), vars=c(2,2), contrast=c(-1,1))


## Save simulations
save(ve.ns20.2n.nullT,ve.ns50.2n.nullT,ve.ns100.2n.nullT,ve.ns20.2n.smalld,ve.ns50.2n.smalld,ve.ns100.2n.smalld,ve.ns20.2n.midd,ve.ns50.2n.midd,ve.ns100.2n.midd, ve.ns20.2n.bigd,ve.ns50.2n.bigd,ve.ns100.2n.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/veN2Seed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/veN2Seed2184.Rdata')

ven2sim1 <- list(ve.ns20.2n.nullT=ve.ns20.2n.nullT$sim.results$sim1$sim.data,
                 ve.ns50.2n.nullT=ve.ns50.2n.nullT$sim.results$sim1$sim.data,
                 ve.ns100.2n.nullT=ve.ns100.2n.nullT$sim.results$sim1$sim.data,
                 ve.ns20.2n.smalld=ve.ns20.2n.smalld$sim.results$sim1$sim.data,
                 ve.ns50.2n.smalld=ve.ns50.2n.smalld$sim.results$sim1$sim.data,
                 ve.ns100.2n.smalld=ve.ns100.2n.smalld$sim.results$sim1$sim.data,
                 ve.ns20.2n.midd=ve.ns20.2n.midd$sim.results$sim1$sim.data,
                 ve.ns50.2n.midd=ve.ns50.2n.midd$sim.results$sim1$sim.data,
                 ve.ns100.2n.midd=ve.ns100.2n.midd$sim.results$sim1$sim.data, 
                 ve.ns20.2n.bigd=ve.ns20.2n.bigd$sim.results$sim1$sim.data,
                 ve.ns50.2n.bigd=ve.ns50.2n.bigd$sim.results$sim1$sim.data,
                 ve.ns100.2n.bigd=ve.ns100.2n.bigd$sim.results$sim1$sim.data)
save(ven2sim1, file='/users/joshwondra/R-projects/Welch rule/veN2Seed2184-sim1.Rdata')


##### Double Ns, equal variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.ve.2n <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.ve.2n[,1,1] <- c(ve.ns20.2n.nullT$classic.reject, ve.ns50.2n.nullT$classic.reject, ve.ns100.2n.nullT$classic.reject)
reject.null.ve.2n[,1,2] <- c(ve.ns20.2n.nullT$welch.reject, ve.ns50.2n.nullT$welch.reject, ve.ns100.2n.nullT$welch.reject)

# small effect
reject.null.ve.2n[,2,1] <- c(ve.ns20.2n.smalld$classic.reject, ve.ns50.2n.smalld$classic.reject, ve.ns100.2n.smalld$classic.reject)
reject.null.ve.2n[,2,2] <- c(ve.ns20.2n.smalld$welch.reject, ve.ns50.2n.smalld$welch.reject, ve.ns100.2n.smalld$welch.reject)

# medium effect
reject.null.ve.2n[,3,1] <- c(ve.ns20.2n.midd$classic.reject, ve.ns50.2n.midd$classic.reject, ve.ns100.2n.midd$classic.reject)
reject.null.ve.2n[,3,2] <- c(ve.ns20.2n.midd$welch.reject, ve.ns50.2n.midd$welch.reject, ve.ns100.2n.midd$welch.reject)

# large effect
reject.null.ve.2n[,4,1] <- c(ve.ns20.2n.bigd$classic.reject, ve.ns50.2n.bigd$classic.reject, ve.ns100.2n.bigd$classic.reject)
reject.null.ve.2n[,4,2] <- c(ve.ns20.2n.bigd$welch.reject, ve.ns50.2n.bigd$welch.reject, ve.ns100.2n.bigd$welch.reject)


##### Double Ns, equal variances coverage rate #####

## store observed coverage rates
obs.coverage.ve.2n <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.ve.2n[,1,1] <- c(ve.ns20.2n.nullT$classic.coverage, ve.ns50.2n.nullT$classic.coverage, ve.ns100.2n.nullT$classic.coverage)
obs.coverage.ve.2n[,1,2] <- c(ve.ns20.2n.nullT$welch.coverage, ve.ns50.2n.nullT$welch.coverage, ve.ns100.2n.nullT$welch.coverage)

# small effect
obs.coverage.ve.2n[,2,1] <- c(ve.ns20.2n.smalld$classic.coverage, ve.ns50.2n.smalld$classic.coverage, ve.ns100.2n.smalld$classic.coverage)
obs.coverage.ve.2n[,2,2] <- c(ve.ns20.2n.smalld$welch.coverage, ve.ns50.2n.smalld$welch.coverage, ve.ns100.2n.smalld$welch.coverage)

# medium effect
obs.coverage.ve.2n[,3,1] <- c(ve.ns20.2n.midd$classic.coverage, ve.ns50.2n.midd$classic.coverage, ve.ns100.2n.midd$classic.coverage)
obs.coverage.ve.2n[,3,2] <- c(ve.ns20.2n.midd$welch.coverage, ve.ns50.2n.midd$welch.coverage, ve.ns100.2n.midd$welch.coverage)

# big effect
obs.coverage.ve.2n[,4,1] <- c(ve.ns20.2n.bigd$classic.coverage, ve.ns50.2n.bigd$classic.coverage, ve.ns100.2n.bigd$classic.coverage)
obs.coverage.ve.2n[,4,2] <- c(ve.ns20.2n.bigd$welch.coverage, ve.ns50.2n.bigd$welch.coverage, ve.ns100.2n.bigd$welch.coverage)



##### Double Ns, equal variances df proportion #####

## store df proportion
df.ratio.ve.2n <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.ve.2n[,1,1] <- c(ve.ns20.2n.nullT$df.ratio.avg, ve.ns50.2n.nullT$df.ratio.avg, ve.ns100.2n.nullT$df.ratio.avg)  # true null
df.ratio.ve.2n[,2,1] <- c(ve.ns20.2n.smalld$df.ratio.avg, ve.ns50.2n.smalld$df.ratio.avg, ve.ns100.2n.smalld$df.ratio.avg)  # small d
df.ratio.ve.2n[,3,1] <- c(ve.ns20.2n.midd$df.ratio.avg, ve.ns50.2n.midd$df.ratio.avg, ve.ns100.2n.midd$df.ratio.avg)  # mid d
df.ratio.ve.2n[,4,1] <- c(ve.ns20.2n.bigd$df.ratio.avg, ve.ns50.2n.bigd$df.ratio.avg, ve.ns100.2n.bigd$df.ratio.avg)  # big d




##### Save tables from double Ns, equal variances simulations #####
save(df.ratio.ve.2n, obs.coverage.ve.2n, reject.null.ve.2n, file='/users/joshwondra/R-projects/Welch rule/ve2nSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/ve2nSeed2184Tables.Rdata')






##### Different Ns and variances, small group small variances #####

##### 1.5x Ns, double variances simulations #####

# true null
set.seed(2184)
v2.ns20.1.5n.ssv.nullT <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns50.1.5n.ssv.nullT <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns100.1.5n.ssv.nullT <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6), vars=c(2,4), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(2),30,sqrt(4))
cohen.diff(.2,50,sqrt(2),75,sqrt(4))
cohen.diff(.2,100,sqrt(2),150,sqrt(4))
set.seed(2184)
v2.ns20.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))
v2.ns50.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))
v2.ns100.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),30,sqrt(4))
cohen.diff(.5,50,sqrt(2),75,sqrt(4))
cohen.diff(.5,100,sqrt(2),150,sqrt(4))
set.seed(2184)
v2.ns20.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))
v2.ns50.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))
v2.ns100.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),30,sqrt(4))
cohen.diff(.8,50,sqrt(2),75,sqrt(4))
cohen.diff(.8,100,sqrt(2),150,sqrt(4))
set.seed(2184)
v2.ns20.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))
v2.ns50.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))
v2.ns100.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))


## Save simulations
save(v2.ns20.1.5n.ssv.nullT,v2.ns50.1.5n.ssv.nullT,v2.ns100.1.5n.ssv.nullT,v2.ns20.1.5n.ssv.smalld,v2.ns50.1.5n.ssv.smalld,v2.ns100.1.5n.ssv.smalld,v2.ns20.1.5n.ssv.midd,v2.ns50.1.5n.ssv.midd,v2.ns100.1.5n.ssv.midd, v2.ns20.1.5n.ssv.bigd,v2.ns50.1.5n.ssv.bigd,v2.ns100.1.5n.ssv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N15ssvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N15ssvSeed2184.Rdata')

v2n15ssvsim1 <- list(v2.ns20.1.5n.ssv.nullT=v2.ns20.1.5n.ssv.nullT$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.ssv.nullT=v2.ns50.1.5n.ssv.nullT$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.ssv.nullT=v2.ns100.1.5n.ssv.nullT$sim.results$sim1$sim.data,
                     v2.ns20.1.5n.ssv.smalld=v2.ns20.1.5n.ssv.smalld$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.ssv.smalld=v2.ns50.1.5n.ssv.smalld$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.ssv.smalld=v2.ns100.1.5n.ssv.smalld$sim.results$sim1$sim.data,
                     v2.ns20.1.5n.ssv.midd=v2.ns20.1.5n.ssv.midd$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.ssv.midd=v2.ns50.1.5n.ssv.midd$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.ssv.midd=v2.ns100.1.5n.ssv.midd$sim.results$sim1$sim.data,
                     v2.ns20.1.5n.ssv.bigd=v2.ns20.1.5n.ssv.bigd$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.ssv.bigd=v2.ns50.1.5n.ssv.bigd$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.ssv.bigd=v2.ns100.1.5n.ssv.bigd$sim.results$sim1$sim.data)
save(v2n15ssvsim1, file='/users/joshwondra/R-projects/Welch rule/v2N15ssvSeed2184-sim1.Rdata')


##### 1.5x Ns, double variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v2.1.5n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v2.1.5n.ssv[,1,1] <- c(v2.ns20.1.5n.ssv.nullT$classic.reject, v2.ns50.1.5n.ssv.nullT$classic.reject, v2.ns100.1.5n.ssv.nullT$classic.reject)
reject.null.v2.1.5n.ssv[,1,2] <- c(v2.ns20.1.5n.ssv.nullT$welch.reject, v2.ns50.1.5n.ssv.nullT$welch.reject, v2.ns100.1.5n.ssv.nullT$welch.reject)

# small effect
reject.null.v2.1.5n.ssv[,2,1] <- c(v2.ns20.1.5n.ssv.smalld$classic.reject, v2.ns50.1.5n.ssv.smalld$classic.reject, v2.ns100.1.5n.ssv.smalld$classic.reject)
reject.null.v2.1.5n.ssv[,2,2] <- c(v2.ns20.1.5n.ssv.smalld$welch.reject, v2.ns50.1.5n.ssv.smalld$welch.reject, v2.ns100.1.5n.ssv.smalld$welch.reject)

# medium effect
reject.null.v2.1.5n.ssv[,3,1] <- c(v2.ns20.1.5n.ssv.midd$classic.reject, v2.ns50.1.5n.ssv.midd$classic.reject, v2.ns100.1.5n.ssv.midd$classic.reject)
reject.null.v2.1.5n.ssv[,3,2] <- c(v2.ns20.1.5n.ssv.midd$welch.reject, v2.ns50.1.5n.ssv.midd$welch.reject, v2.ns100.1.5n.ssv.midd$welch.reject)

# large effect
reject.null.v2.1.5n.ssv[,4,1] <- c(v2.ns20.1.5n.ssv.bigd$classic.reject, v2.ns50.1.5n.ssv.bigd$classic.reject, v2.ns100.1.5n.ssv.bigd$classic.reject)
reject.null.v2.1.5n.ssv[,4,2] <- c(v2.ns20.1.5n.ssv.bigd$welch.reject, v2.ns50.1.5n.ssv.bigd$welch.reject, v2.ns100.1.5n.ssv.bigd$welch.reject)


##### 1.5x Ns, double variances coverage rate #####

## store observed coverage rates
obs.coverage.v2.1.5n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v2.1.5n.ssv[,1,1] <- c(v2.ns20.1.5n.ssv.nullT$classic.coverage, v2.ns50.1.5n.ssv.nullT$classic.coverage, v2.ns100.1.5n.ssv.nullT$classic.coverage)
obs.coverage.v2.1.5n.ssv[,1,2] <- c(v2.ns20.1.5n.ssv.nullT$welch.coverage, v2.ns50.1.5n.ssv.nullT$welch.coverage, v2.ns100.1.5n.ssv.nullT$welch.coverage)

# small effect
obs.coverage.v2.1.5n.ssv[,2,1] <- c(v2.ns20.1.5n.ssv.smalld$classic.coverage, v2.ns50.1.5n.ssv.smalld$classic.coverage, v2.ns100.1.5n.ssv.smalld$classic.coverage)
obs.coverage.v2.1.5n.ssv[,2,2] <- c(v2.ns20.1.5n.ssv.smalld$welch.coverage, v2.ns50.1.5n.ssv.smalld$welch.coverage, v2.ns100.1.5n.ssv.smalld$welch.coverage)

# medium effect
obs.coverage.v2.1.5n.ssv[,3,1] <- c(v2.ns20.1.5n.ssv.midd$classic.coverage, v2.ns50.1.5n.ssv.midd$classic.coverage, v2.ns100.1.5n.ssv.midd$classic.coverage)
obs.coverage.v2.1.5n.ssv[,3,2] <- c(v2.ns20.1.5n.ssv.midd$welch.coverage, v2.ns50.1.5n.ssv.midd$welch.coverage, v2.ns100.1.5n.ssv.midd$welch.coverage)

# big effect
obs.coverage.v2.1.5n.ssv[,4,1] <- c(v2.ns20.1.5n.ssv.bigd$classic.coverage, v2.ns50.1.5n.ssv.bigd$classic.coverage, v2.ns100.1.5n.ssv.bigd$classic.coverage)
obs.coverage.v2.1.5n.ssv[,4,2] <- c(v2.ns20.1.5n.ssv.bigd$welch.coverage, v2.ns50.1.5n.ssv.bigd$welch.coverage, v2.ns100.1.5n.ssv.bigd$welch.coverage)



##### 1.5x Ns, double variances df proportion #####

## store df proportion
df.ratio.v2.1.5n.ssv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v2.1.5n.ssv[,1,1] <- c(v2.ns20.1.5n.ssv.nullT$df.ratio.avg, v2.ns50.1.5n.ssv.nullT$df.ratio.avg, v2.ns100.1.5n.ssv.nullT$df.ratio.avg)  # true null
df.ratio.v2.1.5n.ssv[,2,1] <- c(v2.ns20.1.5n.ssv.smalld$df.ratio.avg, v2.ns50.1.5n.ssv.smalld$df.ratio.avg, v2.ns100.1.5n.ssv.smalld$df.ratio.avg)  # small d
df.ratio.v2.1.5n.ssv[,3,1] <- c(v2.ns20.1.5n.ssv.midd$df.ratio.avg, v2.ns50.1.5n.ssv.midd$df.ratio.avg, v2.ns100.1.5n.ssv.midd$df.ratio.avg)  # mid d
df.ratio.v2.1.5n.ssv[,4,1] <- c(v2.ns20.1.5n.ssv.bigd$df.ratio.avg, v2.ns50.1.5n.ssv.bigd$df.ratio.avg, v2.ns100.1.5n.ssv.bigd$df.ratio.avg)  # big d



##### Save tables from 1.5x Ns, double variances simulations #####
save(df.ratio.v2.1.5n.ssv, obs.coverage.v2.1.5n.ssv, reject.null.v2.1.5n.ssv, file='/users/joshwondra/R-projects/Welch rule/v21andhalfnSSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v21andhalfnSSVSeed2184Tables.Rdata')











##### 1.5x Ns, 5x variances simulations #####

# true null
set.seed(2184)
v5.ns20.1.5n.ssv.nullT <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.nullT <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.nullT <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6), vars=c(2,10), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(2),30,sqrt(10))
cohen.diff(.2,50,sqrt(2),75,sqrt(10))
cohen.diff(.2,100,sqrt(2),150,sqrt(10))
set.seed(2184)
v5.ns20.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),30,sqrt(10))
cohen.diff(.5,50,sqrt(2),75,sqrt(10))
cohen.diff(.5,100,sqrt(2),150,sqrt(10))
set.seed(2184)
v5.ns20.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),30,sqrt(10))
cohen.diff(.8,50,sqrt(2),75,sqrt(10))
cohen.diff(.8,100,sqrt(2),150,sqrt(10))
set.seed(2184)
v5.ns20.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))


## Save simulations
save(v5.ns20.1.5n.ssv.nullT,v5.ns50.1.5n.ssv.nullT,v5.ns100.1.5n.ssv.nullT,v5.ns20.1.5n.ssv.smalld,v5.ns50.1.5n.ssv.smalld,v5.ns100.1.5n.ssv.smalld,v5.ns20.1.5n.ssv.midd,v5.ns50.1.5n.ssv.midd,v5.ns100.1.5n.ssv.midd, v5.ns20.1.5n.ssv.bigd,v5.ns50.1.5n.ssv.bigd,v5.ns100.1.5n.ssv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N15ssvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N15ssvSeed2184.Rdata')

v5n15ssvsim1 <- list(v5.ns20.1.5n.ssv.nullT=v5.ns20.1.5n.ssv.nullT$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.ssv.nullT=v5.ns50.1.5n.ssv.nullT$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.ssv.nullT=v5.ns100.1.5n.ssv.nullT$sim.results$sim1$sim.data,
                     v5.ns20.1.5n.ssv.smalld=v5.ns20.1.5n.ssv.smalld$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.ssv.smalld=v5.ns50.1.5n.ssv.smalld$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.ssv.smalld=v5.ns100.1.5n.ssv.smalld$sim.results$sim1$sim.data,
                     v5.ns20.1.5n.ssv.midd=v5.ns20.1.5n.ssv.midd$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.ssv.midd=v5.ns50.1.5n.ssv.midd$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.ssv.midd=v5.ns100.1.5n.ssv.midd$sim.results$sim1$sim.data,
                     v5.ns20.1.5n.ssv.bigd=v5.ns20.1.5n.ssv.bigd$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.ssv.bigd=v5.ns50.1.5n.ssv.bigd$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.ssv.bigd=v5.ns100.1.5n.ssv.bigd$sim.results$sim1$sim.data)
save(v5n15ssvsim1, file='/users/joshwondra/R-projects/Welch rule/v5N15ssvSeed2184-sim1.Rdata')


##### 1.5x Ns, 5x variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v5.1.5n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v5.1.5n.ssv[,1,1] <- c(v5.ns20.1.5n.ssv.nullT$classic.reject, v5.ns50.1.5n.ssv.nullT$classic.reject, v5.ns100.1.5n.ssv.nullT$classic.reject)
reject.null.v5.1.5n.ssv[,1,2] <- c(v5.ns20.1.5n.ssv.nullT$welch.reject, v5.ns50.1.5n.ssv.nullT$welch.reject, v5.ns100.1.5n.ssv.nullT$welch.reject)

# small effect
reject.null.v5.1.5n.ssv[,2,1] <- c(v5.ns20.1.5n.ssv.smalld$classic.reject, v5.ns50.1.5n.ssv.smalld$classic.reject, v5.ns100.1.5n.ssv.smalld$classic.reject)
reject.null.v5.1.5n.ssv[,2,2] <- c(v5.ns20.1.5n.ssv.smalld$welch.reject, v5.ns50.1.5n.ssv.smalld$welch.reject, v5.ns100.1.5n.ssv.smalld$welch.reject)

# medium effect
reject.null.v5.1.5n.ssv[,3,1] <- c(v5.ns20.1.5n.ssv.midd$classic.reject, v5.ns50.1.5n.ssv.midd$classic.reject, v5.ns100.1.5n.ssv.midd$classic.reject)
reject.null.v5.1.5n.ssv[,3,2] <- c(v5.ns20.1.5n.ssv.midd$welch.reject, v5.ns50.1.5n.ssv.midd$welch.reject, v5.ns100.1.5n.ssv.midd$welch.reject)

# large effect
reject.null.v5.1.5n.ssv[,4,1] <- c(v5.ns20.1.5n.ssv.bigd$classic.reject, v5.ns50.1.5n.ssv.bigd$classic.reject, v5.ns100.1.5n.ssv.bigd$classic.reject)
reject.null.v5.1.5n.ssv[,4,2] <- c(v5.ns20.1.5n.ssv.bigd$welch.reject, v5.ns50.1.5n.ssv.bigd$welch.reject, v5.ns100.1.5n.ssv.bigd$welch.reject)


##### 1.5x Ns, 5x variances coverage rate #####

## store observed coverage rates
obs.coverage.v5.1.5n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v5.1.5n.ssv[,1,1] <- c(v5.ns20.1.5n.ssv.nullT$classic.coverage, v5.ns50.1.5n.ssv.nullT$classic.coverage, v5.ns100.1.5n.ssv.nullT$classic.coverage)
obs.coverage.v5.1.5n.ssv[,1,2] <- c(v5.ns20.1.5n.ssv.nullT$welch.coverage, v5.ns50.1.5n.ssv.nullT$welch.coverage, v5.ns100.1.5n.ssv.nullT$welch.coverage)

# small effect
obs.coverage.v5.1.5n.ssv[,2,1] <- c(v5.ns20.1.5n.ssv.smalld$classic.coverage, v5.ns50.1.5n.ssv.smalld$classic.coverage, v5.ns100.1.5n.ssv.smalld$classic.coverage)
obs.coverage.v5.1.5n.ssv[,2,2] <- c(v5.ns20.1.5n.ssv.smalld$welch.coverage, v5.ns50.1.5n.ssv.smalld$welch.coverage, v5.ns100.1.5n.ssv.smalld$welch.coverage)

# medium effect
obs.coverage.v5.1.5n.ssv[,3,1] <- c(v5.ns20.1.5n.ssv.midd$classic.coverage, v5.ns50.1.5n.ssv.midd$classic.coverage, v5.ns100.1.5n.ssv.midd$classic.coverage)
obs.coverage.v5.1.5n.ssv[,3,2] <- c(v5.ns20.1.5n.ssv.midd$welch.coverage, v5.ns50.1.5n.ssv.midd$welch.coverage, v5.ns100.1.5n.ssv.midd$welch.coverage)

# big effect
obs.coverage.v5.1.5n.ssv[,4,1] <- c(v5.ns20.1.5n.ssv.bigd$classic.coverage, v5.ns50.1.5n.ssv.bigd$classic.coverage, v5.ns100.1.5n.ssv.bigd$classic.coverage)
obs.coverage.v5.1.5n.ssv[,4,2] <- c(v5.ns20.1.5n.ssv.bigd$welch.coverage, v5.ns50.1.5n.ssv.bigd$welch.coverage, v5.ns100.1.5n.ssv.bigd$welch.coverage)



##### 1.5x Ns, 5x variances df proportion #####

## store df proportion
df.ratio.v5.1.5n.ssv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v5.1.5n.ssv[,1,1] <- c(v5.ns20.1.5n.ssv.nullT$df.ratio.avg, v5.ns50.1.5n.ssv.nullT$df.ratio.avg, v5.ns100.1.5n.ssv.nullT$df.ratio.avg)  # true null
df.ratio.v5.1.5n.ssv[,2,1] <- c(v5.ns20.1.5n.ssv.smalld$df.ratio.avg, v5.ns50.1.5n.ssv.smalld$df.ratio.avg, v5.ns100.1.5n.ssv.smalld$df.ratio.avg)  # small d
df.ratio.v5.1.5n.ssv[,3,1] <- c(v5.ns20.1.5n.ssv.midd$df.ratio.avg, v5.ns50.1.5n.ssv.midd$df.ratio.avg, v5.ns100.1.5n.ssv.midd$df.ratio.avg)  # mid d
df.ratio.v5.1.5n.ssv[,4,1] <- c(v5.ns20.1.5n.ssv.bigd$df.ratio.avg, v5.ns50.1.5n.ssv.bigd$df.ratio.avg, v5.ns100.1.5n.ssv.bigd$df.ratio.avg)  # big d





##### Save tables from 1.5x Ns, 5x variances simulations #####
save(df.ratio.v5.1.5n.ssv, obs.coverage.v5.1.5n.ssv, reject.null.v5.1.5n.ssv, file='/users/joshwondra/R-projects/Welch rule/v51andhalfnSSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v51andhalfnSSVSeed2184Tables.Rdata')





##### Double Ns, double variances simulations #####

# true null
set.seed(2184)
v2.ns20.2n.ssv.nullT <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.nullT <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.nullT <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6), vars=c(2,4), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(2),40,sqrt(4))
cohen.diff(.2,50,sqrt(2),100,sqrt(4))
cohen.diff(.2,100,sqrt(2),200,sqrt(4))
set.seed(2184)
v2.ns20.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.28), vars=c(2,4), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),40,sqrt(4))
cohen.diff(.5,50,sqrt(2),100,sqrt(4))
cohen.diff(.5,100,sqrt(2),200,sqrt(4))
set.seed(2184)
v2.ns20.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.71), vars=c(2,4), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),40,sqrt(4))
cohen.diff(.8,50,sqrt(2),100,sqrt(4))
cohen.diff(.8,100,sqrt(2),200,sqrt(4))
set.seed(2184)
v2.ns20.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.13), vars=c(2,4), contrast=c(-1,1))

# Save simulations
save(v2.ns20.2n.ssv.nullT,v2.ns50.2n.ssv.nullT,v2.ns100.2n.ssv.nullT,v2.ns20.2n.ssv.smalld,v2.ns50.2n.ssv.smalld,v2.ns100.2n.ssv.smalld,v2.ns20.2n.ssv.midd,v2.ns50.2n.ssv.midd,v2.ns100.2n.ssv.midd, v2.ns20.2n.ssv.bigd,v2.ns50.2n.ssv.bigd,v2.ns100.2n.ssv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N2ssvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N2ssvSeed2184.Rdata')

v2n2ssvsim1 <- list(v2.ns20.2n.ssv.nullT=v2.ns20.2n.ssv.nullT$sim.results$sim1$sim.data,
                    v2.ns50.2n.ssv.nullT=v2.ns50.2n.ssv.nullT$sim.results$sim1$sim.data,
                    v2.ns100.2n.ssv.nullT=v2.ns100.2n.ssv.nullT$sim.results$sim1$sim.data,
                    v2.ns20.2n.ssv.smalld=v2.ns20.2n.ssv.smalld$sim.results$sim1$sim.data,
                    v2.ns50.2n.ssv.smalld=v2.ns50.2n.ssv.smalld$sim.results$sim1$sim.data,
                    v2.ns100.2n.ssv.smalld=v2.ns100.2n.ssv.smalld$sim.results$sim1$sim.data,
                    v2.ns20.2n.ssv.midd=v2.ns20.2n.ssv.midd$sim.results$sim1$sim.data,
                    v2.ns50.2n.ssv.midd=v2.ns50.2n.ssv.midd$sim.results$sim1$sim.data,
                    v2.ns100.2n.ssv.midd=v2.ns100.2n.ssv.midd$sim.results$sim1$sim.data,
                    v2.ns20.2n.ssv.bigd=v2.ns20.2n.ssv.bigd$sim.results$sim1$sim.data,
                    v2.ns50.2n.ssv.bigd=v2.ns50.2n.ssv.bigd$sim.results$sim1$sim.data,
                    v2.ns100.2n.ssv.bigd=v2.ns100.2n.ssv.bigd$sim.results$sim1$sim.data)
save(v2n2ssvsim1, file='/users/joshwondra/R-projects/Welch rule/v2N2ssvSeed2184-sim1.Rdata')



##### Double Ns, double variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v2.2n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v2.2n.ssv[,1,1] <- c(v2.ns20.2n.ssv.nullT$classic.reject, v2.ns50.2n.ssv.nullT$classic.reject, v2.ns100.2n.ssv.nullT$classic.reject)
reject.null.v2.2n.ssv[,1,2] <- c(v2.ns20.2n.ssv.nullT$welch.reject, v2.ns50.2n.ssv.nullT$welch.reject, v2.ns100.2n.ssv.nullT$welch.reject)

# small effect
reject.null.v2.2n.ssv[,2,1] <- c(v2.ns20.2n.ssv.smalld$classic.reject, v2.ns50.2n.ssv.smalld$classic.reject, v2.ns100.2n.ssv.smalld$classic.reject)
reject.null.v2.2n.ssv[,2,2] <- c(v2.ns20.2n.ssv.smalld$welch.reject, v2.ns50.2n.ssv.smalld$welch.reject, v2.ns100.2n.ssv.smalld$welch.reject)

# medium effect
reject.null.v2.2n.ssv[,3,1] <- c(v2.ns20.2n.ssv.midd$classic.reject, v2.ns50.2n.ssv.midd$classic.reject, v2.ns100.2n.ssv.midd$classic.reject)
reject.null.v2.2n.ssv[,3,2] <- c(v2.ns20.2n.ssv.midd$welch.reject, v2.ns50.2n.ssv.midd$welch.reject, v2.ns100.2n.ssv.midd$welch.reject)

# large effect
reject.null.v2.2n.ssv[,4,1] <- c(v2.ns20.2n.ssv.bigd$classic.reject, v2.ns50.2n.ssv.bigd$classic.reject, v2.ns100.2n.ssv.bigd$classic.reject)
reject.null.v2.2n.ssv[,4,2] <- c(v2.ns20.2n.ssv.bigd$welch.reject, v2.ns50.2n.ssv.bigd$welch.reject, v2.ns100.2n.ssv.bigd$welch.reject)


##### Double Ns, double variances coverage rate #####

## store observed coverage rates
obs.coverage.v2.2n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v2.2n.ssv[,1,1] <- c(v2.ns20.2n.ssv.nullT$classic.coverage, v2.ns50.2n.ssv.nullT$classic.coverage, v2.ns100.2n.ssv.nullT$classic.coverage)
obs.coverage.v2.2n.ssv[,1,2] <- c(v2.ns20.2n.ssv.nullT$welch.coverage, v2.ns50.2n.ssv.nullT$welch.coverage, v2.ns100.2n.ssv.nullT$welch.coverage)

# small effect
obs.coverage.v2.2n.ssv[,2,1] <- c(v2.ns20.2n.ssv.smalld$classic.coverage, v2.ns50.2n.ssv.smalld$classic.coverage, v2.ns100.2n.ssv.smalld$classic.coverage)
obs.coverage.v2.2n.ssv[,2,2] <- c(v2.ns20.2n.ssv.smalld$welch.coverage, v2.ns50.2n.ssv.smalld$welch.coverage, v2.ns100.2n.ssv.smalld$welch.coverage)

# medium effect
obs.coverage.v2.2n.ssv[,3,1] <- c(v2.ns20.2n.ssv.midd$classic.coverage, v2.ns50.2n.ssv.midd$classic.coverage, v2.ns100.2n.ssv.midd$classic.coverage)
obs.coverage.v2.2n.ssv[,3,2] <- c(v2.ns20.2n.ssv.midd$welch.coverage, v2.ns50.2n.ssv.midd$welch.coverage, v2.ns100.2n.ssv.midd$welch.coverage)

# big effect
obs.coverage.v2.2n.ssv[,4,1] <- c(v2.ns20.2n.ssv.bigd$classic.coverage, v2.ns50.2n.ssv.bigd$classic.coverage, v2.ns100.2n.ssv.bigd$classic.coverage)
obs.coverage.v2.2n.ssv[,4,2] <- c(v2.ns20.2n.ssv.bigd$welch.coverage, v2.ns50.2n.ssv.bigd$welch.coverage, v2.ns100.2n.ssv.bigd$welch.coverage)



##### Double Ns, double variances df proportion #####

## store df proportion
df.ratio.v2.2n.ssv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v2.2n.ssv[,1,1] <- c(v2.ns20.2n.ssv.nullT$df.ratio.avg, v2.ns50.2n.ssv.nullT$df.ratio.avg, v2.ns100.2n.ssv.nullT$df.ratio.avg)  # true null
df.ratio.v2.2n.ssv[,2,1] <- c(v2.ns20.2n.ssv.smalld$df.ratio.avg, v2.ns50.2n.ssv.smalld$df.ratio.avg, v2.ns100.2n.ssv.smalld$df.ratio.avg)  # small d
df.ratio.v2.2n.ssv[,3,1] <- c(v2.ns20.2n.ssv.midd$df.ratio.avg, v2.ns50.2n.ssv.midd$df.ratio.avg, v2.ns100.2n.ssv.midd$df.ratio.avg)  # mid d
df.ratio.v2.2n.ssv[,4,1] <- c(v2.ns20.2n.ssv.bigd$df.ratio.avg, v2.ns50.2n.ssv.bigd$df.ratio.avg, v2.ns100.2n.ssv.bigd$df.ratio.avg)  # big d



##### Save tables from double Ns, double variances simulations #####
save(df.ratio.v2.2n.ssv, obs.coverage.v2.2n.ssv, reject.null.v2.2n.ssv, file='/users/joshwondra/R-projects/Welch rule/v22nSSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v22nSSVSeed2184Tables.Rdata')











##### Double Ns, 5x variances simulations #####

# true null
set.seed(2184)
v5.ns20.2n.ssv.nullT <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.nullT <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.nullT <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6), vars=c(2,10), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(2),40,sqrt(10))
cohen.diff(.2,50,sqrt(2),100,sqrt(10))
cohen.diff(.2,100,sqrt(2),200,sqrt(10))
set.seed(2184)
v5.ns20.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.28), vars=c(2,10), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),40,sqrt(10))
cohen.diff(.5,50,sqrt(2),100,sqrt(10))
cohen.diff(.5,100,sqrt(2),200,sqrt(10))
set.seed(2184)
v5.ns20.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.71), vars=c(2,10), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),40,sqrt(10))
cohen.diff(.8,50,sqrt(2),100,sqrt(10))
cohen.diff(.8,100,sqrt(2),200,sqrt(10))
set.seed(2184)
v5.ns20.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.13), vars=c(2,10), contrast=c(-1,1))

# Save simulations
save(v5.ns20.2n.ssv.nullT,v5.ns50.2n.ssv.nullT,v5.ns100.2n.ssv.nullT,v5.ns20.2n.ssv.smalld,v5.ns50.2n.ssv.smalld,v5.ns100.2n.ssv.smalld,v5.ns20.2n.ssv.midd,v5.ns50.2n.ssv.midd,v5.ns100.2n.ssv.midd, v5.ns20.2n.ssv.bigd,v5.ns50.2n.ssv.bigd,v5.ns100.2n.ssv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N2ssvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N2ssvSeed2184.Rdata')

v5n2ssvsim1 <- list(v5.ns20.2n.ssv.nullT=v5.ns20.2n.ssv.nullT$sim.results$sim1$sim.data,
                    v5.ns50.2n.ssv.nullT=v5.ns50.2n.ssv.nullT$sim.results$sim1$sim.data,
                    v5.ns100.2n.ssv.nullT=v5.ns100.2n.ssv.nullT$sim.results$sim1$sim.data,
                    v5.ns20.2n.ssv.smalld=v5.ns20.2n.ssv.smalld$sim.results$sim1$sim.data,
                    v5.ns50.2n.ssv.smalld=v5.ns50.2n.ssv.smalld$sim.results$sim1$sim.data,
                    v5.ns100.2n.ssv.smalld=v5.ns100.2n.ssv.smalld$sim.results$sim1$sim.data,
                    v5.ns20.2n.ssv.midd=v5.ns20.2n.ssv.midd$sim.results$sim1$sim.data,
                    v5.ns50.2n.ssv.midd=v5.ns50.2n.ssv.midd$sim.results$sim1$sim.data,
                    v5.ns100.2n.ssv.midd=v5.ns100.2n.ssv.midd$sim.results$sim1$sim.data,
                    v5.ns20.2n.ssv.bigd=v5.ns20.2n.ssv.bigd$sim.results$sim1$sim.data,
                    v5.ns50.2n.ssv.bigd=v5.ns50.2n.ssv.bigd$sim.results$sim1$sim.data,
                    v5.ns100.2n.ssv.bigd=v5.ns100.2n.ssv.bigd$sim.results$sim1$sim.data)
save(v5n2ssvsim1, file='/users/joshwondra/R-projects/Welch rule/v5N2ssvSeed2184-sim1.Rdata')



##### Double Ns, 5x variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v5.2n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v5.2n.ssv[,1,1] <- c(v5.ns20.2n.ssv.nullT$classic.reject, v5.ns50.2n.ssv.nullT$classic.reject, v5.ns100.2n.ssv.nullT$classic.reject)
reject.null.v5.2n.ssv[,1,2] <- c(v5.ns20.2n.ssv.nullT$welch.reject, v5.ns50.2n.ssv.nullT$welch.reject, v5.ns100.2n.ssv.nullT$welch.reject)

# small effect
reject.null.v5.2n.ssv[,2,1] <- c(v5.ns20.2n.ssv.smalld$classic.reject, v5.ns50.2n.ssv.smalld$classic.reject, v5.ns100.2n.ssv.smalld$classic.reject)
reject.null.v5.2n.ssv[,2,2] <- c(v5.ns20.2n.ssv.smalld$welch.reject, v5.ns50.2n.ssv.smalld$welch.reject, v5.ns100.2n.ssv.smalld$welch.reject)

# medium effect
reject.null.v5.2n.ssv[,3,1] <- c(v5.ns20.2n.ssv.midd$classic.reject, v5.ns50.2n.ssv.midd$classic.reject, v5.ns100.2n.ssv.midd$classic.reject)
reject.null.v5.2n.ssv[,3,2] <- c(v5.ns20.2n.ssv.midd$welch.reject, v5.ns50.2n.ssv.midd$welch.reject, v5.ns100.2n.ssv.midd$welch.reject)

# large effect
reject.null.v5.2n.ssv[,4,1] <- c(v5.ns20.2n.ssv.bigd$classic.reject, v5.ns50.2n.ssv.bigd$classic.reject, v5.ns100.2n.ssv.bigd$classic.reject)
reject.null.v5.2n.ssv[,4,2] <- c(v5.ns20.2n.ssv.bigd$welch.reject, v5.ns50.2n.ssv.bigd$welch.reject, v5.ns100.2n.ssv.bigd$welch.reject)


##### Double Ns, 5x variances coverage rate #####

## store observed coverage rates
obs.coverage.v5.2n.ssv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v5.2n.ssv[,1,1] <- c(v5.ns20.2n.ssv.nullT$classic.coverage, v5.ns50.2n.ssv.nullT$classic.coverage, v5.ns100.2n.ssv.nullT$classic.coverage)
obs.coverage.v5.2n.ssv[,1,2] <- c(v5.ns20.2n.ssv.nullT$welch.coverage, v5.ns50.2n.ssv.nullT$welch.coverage, v5.ns100.2n.ssv.nullT$welch.coverage)

# small effect
obs.coverage.v5.2n.ssv[,2,1] <- c(v5.ns20.2n.ssv.smalld$classic.coverage, v5.ns50.2n.ssv.smalld$classic.coverage, v5.ns100.2n.ssv.smalld$classic.coverage)
obs.coverage.v5.2n.ssv[,2,2] <- c(v5.ns20.2n.ssv.smalld$welch.coverage, v5.ns50.2n.ssv.smalld$welch.coverage, v5.ns100.2n.ssv.smalld$welch.coverage)

# medium effect
obs.coverage.v5.2n.ssv[,3,1] <- c(v5.ns20.2n.ssv.midd$classic.coverage, v5.ns50.2n.ssv.midd$classic.coverage, v5.ns100.2n.ssv.midd$classic.coverage)
obs.coverage.v5.2n.ssv[,3,2] <- c(v5.ns20.2n.ssv.midd$welch.coverage, v5.ns50.2n.ssv.midd$welch.coverage, v5.ns100.2n.ssv.midd$welch.coverage)

# big effect
obs.coverage.v5.2n.ssv[,4,1] <- c(v5.ns20.2n.ssv.bigd$classic.coverage, v5.ns50.2n.ssv.bigd$classic.coverage, v5.ns100.2n.ssv.bigd$classic.coverage)
obs.coverage.v5.2n.ssv[,4,2] <- c(v5.ns20.2n.ssv.bigd$welch.coverage, v5.ns50.2n.ssv.bigd$welch.coverage, v5.ns100.2n.ssv.bigd$welch.coverage)


##### Double Ns, 5x variances df proportion #####

## store df proportion
df.ratio.v5.2n.ssv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v5.2n.ssv[,1,1] <- c(v5.ns20.2n.ssv.nullT$df.ratio.avg, v5.ns50.2n.ssv.nullT$df.ratio.avg, v5.ns100.2n.ssv.nullT$df.ratio.avg)  # true null
df.ratio.v5.2n.ssv[,2,1] <- c(v5.ns20.2n.ssv.smalld$df.ratio.avg, v5.ns50.2n.ssv.smalld$df.ratio.avg, v5.ns100.2n.ssv.smalld$df.ratio.avg)  # small d
df.ratio.v5.2n.ssv[,3,1] <- c(v5.ns20.2n.ssv.midd$df.ratio.avg, v5.ns50.2n.ssv.midd$df.ratio.avg, v5.ns100.2n.ssv.midd$df.ratio.avg)  # mid d
df.ratio.v5.2n.ssv[,4,1] <- c(v5.ns20.2n.ssv.bigd$df.ratio.avg, v5.ns50.2n.ssv.bigd$df.ratio.avg, v5.ns100.2n.ssv.bigd$df.ratio.avg)  # big d




##### Save tables from double Ns, 5x variances simulations #####
save(df.ratio.v5.2n.ssv, obs.coverage.v5.2n.ssv, reject.null.v5.2n.ssv, file='/users/joshwondra/R-projects/Welch rule/v52nSSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v52nSSVSeed2184Tables.Rdata')







##### Different Ns and variances, big group small variances #####

##### 1.5x Ns, double variances simulations #####

# true null
set.seed(2184)
v2.ns20.1.5n.bsv.nullT <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6), vars=c(4,2), contrast=c(-1,1))
v2.ns50.1.5n.bsv.nullT <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6), vars=c(4,2), contrast=c(-1,1))
v2.ns100.1.5n.bsv.nullT <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6), vars=c(4,2), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(4),30,sqrt(2))
cohen.diff(.2,50,sqrt(4),75,sqrt(2))
cohen.diff(.2,100,sqrt(4),150,sqrt(2))
set.seed(2184)
v2.ns20.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.28), vars=c(4,2), contrast=c(-1,1))
v2.ns50.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.28), vars=c(4,2), contrast=c(-1,1))
v2.ns100.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.28), vars=c(4,2), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(4),30,sqrt(2))
cohen.diff(.5,50,sqrt(4),75,sqrt(2))
cohen.diff(.5,100,sqrt(4),150,sqrt(2))
set.seed(2184)
v2.ns20.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.71), vars=c(4,2), contrast=c(-1,1))
v2.ns50.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.71), vars=c(4,2), contrast=c(-1,1))
v2.ns100.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.71), vars=c(4,2), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(4),30,sqrt(2))
cohen.diff(.8,50,sqrt(4),75,sqrt(2))
cohen.diff(.8,100,sqrt(4),150,sqrt(2))
set.seed(2184)
v2.ns20.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.13), vars=c(4,2), contrast=c(-1,1))
v2.ns50.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.13), vars=c(4,2), contrast=c(-1,1))
v2.ns100.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.13), vars=c(4,2), contrast=c(-1,1))

# Save simulations
save(v2.ns20.1.5n.bsv.nullT,v2.ns50.1.5n.bsv.nullT,v2.ns100.1.5n.bsv.nullT,v2.ns20.1.5n.bsv.smalld,v2.ns50.1.5n.bsv.smalld,v2.ns100.1.5n.bsv.smalld,v2.ns20.1.5n.bsv.midd,v2.ns50.1.5n.bsv.midd,v2.ns100.1.5n.bsv.midd, v2.ns20.1.5n.bsv.bigd,v2.ns50.1.5n.bsv.bigd,v2.ns100.1.5n.bsv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N15bsvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N15bsvSeed2184.Rdata')

v2n15bsvsim1 <- list(v2.ns20.1.5n.bsv.nullT=v2.ns20.1.5n.bsv.nullT$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.bsv.nullT=v2.ns50.1.5n.bsv.nullT$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.bsv.nullT=v2.ns100.1.5n.bsv.nullT$sim.results$sim1$sim.data,
                     v2.ns20.1.5n.bsv.smalld=v2.ns20.1.5n.bsv.smalld$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.bsv.smalld=v2.ns50.1.5n.bsv.smalld$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.bsv.smalld=v2.ns100.1.5n.bsv.smalld$sim.results$sim1$sim.data,
                     v2.ns20.1.5n.bsv.midd=v2.ns20.1.5n.bsv.midd$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.bsv.midd=v2.ns50.1.5n.bsv.midd$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.bsv.midd=v2.ns100.1.5n.bsv.midd$sim.results$sim1$sim.data,
                     v2.ns20.1.5n.bsv.bigd=v2.ns20.1.5n.bsv.bigd$sim.results$sim1$sim.data,
                     v2.ns50.1.5n.bsv.bigd=v2.ns50.1.5n.bsv.bigd$sim.results$sim1$sim.data,
                     v2.ns100.1.5n.bsv.bigd=v2.ns100.1.5n.bsv.bigd$sim.results$sim1$sim.data)
save(v2n15bsvsim1, file='/users/joshwondra/R-projects/Welch rule/v2N15bsvSeed2184-sim1.Rdata')



##### 1.5x Ns, double variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v2.1.5n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v2.1.5n.bsv[,1,1] <- c(v2.ns20.1.5n.bsv.nullT$classic.reject, v2.ns50.1.5n.bsv.nullT$classic.reject, v2.ns100.1.5n.bsv.nullT$classic.reject)
reject.null.v2.1.5n.bsv[,1,2] <- c(v2.ns20.1.5n.bsv.nullT$welch.reject, v2.ns50.1.5n.bsv.nullT$welch.reject, v2.ns100.1.5n.bsv.nullT$welch.reject)

# small effect
reject.null.v2.1.5n.bsv[,2,1] <- c(v2.ns20.1.5n.bsv.smalld$classic.reject, v2.ns50.1.5n.bsv.smalld$classic.reject, v2.ns100.1.5n.bsv.smalld$classic.reject)
reject.null.v2.1.5n.bsv[,2,2] <- c(v2.ns20.1.5n.bsv.smalld$welch.reject, v2.ns50.1.5n.bsv.smalld$welch.reject, v2.ns100.1.5n.bsv.smalld$welch.reject)

# medium effect
reject.null.v2.1.5n.bsv[,3,1] <- c(v2.ns20.1.5n.bsv.midd$classic.reject, v2.ns50.1.5n.bsv.midd$classic.reject, v2.ns100.1.5n.bsv.midd$classic.reject)
reject.null.v2.1.5n.bsv[,3,2] <- c(v2.ns20.1.5n.bsv.midd$welch.reject, v2.ns50.1.5n.bsv.midd$welch.reject, v2.ns100.1.5n.bsv.midd$welch.reject)

# large effect
reject.null.v2.1.5n.bsv[,4,1] <- c(v2.ns20.1.5n.bsv.bigd$classic.reject, v2.ns50.1.5n.bsv.bigd$classic.reject, v2.ns100.1.5n.bsv.bigd$classic.reject)
reject.null.v2.1.5n.bsv[,4,2] <- c(v2.ns20.1.5n.bsv.bigd$welch.reject, v2.ns50.1.5n.bsv.bigd$welch.reject, v2.ns100.1.5n.bsv.bigd$welch.reject)


##### 1.5x Ns, double variances coverage rate #####

## store observed coverage rates
obs.coverage.v2.1.5n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v2.1.5n.bsv[,1,1] <- c(v2.ns20.1.5n.bsv.nullT$classic.coverage, v2.ns50.1.5n.bsv.nullT$classic.coverage, v2.ns100.1.5n.bsv.nullT$classic.coverage)
obs.coverage.v2.1.5n.bsv[,1,2] <- c(v2.ns20.1.5n.bsv.nullT$welch.coverage, v2.ns50.1.5n.bsv.nullT$welch.coverage, v2.ns100.1.5n.bsv.nullT$welch.coverage)

# small effect
obs.coverage.v2.1.5n.bsv[,2,1] <- c(v2.ns20.1.5n.bsv.smalld$classic.coverage, v2.ns50.1.5n.bsv.smalld$classic.coverage, v2.ns100.1.5n.bsv.smalld$classic.coverage)
obs.coverage.v2.1.5n.bsv[,2,2] <- c(v2.ns20.1.5n.bsv.smalld$welch.coverage, v2.ns50.1.5n.bsv.smalld$welch.coverage, v2.ns100.1.5n.bsv.smalld$welch.coverage)

# medium effect
obs.coverage.v2.1.5n.bsv[,3,1] <- c(v2.ns20.1.5n.bsv.midd$classic.coverage, v2.ns50.1.5n.bsv.midd$classic.coverage, v2.ns100.1.5n.bsv.midd$classic.coverage)
obs.coverage.v2.1.5n.bsv[,3,2] <- c(v2.ns20.1.5n.bsv.midd$welch.coverage, v2.ns50.1.5n.bsv.midd$welch.coverage, v2.ns100.1.5n.bsv.midd$welch.coverage)

# big effect
obs.coverage.v2.1.5n.bsv[,4,1] <- c(v2.ns20.1.5n.bsv.bigd$classic.coverage, v2.ns50.1.5n.bsv.bigd$classic.coverage, v2.ns100.1.5n.bsv.bigd$classic.coverage)
obs.coverage.v2.1.5n.bsv[,4,2] <- c(v2.ns20.1.5n.bsv.bigd$welch.coverage, v2.ns50.1.5n.bsv.bigd$welch.coverage, v2.ns100.1.5n.bsv.bigd$welch.coverage)



##### 1.5x Ns, double variances df proportion #####

## store df proportion
df.ratio.v2.1.5n.bsv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v2.1.5n.bsv[,1,1] <- c(v2.ns20.1.5n.bsv.nullT$df.ratio.avg, v2.ns50.1.5n.bsv.nullT$df.ratio.avg, v2.ns100.1.5n.bsv.nullT$df.ratio.avg)  # true null
df.ratio.v2.1.5n.bsv[,2,1] <- c(v2.ns20.1.5n.bsv.smalld$df.ratio.avg, v2.ns50.1.5n.bsv.smalld$df.ratio.avg, v2.ns100.1.5n.bsv.smalld$df.ratio.avg)  # small d
df.ratio.v2.1.5n.bsv[,3,1] <- c(v2.ns20.1.5n.bsv.midd$df.ratio.avg, v2.ns50.1.5n.bsv.midd$df.ratio.avg, v2.ns100.1.5n.bsv.midd$df.ratio.avg)  # mid d
df.ratio.v2.1.5n.bsv[,4,1] <- c(v2.ns20.1.5n.bsv.bigd$df.ratio.avg, v2.ns50.1.5n.bsv.bigd$df.ratio.avg, v2.ns100.1.5n.bsv.bigd$df.ratio.avg)  # big d




##### Save tables from 1.5x Ns, double variances simulations #####
save(df.ratio.v2.1.5n.bsv, obs.coverage.v2.1.5n.bsv, reject.null.v2.1.5n.bsv, file='/users/joshwondra/R-projects/Welch rule/v21andhalfnBSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v21andhalfnBSVSeed2184Tables.Rdata')











##### 1.5x Ns, 5x variances simulations #####

# true null
set.seed(2184)
v5.ns20.1.5n.bsv.nullT <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6), vars=c(10,2), contrast=c(-1,1))
v5.ns50.1.5n.bsv.nullT <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6), vars=c(10,2), contrast=c(-1,1))
v5.ns100.1.5n.bsv.nullT <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6), vars=c(10,2), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(10),30,sqrt(2))
cohen.diff(.2,50,sqrt(10),75,sqrt(2))
cohen.diff(.2,100,sqrt(10),150,sqrt(2))
set.seed(2184)
v5.ns20.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.28), vars=c(10,2), contrast=c(-1,1))
v5.ns50.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.28), vars=c(10,2), contrast=c(-1,1))
v5.ns100.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.28), vars=c(10,2), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(10),30,sqrt(2))
cohen.diff(.5,50,sqrt(10),75,sqrt(2))
cohen.diff(.5,100,sqrt(10),150,sqrt(2))
set.seed(2184)
v5.ns20.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.71), vars=c(10,2), contrast=c(-1,1))
v5.ns50.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.71), vars=c(10,2), contrast=c(-1,1))
v5.ns100.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.71), vars=c(10,2), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(10),30,sqrt(2))
cohen.diff(.8,50,sqrt(10),75,sqrt(2))
cohen.diff(.8,100,sqrt(10),150,sqrt(2))
set.seed(2184)
v5.ns20.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.13), vars=c(10,2), contrast=c(-1,1))
v5.ns50.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.13), vars=c(10,2), contrast=c(-1,1))
v5.ns100.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.13), vars=c(10,2), contrast=c(-1,1))

# Save simulations
# Save simulations
save(v5.ns20.1.5n.bsv.nullT,v5.ns50.1.5n.bsv.nullT,v5.ns100.1.5n.bsv.nullT,v5.ns20.1.5n.bsv.smalld,v5.ns50.1.5n.bsv.smalld,v5.ns100.1.5n.bsv.smalld,v5.ns20.1.5n.bsv.midd,v5.ns50.1.5n.bsv.midd,v5.ns100.1.5n.bsv.midd, v5.ns20.1.5n.bsv.bigd,v5.ns50.1.5n.bsv.bigd,v5.ns100.1.5n.bsv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N15bsvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N15bsvSeed2184.Rdata')

v5n15bsvsim1 <- list(v5.ns20.1.5n.bsv.nullT=v5.ns20.1.5n.bsv.nullT$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.bsv.nullT=v5.ns50.1.5n.bsv.nullT$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.bsv.nullT=v5.ns100.1.5n.bsv.nullT$sim.results$sim1$sim.data,
                     v5.ns20.1.5n.bsv.smalld=v5.ns20.1.5n.bsv.smalld$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.bsv.smalld=v5.ns50.1.5n.bsv.smalld$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.bsv.smalld=v5.ns100.1.5n.bsv.smalld$sim.results$sim1$sim.data,
                     v5.ns20.1.5n.bsv.midd=v5.ns20.1.5n.bsv.midd$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.bsv.midd=v5.ns50.1.5n.bsv.midd$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.bsv.midd=v5.ns100.1.5n.bsv.midd$sim.results$sim1$sim.data,
                     v5.ns20.1.5n.bsv.bigd=v5.ns20.1.5n.bsv.bigd$sim.results$sim1$sim.data,
                     v5.ns50.1.5n.bsv.bigd=v5.ns50.1.5n.bsv.bigd$sim.results$sim1$sim.data,
                     v5.ns100.1.5n.bsv.bigd=v5.ns100.1.5n.bsv.bigd$sim.results$sim1$sim.data)
save(v5n15bsvsim1, file='/users/joshwondra/R-projects/Welch rule/v5N15bsvSeed2184-sim1.Rdata')



##### 1.5x Ns, 5x variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v5.1.5n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v5.1.5n.bsv[,1,1] <- c(v5.ns20.1.5n.bsv.nullT$classic.reject, v5.ns50.1.5n.bsv.nullT$classic.reject, v5.ns100.1.5n.bsv.nullT$classic.reject)
reject.null.v5.1.5n.bsv[,1,2] <- c(v5.ns20.1.5n.bsv.nullT$welch.reject, v5.ns50.1.5n.bsv.nullT$welch.reject, v5.ns100.1.5n.bsv.nullT$welch.reject)

# small effect
reject.null.v5.1.5n.bsv[,2,1] <- c(v5.ns20.1.5n.bsv.smalld$classic.reject, v5.ns50.1.5n.bsv.smalld$classic.reject, v5.ns100.1.5n.bsv.smalld$classic.reject)
reject.null.v5.1.5n.bsv[,2,2] <- c(v5.ns20.1.5n.bsv.smalld$welch.reject, v5.ns50.1.5n.bsv.smalld$welch.reject, v5.ns100.1.5n.bsv.smalld$welch.reject)

# medium effect
reject.null.v5.1.5n.bsv[,3,1] <- c(v5.ns20.1.5n.bsv.midd$classic.reject, v5.ns50.1.5n.bsv.midd$classic.reject, v5.ns100.1.5n.bsv.midd$classic.reject)
reject.null.v5.1.5n.bsv[,3,2] <- c(v5.ns20.1.5n.bsv.midd$welch.reject, v5.ns50.1.5n.bsv.midd$welch.reject, v5.ns100.1.5n.bsv.midd$welch.reject)

# large effect
reject.null.v5.1.5n.bsv[,4,1] <- c(v5.ns20.1.5n.bsv.bigd$classic.reject, v5.ns50.1.5n.bsv.bigd$classic.reject, v5.ns100.1.5n.bsv.bigd$classic.reject)
reject.null.v5.1.5n.bsv[,4,2] <- c(v5.ns20.1.5n.bsv.bigd$welch.reject, v5.ns50.1.5n.bsv.bigd$welch.reject, v5.ns100.1.5n.bsv.bigd$welch.reject)


##### 1.5x Ns, 5x variances coverage rate #####

## store observed coverage rates
obs.coverage.v5.1.5n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,30", "N=50,75", "N=100,150"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v5.1.5n.bsv[,1,1] <- c(v5.ns20.1.5n.bsv.nullT$classic.coverage, v5.ns50.1.5n.bsv.nullT$classic.coverage, v5.ns100.1.5n.bsv.nullT$classic.coverage)
obs.coverage.v5.1.5n.bsv[,1,2] <- c(v5.ns20.1.5n.bsv.nullT$welch.coverage, v5.ns50.1.5n.bsv.nullT$welch.coverage, v5.ns100.1.5n.bsv.nullT$welch.coverage)

# small effect
obs.coverage.v5.1.5n.bsv[,2,1] <- c(v5.ns20.1.5n.bsv.smalld$classic.coverage, v5.ns50.1.5n.bsv.smalld$classic.coverage, v5.ns100.1.5n.bsv.smalld$classic.coverage)
obs.coverage.v5.1.5n.bsv[,2,2] <- c(v5.ns20.1.5n.bsv.smalld$welch.coverage, v5.ns50.1.5n.bsv.smalld$welch.coverage, v5.ns100.1.5n.bsv.smalld$welch.coverage)

# medium effect
obs.coverage.v5.1.5n.bsv[,3,1] <- c(v5.ns20.1.5n.bsv.midd$classic.coverage, v5.ns50.1.5n.bsv.midd$classic.coverage, v5.ns100.1.5n.bsv.midd$classic.coverage)
obs.coverage.v5.1.5n.bsv[,3,2] <- c(v5.ns20.1.5n.bsv.midd$welch.coverage, v5.ns50.1.5n.bsv.midd$welch.coverage, v5.ns100.1.5n.bsv.midd$welch.coverage)

# big effect
obs.coverage.v5.1.5n.bsv[,4,1] <- c(v5.ns20.1.5n.bsv.bigd$classic.coverage, v5.ns50.1.5n.bsv.bigd$classic.coverage, v5.ns100.1.5n.bsv.bigd$classic.coverage)
obs.coverage.v5.1.5n.bsv[,4,2] <- c(v5.ns20.1.5n.bsv.bigd$welch.coverage, v5.ns50.1.5n.bsv.bigd$welch.coverage, v5.ns100.1.5n.bsv.bigd$welch.coverage)



##### 1.5x Ns, 5x variances df proportion #####

## store df proportion
df.ratio.v5.1.5n.bsv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v5.1.5n.bsv[,1,1] <- c(v5.ns20.1.5n.bsv.nullT$df.ratio.avg, v5.ns50.1.5n.bsv.nullT$df.ratio.avg, v5.ns100.1.5n.bsv.nullT$df.ratio.avg)  # true null
df.ratio.v5.1.5n.bsv[,2,1] <- c(v5.ns20.1.5n.bsv.smalld$df.ratio.avg, v5.ns50.1.5n.bsv.smalld$df.ratio.avg, v5.ns100.1.5n.bsv.smalld$df.ratio.avg)  # small d
df.ratio.v5.1.5n.bsv[,3,1] <- c(v5.ns20.1.5n.bsv.midd$df.ratio.avg, v5.ns50.1.5n.bsv.midd$df.ratio.avg, v5.ns100.1.5n.bsv.midd$df.ratio.avg)  # mid d
df.ratio.v5.1.5n.bsv[,4,1] <- c(v5.ns20.1.5n.bsv.bigd$df.ratio.avg, v5.ns50.1.5n.bsv.bigd$df.ratio.avg, v5.ns100.1.5n.bsv.bigd$df.ratio.avg)  # big d



##### Save tables from 1.5x Ns, 5x variances simulations #####
save(df.ratio.v5.1.5n.bsv, obs.coverage.v5.1.5n.bsv, reject.null.v5.1.5n.bsv, file='/users/joshwondra/R-projects/Welch rule/v51andhalfnBSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v51andhalfnBSVSeed2184Tables.Rdata')





##### Double Ns, double variances simulations #####

# true null
set.seed(2184)
v2.ns20.2n.bsv.nullT <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6), vars=c(4,2), contrast=c(-1,1))
v2.ns50.2n.bsv.nullT <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6), vars=c(4,2), contrast=c(-1,1))
v2.ns100.2n.bsv.nullT <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6), vars=c(4,2), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(4),40,sqrt(2))
cohen.diff(.2,50,sqrt(4),100,sqrt(2))
cohen.diff(.2,100,sqrt(4),200,sqrt(2))
set.seed(2184)
v2.ns20.2n.bsv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.28), vars=c(4,2), contrast=c(-1,1))
v2.ns50.2n.bsv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.28), vars=c(4,2), contrast=c(-1,1))
v2.ns100.2n.bsv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.28), vars=c(4,2), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(4),40,sqrt(2))
cohen.diff(.5,50,sqrt(4),100,sqrt(2))
cohen.diff(.5,100,sqrt(4),200,sqrt(2))
set.seed(2184)
v2.ns20.2n.bsv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.71), vars=c(4,2), contrast=c(-1,1))
v2.ns50.2n.bsv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.71), vars=c(4,2), contrast=c(-1,1))
v2.ns100.2n.bsv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.71), vars=c(4,2), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(4),40,sqrt(2))
cohen.diff(.8,50,sqrt(4),100,sqrt(2))
cohen.diff(.8,100,sqrt(4),200,sqrt(2))
set.seed(2184)
v2.ns20.2n.bsv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.13), vars=c(4,2), contrast=c(-1,1))
v2.ns50.2n.bsv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.13), vars=c(4,2), contrast=c(-1,1))
v2.ns100.2n.bsv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.13), vars=c(4,2), contrast=c(-1,1))


# Save simulations
save(v2.ns20.2n.bsv.nullT,v2.ns50.2n.bsv.nullT,v2.ns100.2n.bsv.nullT,v2.ns20.2n.bsv.smalld,v2.ns50.2n.bsv.smalld,v2.ns100.2n.bsv.smalld,v2.ns20.2n.bsv.midd,v2.ns50.2n.bsv.midd,v2.ns100.2n.bsv.midd, v2.ns20.2n.bsv.bigd,v2.ns50.2n.bsv.bigd,v2.ns100.2n.bsv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N2bsvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v2N2bsvSeed2184.Rdata')

v2n2bsvsim1 <- list(v2.ns20.2n.bsv.nullT=v2.ns20.2n.bsv.nullT$sim.results$sim1$sim.data,
                    v2.ns50.2n.bsv.nullT=v2.ns50.2n.bsv.nullT$sim.results$sim1$sim.data,
                    v2.ns100.2n.bsv.nullT=v2.ns100.2n.bsv.nullT$sim.results$sim1$sim.data,
                    v2.ns20.2n.bsv.smalld=v2.ns20.2n.bsv.smalld$sim.results$sim1$sim.data,
                    v2.ns50.2n.bsv.smalld=v2.ns50.2n.bsv.smalld$sim.results$sim1$sim.data,
                    v2.ns100.2n.bsv.smalld=v2.ns100.2n.bsv.smalld$sim.results$sim1$sim.data,
                    v2.ns20.2n.bsv.midd=v2.ns20.2n.bsv.midd$sim.results$sim1$sim.data,
                    v2.ns50.2n.bsv.midd=v2.ns50.2n.bsv.midd$sim.results$sim1$sim.data,
                    v2.ns100.2n.bsv.midd=v2.ns100.2n.bsv.midd$sim.results$sim1$sim.data,
                    v2.ns20.2n.bsv.bigd=v2.ns20.2n.bsv.bigd$sim.results$sim1$sim.data,
                    v2.ns50.2n.bsv.bigd=v2.ns50.2n.bsv.bigd$sim.results$sim1$sim.data,
                    v2.ns100.2n.bsv.bigd=v2.ns100.2n.bsv.bigd$sim.results$sim1$sim.data)
save(v2n2bsvsim1, file='/users/joshwondra/R-projects/Welch rule/v2N2bsvSeed2184-sim1.Rdata')


##### Double Ns, double variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v2.2n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v2.2n.bsv[,1,1] <- c(v2.ns20.2n.bsv.nullT$classic.reject, v2.ns50.2n.bsv.nullT$classic.reject, v2.ns100.2n.bsv.nullT$classic.reject)
reject.null.v2.2n.bsv[,1,2] <- c(v2.ns20.2n.bsv.nullT$welch.reject, v2.ns50.2n.bsv.nullT$welch.reject, v2.ns100.2n.bsv.nullT$welch.reject)

# small effect
reject.null.v2.2n.bsv[,2,1] <- c(v2.ns20.2n.bsv.smalld$classic.reject, v2.ns50.2n.bsv.smalld$classic.reject, v2.ns100.2n.bsv.smalld$classic.reject)
reject.null.v2.2n.bsv[,2,2] <- c(v2.ns20.2n.bsv.smalld$welch.reject, v2.ns50.2n.bsv.smalld$welch.reject, v2.ns100.2n.bsv.smalld$welch.reject)

# medium effect
reject.null.v2.2n.bsv[,3,1] <- c(v2.ns20.2n.bsv.midd$classic.reject, v2.ns50.2n.bsv.midd$classic.reject, v2.ns100.2n.bsv.midd$classic.reject)
reject.null.v2.2n.bsv[,3,2] <- c(v2.ns20.2n.bsv.midd$welch.reject, v2.ns50.2n.bsv.midd$welch.reject, v2.ns100.2n.bsv.midd$welch.reject)

# large effect
reject.null.v2.2n.bsv[,4,1] <- c(v2.ns20.2n.bsv.bigd$classic.reject, v2.ns50.2n.bsv.bigd$classic.reject, v2.ns100.2n.bsv.bigd$classic.reject)
reject.null.v2.2n.bsv[,4,2] <- c(v2.ns20.2n.bsv.bigd$welch.reject, v2.ns50.2n.bsv.bigd$welch.reject, v2.ns100.2n.bsv.bigd$welch.reject)


##### Double Ns, double variances coverage rate #####

## store observed coverage rates
obs.coverage.v2.2n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v2.2n.bsv[,1,1] <- c(v2.ns20.2n.bsv.nullT$classic.coverage, v2.ns50.2n.bsv.nullT$classic.coverage, v2.ns100.2n.bsv.nullT$classic.coverage)
obs.coverage.v2.2n.bsv[,1,2] <- c(v2.ns20.2n.bsv.nullT$welch.coverage, v2.ns50.2n.bsv.nullT$welch.coverage, v2.ns100.2n.bsv.nullT$welch.coverage)

# small effect
obs.coverage.v2.2n.bsv[,2,1] <- c(v2.ns20.2n.bsv.smalld$classic.coverage, v2.ns50.2n.bsv.smalld$classic.coverage, v2.ns100.2n.bsv.smalld$classic.coverage)
obs.coverage.v2.2n.bsv[,2,2] <- c(v2.ns20.2n.bsv.smalld$welch.coverage, v2.ns50.2n.bsv.smalld$welch.coverage, v2.ns100.2n.bsv.smalld$welch.coverage)

# medium effect
obs.coverage.v2.2n.bsv[,3,1] <- c(v2.ns20.2n.bsv.midd$classic.coverage, v2.ns50.2n.bsv.midd$classic.coverage, v2.ns100.2n.bsv.midd$classic.coverage)
obs.coverage.v2.2n.bsv[,3,2] <- c(v2.ns20.2n.bsv.midd$welch.coverage, v2.ns50.2n.bsv.midd$welch.coverage, v2.ns100.2n.bsv.midd$welch.coverage)

# big effect
obs.coverage.v2.2n.bsv[,4,1] <- c(v2.ns20.2n.bsv.bigd$classic.coverage, v2.ns50.2n.bsv.bigd$classic.coverage, v2.ns100.2n.bsv.bigd$classic.coverage)
obs.coverage.v2.2n.bsv[,4,2] <- c(v2.ns20.2n.bsv.bigd$welch.coverage, v2.ns50.2n.bsv.bigd$welch.coverage, v2.ns100.2n.bsv.bigd$welch.coverage)



##### Double Ns, double variances df proportion #####

## store df proportion
df.ratio.v2.2n.bsv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v2.2n.bsv[,1,1] <- c(v2.ns20.2n.bsv.nullT$df.ratio.avg, v2.ns50.2n.bsv.nullT$df.ratio.avg, v2.ns100.2n.bsv.nullT$df.ratio.avg)  # true null
df.ratio.v2.2n.bsv[,2,1] <- c(v2.ns20.2n.bsv.smalld$df.ratio.avg, v2.ns50.2n.bsv.smalld$df.ratio.avg, v2.ns100.2n.bsv.smalld$df.ratio.avg)  # small d
df.ratio.v2.2n.bsv[,3,1] <- c(v2.ns20.2n.bsv.midd$df.ratio.avg, v2.ns50.2n.bsv.midd$df.ratio.avg, v2.ns100.2n.bsv.midd$df.ratio.avg)  # mid d
df.ratio.v2.2n.bsv[,4,1] <- c(v2.ns20.2n.bsv.bigd$df.ratio.avg, v2.ns50.2n.bsv.bigd$df.ratio.avg, v2.ns100.2n.bsv.bigd$df.ratio.avg)  # big d



##### Save tables from double Ns, double variances simulations #####
save(df.ratio.v2.2n.bsv, obs.coverage.v2.2n.bsv, reject.null.v2.2n.bsv, file='/users/joshwondra/R-projects/Welch rule/v22nBSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v22nBSVSeed2184Tables.Rdata')







##### Double Ns, 5x variances simulations #####

# true null
set.seed(2184)
v5.ns20.2n.bsv.nullT <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6), vars=c(10,2), contrast=c(-1,1))
v5.ns50.2n.bsv.nullT <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6), vars=c(10,2), contrast=c(-1,1))
v5.ns100.2n.bsv.nullT <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6), vars=c(10,2), contrast=c(-1,1))

# small d
cohen.diff(.2,20,sqrt(10),40,sqrt(2))
cohen.diff(.2,50,sqrt(10),100,sqrt(2))
cohen.diff(.2,100,sqrt(10),200,sqrt(2))
set.seed(2184)
v5.ns20.2n.bsv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.28), vars=c(10,2), contrast=c(-1,1))
v5.ns50.2n.bsv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.28), vars=c(10,2), contrast=c(-1,1))
v5.ns100.2n.bsv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.28), vars=c(10,2), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(10),40,sqrt(2))
cohen.diff(.5,50,sqrt(10),100,sqrt(2))
cohen.diff(.5,100,sqrt(10),200,sqrt(2))
set.seed(2184)
v5.ns20.2n.bsv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.71), vars=c(10,2), contrast=c(-1,1))
v5.ns50.2n.bsv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.71), vars=c(10,2), contrast=c(-1,1))
v5.ns100.2n.bsv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.71), vars=c(10,2), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(10),40,sqrt(2))
cohen.diff(.8,50,sqrt(10),100,sqrt(2))
cohen.diff(.8,100,sqrt(10),200,sqrt(2))
set.seed(2184)
v5.ns20.2n.bsv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.13), vars=c(10,2), contrast=c(-1,1))
v5.ns50.2n.bsv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.13), vars=c(10,2), contrast=c(-1,1))
v5.ns100.2n.bsv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.13), vars=c(10,2), contrast=c(-1,1))


# Save simulations
save(v5.ns20.2n.bsv.nullT,v5.ns50.2n.bsv.nullT,v5.ns100.2n.bsv.nullT,v5.ns20.2n.bsv.smalld,v5.ns50.2n.bsv.smalld,v5.ns100.2n.bsv.smalld,v5.ns20.2n.bsv.midd,v5.ns50.2n.bsv.midd,v5.ns100.2n.bsv.midd, v5.ns20.2n.bsv.bigd,v5.ns50.2n.bsv.bigd,v5.ns100.2n.bsv.bigd, file='/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N2bsvSeed2184.Rdata')
load('/users/joshwondra/Dropbox/Research/Current Projects/Welch/Data/v5N2bsvSeed2184.Rdata')

v5n2bsvsim1 <- list(v5.ns20.2n.bsv.nullT=v5.ns20.2n.bsv.nullT$sim.results$sim1$sim.data,
                    v5.ns50.2n.bsv.nullT=v5.ns50.2n.bsv.nullT$sim.results$sim1$sim.data,
                    v5.ns100.2n.bsv.nullT=v5.ns100.2n.bsv.nullT$sim.results$sim1$sim.data,
                    v5.ns20.2n.bsv.smalld=v5.ns20.2n.bsv.smalld$sim.results$sim1$sim.data,
                    v5.ns50.2n.bsv.smalld=v5.ns50.2n.bsv.smalld$sim.results$sim1$sim.data,
                    v5.ns100.2n.bsv.smalld=v5.ns100.2n.bsv.smalld$sim.results$sim1$sim.data,
                    v5.ns20.2n.bsv.midd=v5.ns20.2n.bsv.midd$sim.results$sim1$sim.data,
                    v5.ns50.2n.bsv.midd=v5.ns50.2n.bsv.midd$sim.results$sim1$sim.data,
                    v5.ns100.2n.bsv.midd=v5.ns100.2n.bsv.midd$sim.results$sim1$sim.data,
                    v5.ns20.2n.bsv.bigd=v5.ns20.2n.bsv.bigd$sim.results$sim1$sim.data,
                    v5.ns50.2n.bsv.bigd=v5.ns50.2n.bsv.bigd$sim.results$sim1$sim.data,
                    v5.ns100.2n.bsv.bigd=v5.ns100.2n.bsv.bigd$sim.results$sim1$sim.data)
save(v5n2bsvsim1, file='/users/joshwondra/R-projects/Welch rule/v5N2bsvSeed2184-sim1.Rdata')



##### Double Ns, 5x variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v5.2n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v5.2n.bsv[,1,1] <- c(v5.ns20.2n.bsv.nullT$classic.reject, v5.ns50.2n.bsv.nullT$classic.reject, v5.ns100.2n.bsv.nullT$classic.reject)
reject.null.v5.2n.bsv[,1,2] <- c(v5.ns20.2n.bsv.nullT$welch.reject, v5.ns50.2n.bsv.nullT$welch.reject, v5.ns100.2n.bsv.nullT$welch.reject)

# small effect
reject.null.v5.2n.bsv[,2,1] <- c(v5.ns20.2n.bsv.smalld$classic.reject, v5.ns50.2n.bsv.smalld$classic.reject, v5.ns100.2n.bsv.smalld$classic.reject)
reject.null.v5.2n.bsv[,2,2] <- c(v5.ns20.2n.bsv.smalld$welch.reject, v5.ns50.2n.bsv.smalld$welch.reject, v5.ns100.2n.bsv.smalld$welch.reject)

# medium effect
reject.null.v5.2n.bsv[,3,1] <- c(v5.ns20.2n.bsv.midd$classic.reject, v5.ns50.2n.bsv.midd$classic.reject, v5.ns100.2n.bsv.midd$classic.reject)
reject.null.v5.2n.bsv[,3,2] <- c(v5.ns20.2n.bsv.midd$welch.reject, v5.ns50.2n.bsv.midd$welch.reject, v5.ns100.2n.bsv.midd$welch.reject)

# large effect
reject.null.v5.2n.bsv[,4,1] <- c(v5.ns20.2n.bsv.bigd$classic.reject, v5.ns50.2n.bsv.bigd$classic.reject, v5.ns100.2n.bsv.bigd$classic.reject)
reject.null.v5.2n.bsv[,4,2] <- c(v5.ns20.2n.bsv.bigd$welch.reject, v5.ns50.2n.bsv.bigd$welch.reject, v5.ns100.2n.bsv.bigd$welch.reject)


##### Double Ns, 5x variances coverage rate #####

## store observed coverage rates
obs.coverage.v5.2n.bsv <- array(dim=c(3,4,2), dimnames=list(c("N=20,40", "N=50,100", "N=100,200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v5.2n.bsv[,1,1] <- c(v5.ns20.2n.bsv.nullT$classic.coverage, v5.ns50.2n.bsv.nullT$classic.coverage, v5.ns100.2n.bsv.nullT$classic.coverage)
obs.coverage.v5.2n.bsv[,1,2] <- c(v5.ns20.2n.bsv.nullT$welch.coverage, v5.ns50.2n.bsv.nullT$welch.coverage, v5.ns100.2n.bsv.nullT$welch.coverage)

# small effect
obs.coverage.v5.2n.bsv[,2,1] <- c(v5.ns20.2n.bsv.smalld$classic.coverage, v5.ns50.2n.bsv.smalld$classic.coverage, v5.ns100.2n.bsv.smalld$classic.coverage)
obs.coverage.v5.2n.bsv[,2,2] <- c(v5.ns20.2n.bsv.smalld$welch.coverage, v5.ns50.2n.bsv.smalld$welch.coverage, v5.ns100.2n.bsv.smalld$welch.coverage)

# medium effect
obs.coverage.v5.2n.bsv[,3,1] <- c(v5.ns20.2n.bsv.midd$classic.coverage, v5.ns50.2n.bsv.midd$classic.coverage, v5.ns100.2n.bsv.midd$classic.coverage)
obs.coverage.v5.2n.bsv[,3,2] <- c(v5.ns20.2n.bsv.midd$welch.coverage, v5.ns50.2n.bsv.midd$welch.coverage, v5.ns100.2n.bsv.midd$welch.coverage)

# big effect
obs.coverage.v5.2n.bsv[,4,1] <- c(v5.ns20.2n.bsv.bigd$classic.coverage, v5.ns50.2n.bsv.bigd$classic.coverage, v5.ns100.2n.bsv.bigd$classic.coverage)
obs.coverage.v5.2n.bsv[,4,2] <- c(v5.ns20.2n.bsv.bigd$welch.coverage, v5.ns50.2n.bsv.bigd$welch.coverage, v5.ns100.2n.bsv.bigd$welch.coverage)



##### Double Ns, 5x variances df proportion #####

## store df proportion
df.ratio.v5.2n.bsv <- array(dim=c(3,4,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))

# true null
df.ratio.v5.2n.bsv[,1,1] <- c(v5.ns20.2n.bsv.nullT$df.ratio.avg, v5.ns50.2n.bsv.nullT$df.ratio.avg, v5.ns100.2n.bsv.nullT$df.ratio.avg)  # true null
df.ratio.v5.2n.bsv[,2,1] <- c(v5.ns20.2n.bsv.smalld$df.ratio.avg, v5.ns50.2n.bsv.smalld$df.ratio.avg, v5.ns100.2n.bsv.smalld$df.ratio.avg)  # small d
df.ratio.v5.2n.bsv[,3,1] <- c(v5.ns20.2n.bsv.midd$df.ratio.avg, v5.ns50.2n.bsv.midd$df.ratio.avg, v5.ns100.2n.bsv.midd$df.ratio.avg)  # mid d
df.ratio.v5.2n.bsv[,4,1] <- c(v5.ns20.2n.bsv.bigd$df.ratio.avg, v5.ns50.2n.bsv.bigd$df.ratio.avg, v5.ns100.2n.bsv.bigd$df.ratio.avg)  # big d



##### Save tables from double Ns, 5x variances simulations #####
save(df.ratio.v5.2n.bsv, obs.coverage.v5.2n.bsv, reject.null.v5.2n.bsv, file='/users/joshwondra/R-projects/Welch rule/v52nBSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v52nBSVSeed2184Tables.Rdata')







##### False Positives and Power #####

round(reject.null.ve.ne, digits=2)
round(reject.null.v2.ne, digits=2)
round(reject.null.v5.ne, digits=2)

round(reject.null.ve.ne, digits=2)
round(reject.null.ve.1.5n, digits=2)
round(reject.null.ve.2n, digits=2)

round(reject.null.ve.1.5n, digits=2)
round(reject.null.v2.1.5n.ssv, digits=2)
round(reject.null.v5.1.5n.ssv, digits=2)

round(reject.null.v2.1.5n.bsv, digits=2)

round(reject.null.ve.2n, digits=2)
round(reject.null.v2.2n.ssv, digits=2)
round(reject.null.v5.2n.ssv, digits=2)

## differences between expectation and observation
round(reject.null.ve.ne-exp.rejects.ne, digits=2)
round(reject.null.v2.ne-exp.rejects.ne, digits=2)
round(reject.null.v5.ne-exp.rejects.ne, digits=2)

round(reject.null.ve.ne-exp.rejects.ne, digits=2)
round(reject.null.ve.1.5n-exp.rejects.1.5n, digits=2)
round(reject.null.ve.2n-exp.rejects.2n, digits=2)

round(reject.null.ve.1.5n-exp.rejects.1.5n, digits=2)
round(reject.null.v2.1.5n.ssv-exp.rejects.1.5n, digits=2)
round(reject.null.v5.1.5n.ssv-exp.rejects.1.5n, digits=2)

round(reject.null.ve.2n-exp.rejects.2n, digits=2)
round(reject.null.v2.2n.ssv-exp.rejects.2n, digits=2)
round(reject.null.v5.2n.ssv-exp.rejects.2n, digits=2)

##### Expected False Positive Rate and Power #####

## Equal Ns
exp.rejects.ne <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))
exp.rejects.ne[,1,] <- rep(.05,3)
exp.rejects.ne[,2,] <- c(.095, .168, .291)
exp.rejects.ne[,3,] <- c(.338, .697, .940)
exp.rejects.ne[,4,] <- c(.693, .977, 1)

## 1.5x Ns
exp.rejects.1.5n <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))
exp.rejects.1.5n[,1,] <- rep(.05,3)
exp.rejects.1.5n[,2,] <- c(.104, .192, .339)
exp.rejects.1.5n[,3,] <- c(.397,.776,.971)
exp.rejects.1.5n[,4,] <- c(.775, .992, 1)

## Double Ns
exp.rejects.2n <- array(dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8")))
exp.rejects.2n[,1,] <- rep(.05,3)
exp.rejects.2n[,2,] <- c(.111,.209,.370)
exp.rejects.2n[,3,] <- c(.435, .818, .983)
exp.rejects.2n[,4,] <- c(.819,.996,1)



##### Coverage Rates #####

round(obs.coverage.ve.ne, digits=2)
round(obs.coverage.v2.ne, digits=2)
round(obs.coverage.v5.ne, digits=2)

round(obs.coverage.ve.ne, digits=2)
round(obs.coverage.ve.1.5n, digits=2)
round(obs.coverage.ve.2n, digits=2)

round(obs.coverage.ve.1.5n, digits=2)
round(obs.coverage.v2.1.5n.ssv, digits=2)
round(obs.coverage.v5.1.5n.ssv, digits=2)

round(obs.coverage.ve.2n, digits=2)
round(obs.coverage.v2.2n.ssv, digits=2)
round(obs.coverage.v5.2n.ssv, digits=2)


#### Expected coverage rates #####

exp.coverage <- array(rep(.95, 24), dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c('classic', 'welch')))




#### Df Ratio tables ####

## store selected dfs

# order of sims is ve.ns20.ne.nullT, ve.ns50.ne.nullT, ve.ns100.ne.nullT; v2.ns50.ne.nullT, v5.ns50.ne.nullT; ve.ns50.1.5n.nullT, ve.ns50.2n.nullT; v2.ns50.1.5n.ssv.nullT, v2.ns50.2n.ssv.nullT, v5.ns50.1.5n.ssv.nullT, v5.ns50.2n.ssv.nullT; v2.ns50.1.5n.bsv.nullT, v2.ns50.2n.bsv.nullT, v5.ns50.1.5n.bsv.nullT, v5.ns50.2n.bsv.nullT
df.var.ratio <- factor(c(1,1,1,2,5,1,1,2,2,5,5,2,2,5,5))
df.ns <- factor(c('20,20','50,50','100,100','50,50','50,50','50,75','50,100','50,75','50,100','50,75','50,100','75,50','100,50','75,50','100,50'))
dfs.list <- list(
    ve.ns20.ne.nullT=unlist(lapply(ve.ns20.ne.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    ve.ns50.ne.nullT=unlist(lapply(ve.ns50.ne.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    ve.ns100.ne.nullT=unlist(lapply(ve.ns100.ne.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v2.ns50.ne.nullT=unlist(lapply(v2.ns50.ne.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v5.ns50.ne.nullT=unlist(lapply(v5.ns50.ne.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    ve.ns50.1.5n.nullT=unlist(lapply(ve.ns50.1.5n.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    ve.ns50.2n.nullT=unlist(lapply(ve.ns50.2n.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v2.ns50.1.5n.ssv.nullT=unlist(lapply(v2.ns50.1.5n.ssv.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v2.ns50.2n.ssv.nullT=unlist(lapply(v2.ns50.2n.ssv.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v5.ns50.1.5n.ssv.nullT=unlist(lapply(v5.ns50.1.5n.ssv.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v5.ns50.2n.ssv.nullT=unlist(lapply(v5.ns50.2n.ssv.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v2.ns50.1.5n.bsv.nullT=unlist(lapply(v2.ns50.1.5n.bsv.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v2.ns50.2n.bsv.nullT=unlist(lapply(v2.ns50.2n.bsv.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v5.ns50.1.5n.bsv.nullT=unlist(lapply(v5.ns50.1.5n.bsv.nullT$sim.results, function(x){x$df.welch/x$df.classic})), 
    v5.ns50.2n.bsv.nullT=unlist(lapply(v5.ns50.2n.bsv.nullT$sim.results, function(x){x$df.welch/x$df.classic})))

save(dfs.list, file='/users/joshwondra/R-projects/Welch rule/dfs.list.RData')


df.ratio.ve.ne
df.ratio.ve.1.5n
df.ratio.ve.2n

df.ratio.v2.ne
df.ratio.v2.1.5n.ssv
df.ratio.v2.1.5n.bsv
df.ratio.v2.2n.ssv
df.ratio.v2.2n.bsv

df.ratio.v5.ne
df.ratio.v5.1.5n.ssv
df.ratio.v5.1.5n.bsv
df.ratio.v5.2n.ssv
df.ratio.v5.2n.bsv


##### DISPLAY TABLES #####

round(reject.null.ve.ne, digits=2)
round(obs.coverage.ve.ne, digits=2)
round(df.ratio.ve.ne, digits=2)
round(exp.rejects.ne, digits=2)
round(reject.null.ve.ne-array(c(exp.rejects.ne, exp.rejects.ne), dim=c(3,4,2)), digits=2)

round(reject.null.v2.ne, digits=2)
round(obs.coverage.v2.ne, digits=2)
round(df.ratio.v2.ne, digits=2)
round(exp.rejects.ne, digits=2)
round(reject.null.v2.ne-array(c(exp.rejects.ne, exp.rejects.ne), dim=c(3,4,2)), digits=2)


round(reject.null.v5.ne, digits=2)
round(obs.coverage.v5.ne, digits=2)
round(df.ratio.v5.ne, digits=2)
round(exp.rejects.ne, digits=2)
round(reject.null.v5.ne-array(c(exp.rejects.ne, exp.rejects.ne), dim=c(3,4,2)), digits=2)

round(reject.null.ve.1.5n, digits=2)
round(obs.coverage.ve.1.5n, digits=2)
round(df.ratio.ve.1.5n, digits=2)
round(exp.rejects.1.5n, digits=2)
round(reject.null.ve.1.5n-array(c(exp.rejects.1.5n, exp.rejects.1.5n), dim=c(3,4,2)), digits=2)


round(reject.null.ve.2n, digits=2)
round(obs.coverage.ve.2n, digits=2)
round(df.ratio.ve.2n, digits=2)
round(exp.rejects.2n, digits=2)
round(reject.null.ve.2n-array(c(exp.rejects.2n, exp.rejects.2n), dim=c(3,4,2)), digits=2)


round(reject.null.v2.1.5n.ssv, digits=2)
round(obs.coverage.v2.1.5n.ssv, digits=2)
round(df.ratio.v2.1.5n.ssv, digits=2)
round(exp.rejects.1.5n, digits=2)
round(reject.null.v2.1.5n.ssv-array(c(exp.rejects.1.5n, exp.rejects.1.5n), dim=c(3,4,2)), digits=2)


round(reject.null.v5.1.5n.ssv, digits=2)
round(obs.coverage.v5.1.5n.ssv, digits=2)
round(df.ratio.v5.1.5n.ssv, digits=2)
round(exp.rejects.1.5n, digits=2)
round(reject.null.v5.1.5n.ssv-array(c(exp.rejects.1.5n, exp.rejects.1.5n), dim=c(3,4,2)), digits=2)


round(reject.null.v2.2n.ssv, digits=2)
round(obs.coverage.v2.2n.ssv, digits=2)
round(df.ratio.v2.2n.ssv, digits=2)
round(exp.rejects.2n, digits=2)
round(reject.null.v2.2n.ssv-array(c(exp.rejects.1.5n, exp.rejects.1.5n), dim=c(3,4,2)), digits=2)


round(reject.null.v5.2n.ssv, digits=2)
round(obs.coverage.v5.2n.ssv, digits=2)
round(df.ratio.v5.2n.ssv, digits=2)
round(exp.rejects.2n, digits=2)

round(reject.null.v2.1.5n.bsv, digits=2)
round(obs.coverage.v2.1.5n.bsv, digits=2)
round(df.ratio.v2.1.5n.bsv, digits=2)
round(exp.rejects.1.5n, digits=2)

round(reject.null.v5.1.5n.bsv, digits=2)
round(obs.coverage.v5.1.5n.bsv, digits=2)
round(df.ratio.v5.1.5n.bsv, digits=2)
round(exp.rejects.1.5n, digits=2)

round(reject.null.v2.2n.bsv, digits=2)
round(obs.coverage.v2.2n.bsv, digits=2)
round(df.ratio.v2.2n.bsv, digits=2)
round(exp.rejects.2n, digits=2)

round(reject.null.v5.2n.bsv, digits=2)
round(obs.coverage.v5.2n.bsv, digits=2)
round(df.ratio.v5.2n.bsv, digits=2)
round(exp.rejects.2n, digits=2)





##### DF test #####

set.seed(2184)
ve.ns20.ne.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,2.2), contrast=c(-1,1))
ve.ns50.ne.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,2.2), contrast=c(-1,1))
ve.ns100.ne.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,2.2), contrast=c(-1,1))

# expected df ratio for v2 n1.5, small n 20, big n big variance
num <- (4/30+2/20)^2
denom <- ((4/30)^2/(30-1))+((2/20)^2/(20-1))
num/denom

welch.df <- function(N1, var1, N2, var2){
    num <- (var1/N1+var2/N2)^2
    denom <- ((var1/N1)^2/(N1-1))+((var2/N2)^2/(N2-1))
    result <- num/denom
    return(result)
}

50+75-2 #clssic df
welch.df(N1=50, var1=2, N2=75, var2=2) #105.16
welch.df(N1=50, var1=2, N2=75, var2=4) #122.53
welch.df(N1=50, var1=2, N2=75, var2=10) #110.09
welch.df(N1=50, var1=4, N2=75, var2=2) #81.14
welch.df(N1=50, var1=10, N2=75, var2=2) #62.21




args(t.compare)

t.compare(1,c(20,30),c(0,0),c(4,2),c(-1,1))
nsims<-5
Ns<-c(20,30)
means <- c(0,0)
vars <- c(4,2)
contrast <- c(-1,1)





##### 2 x 2 Interaction #####

##### Equal Ns, equal variances simulations #####

set.seed(2184)
int.ve.ns50.ne.nullT <- t.compare(nsims=10000, Ns=c(50,50,50,50), means=c(6,6,6,6), vars=c(2,2,2,2), contrast=c(-1,1,1,-1))

set.seed(2184)
int.ve.ns50.ne.midd <- t.compare(nsims=10000, Ns=c(50,50,50,50), means=c(6,6.71,6.71,6), vars=c(2,2,2,2), contrast=c(-1,1,1,-1))

##### Equal Ns, equal variances false positives and power #####

## array to store proportion of rejected null hypotheses
int.reject.null.ve.ne <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.reject.null.ve.ne[1,1,1] <- int.ve.ns50.ne.nullT$classic.reject
int.reject.null.ve.ne[1,1,2] <- int.ve.ns50.ne.nullT$welch.reject
int.reject.null.ve.ne[1,2,1] <- int.ve.ns50.ne.midd$classic.reject
int.reject.null.ve.ne[1,2,2] <- int.ve.ns50.ne.midd$welch.reject


##### Equal Ns, equal variances coverage rate #####

## store observed coverage rates
int.obs.coverage.ve.ne <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.obs.coverage.ve.ne[1,1,1] <- int.ve.ns50.ne.nullT$classic.coverage
int.obs.coverage.ve.ne[1,1,2] <- int.ve.ns50.ne.nullT$welch.coverage
int.obs.coverage.ve.ne[1,2,1] <- int.ve.ns50.ne.midd$classic.coverage
int.obs.coverage.ve.ne[1,2,2] <- int.ve.ns50.ne.midd$welch.coverage



##### Save tables from equal Ns, equal variances simulations #####
save(int.reject.null.ve.ne, int.obs.coverage.ve.ne, file='/users/joshwondra/R-projects/Welch rule/interaction2x2-veNeSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/interaction2x2-veNeSeed2184Tables.Rdata')






##### Big group/small variance simulations #####

set.seed(2184)
int.v5.ns50.n2.nullT.bignsmallvar <- t.compare(nsims=10000, Ns=c(100,50,50,50), means=c(6,6,6,6), vars=c(2,10,10,10), contrast=c(-1,1,1,-1))

set.seed(2184)
int.v5.ns50.n2.midd.bignsmallvar <- t.compare(nsims=10000, Ns=c(100,50,50,50), means=c(6,6.71,6.71,6), vars=c(2,10,10,10), contrast=c(-1,1,1,-1))

##### Big group/small variance false positives and power #####

## array to store proportion of rejected null hypotheses
int.reject.null.v5.n2.bignsmallvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.reject.null.v5.n2.bignsmallvar[1,1,1] <- int.v5.ns50.n2.nullT.bignsmallvar$classic.reject
int.reject.null.v5.n2.bignsmallvar[1,1,2] <- int.v5.ns50.n2.nullT.bignsmallvar$welch.reject
int.reject.null.v5.n2.bignsmallvar[1,2,1] <- int.v5.ns50.n2.midd.bignsmallvar$classic.reject
int.reject.null.v5.n2.bignsmallvar[1,2,2] <- int.v5.ns50.n2.midd.bignsmallvar$welch.reject


##### Big group/small variance coverage rate #####

## store observed coverage rates
int.obs.coverage.v5.n2.bignsmallvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.obs.coverage.v5.n2.bignsmallvar[1,1,1] <- int.v5.ns50.n2.nullT.bignsmallvar$classic.coverage
int.obs.coverage.v5.n2.bignsmallvar[1,1,2] <- int.v5.ns50.n2.nullT.bignsmallvar$welch.coverage
int.obs.coverage.v5.n2.bignsmallvar[1,2,1] <- int.v5.ns50.n2.midd.bignsmallvar$classic.coverage
int.obs.coverage.v5.n2.bignsmallvar[1,2,2] <- int.v5.ns50.n2.midd.bignsmallvar$welch.coverage


##### Save tables from big group/small variance simulations #####
save(int.reject.null.v5.n2.bignsmallvar, int.obs.coverage.v5.n2.bignsmallvar, file='/users/joshwondra/R-projects/Welch rule/interaction2x2-biggroupsmallvarSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/interaction2x2-biggroupsmallvarSeed2184Tables.Rdata')









##### Big group/large variance simulations #####

set.seed(2184)
int.v5.ns50.n2.nullT.bignbigvar <- t.compare(nsims=10000, Ns=c(100,50,50,50), means=c(6,6,6,6), vars=c(10,2,2,2), contrast=c(-1,1,1,-1))

set.seed(2184)
int.v5.ns50.n2.midd.bignbigvar <- t.compare(nsims=10000, Ns=c(100,50,50,50), means=c(6,6.71,6.71,6), vars=c(10,2,2,2), contrast=c(-1,1,1,-1))

##### Big group/large variance false positives and power #####

## array to store proportion of rejected null hypotheses
int.reject.null.v5.n2.bignbigvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.reject.null.v5.n2.bignbigvar[1,1,1] <- int.v5.ns50.n2.nullT.bignbigvar$classic.reject
int.reject.null.v5.n2.bignbigvar[1,1,2] <- int.v5.ns50.n2.nullT.bignbigvar$welch.reject
int.reject.null.v5.n2.bignbigvar[1,2,1] <- int.v5.ns50.n2.midd.bignbigvar$classic.reject
int.reject.null.v5.n2.bignbigvar[1,2,2] <- int.v5.ns50.n2.midd.bignbigvar$welch.reject


##### Big group/large variance coverage rate #####

## store observed coverage rates
int.obs.coverage.v5.n2.bignbigvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.obs.coverage.v5.n2.bignbigvar[1,1,1] <- int.v5.ns50.n2.nullT.bignbigvar$classic.coverage
int.obs.coverage.v5.n2.bignbigvar[1,1,2] <- int.v5.ns50.n2.nullT.bignbigvar$welch.coverage
int.obs.coverage.v5.n2.bignbigvar[1,2,1] <- int.v5.ns50.n2.midd.bignbigvar$classic.coverage
int.obs.coverage.v5.n2.bignbigvar[1,2,2] <- int.v5.ns50.n2.midd.bignbigvar$welch.coverage



##### Save tables from big group/big variance simulations #####
save(int.reject.null.v5.n2.bignbigvar, int.obs.coverage.v5.n2.bignbigvar, file='/users/joshwondra/R-projects/Welch rule/interaction2x2-biggroupbigvarSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/interaction2x2-biggroupbigvarSeed2184Tables.Rdata')









##### Small group/small variance simulations #####

set.seed(2184)
int.v5.ns50.n2.nullT.smallnsmallvar <- t.compare(nsims=10000, Ns=c(50,100,100,100), means=c(6,6,6,6), vars=c(2,10,10,10), contrast=c(-1,1,1,-1))

set.seed(2184)
int.v5.ns50.n2.midd.smallnsmallvar <- t.compare(nsims=10000, Ns=c(50,100,100,100), means=c(6,6.71,6.71,6), vars=c(2,10,10,10), contrast=c(-1,1,1,-1))

##### Small group/small variance false positives and power #####

## array to store proportion of rejected null hypotheses
int.reject.null.v5.n2.smallnsmallvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.reject.null.v5.n2.smallnsmallvar[1,1,1] <- int.v5.ns50.n2.nullT.smallnsmallvar$classic.reject
int.reject.null.v5.n2.smallnsmallvar[1,1,2] <- int.v5.ns50.n2.nullT.smallnsmallvar$welch.reject
int.reject.null.v5.n2.smallnsmallvar[1,2,1] <- int.v5.ns50.n2.midd.smallnsmallvar$classic.reject
int.reject.null.v5.n2.smallnsmallvar[1,2,2] <- int.v5.ns50.n2.midd.smallnsmallvar$welch.reject


##### Small group/small variance coverage rate #####

## store observed coverage rates
int.obs.coverage.v5.n2.smallnsmallvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.obs.coverage.v5.n2.smallnsmallvar[1,1,1] <- int.v5.ns50.n2.nullT.smallnsmallvar$classic.coverage
int.obs.coverage.v5.n2.smallnsmallvar[1,1,2] <- int.v5.ns50.n2.nullT.smallnsmallvar$welch.coverage
int.obs.coverage.v5.n2.smallnsmallvar[1,2,1] <- int.v5.ns50.n2.midd.smallnsmallvar$classic.coverage
int.obs.coverage.v5.n2.smallnsmallvar[1,2,2] <- int.v5.ns50.n2.midd.smallnsmallvar$welch.coverage


##### Save tables from small group/small variance simulations #####
save(int.reject.null.v5.n2.smallnsmallvar, int.obs.coverage.v5.n2.smallnsmallvar, file='/users/joshwondra/R-projects/Welch rule/interaction2x2-smallgroupsmallvarSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/interaction2x2-smallgroupsmallvarSeed2184Tables.Rdata')















##### Small group/large variance simulations #####

set.seed(2184)
int.v5.ns50.n2.nullT.smallnbigvar <- t.compare(nsims=10000, Ns=c(50,100,100,100), means=c(6,6,6,6), vars=c(10,2,2,2), contrast=c(-1,1,1,-1))

set.seed(2184)
int.v5.ns50.n2.midd.smallnbigvar <- t.compare(nsims=10000, Ns=c(50,100,100,100), means=c(6,6.71,6.71,6), vars=c(10,2,2,2), contrast=c(-1,1,1,-1))

##### Small group/large variance false positives and power #####

## array to store proportion of rejected null hypotheses
int.reject.null.v5.n2.smallnbigvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.reject.null.v5.n2.smallnbigvar[1,1,1] <- int.v5.ns50.n2.nullT.smallnbigvar$classic.reject
int.reject.null.v5.n2.smallnbigvar[1,1,2] <- int.v5.ns50.n2.nullT.smallnbigvar$welch.reject
int.reject.null.v5.n2.smallnbigvar[1,2,1] <- int.v5.ns50.n2.midd.smallnbigvar$classic.reject
int.reject.null.v5.n2.smallnbigvar[1,2,2] <- int.v5.ns50.n2.midd.smallnbigvar$welch.reject


##### Small group/large variance coverage rate #####

## store observed coverage rates
int.obs.coverage.v5.n2.smallnbigvar <- array(dim=c(1,2,2), dimnames=list(c("N=50"), c("d=0", "d=.5"), c("classic", "welch")))

# true null
int.obs.coverage.v5.n2.smallnbigvar[1,1,1] <- int.v5.ns50.n2.nullT.smallnbigvar$classic.coverage
int.obs.coverage.v5.n2.smallnbigvar[1,1,2] <- int.v5.ns50.n2.nullT.smallnbigvar$welch.coverage
int.obs.coverage.v5.n2.smallnbigvar[1,2,1] <- int.v5.ns50.n2.midd.smallnbigvar$classic.coverage
int.obs.coverage.v5.n2.smallnbigvar[1,2,2] <- int.v5.ns50.n2.midd.smallnbigvar$welch.coverage



##### Save tables from small group/big variance simulations #####
save(int.reject.null.v5.n2.smallnbigvar, int.obs.coverage.v5.n2.smallnbigvar, file='/users/joshwondra/R-projects/Welch rule/interaction2x2-smallgroupbigvarSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/interaction2x2-smallgroupbigvarSeed2184Tables.Rdata')

