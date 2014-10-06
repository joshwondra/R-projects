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
    
    result <- list(ihat=ihat, se.classic=se.classic, t.classic=t.classic, df.classic=df.classic, p.classic=p.classic, se.welch=se.welch, t.welch=t.welch, df.welch=df.welch, p.welch=p.welch)
    return(result)
}

## function to run simulations
# contrast is an n x m matrix with each of the n rows as a contrast and each of the m columns representing a group
t.compare <- function(nsims, Ns, means, vars, contrast) {
    sims <- vector('list',nsims)
    group <- rep(1:length(Ns), Ns)   #vector of length N with group codes
    #dv <- vector('numeric',sum(Ns))
    
    sims <- lapply(sims, function(x){
        dv <- rnorm(n=sum(Ns), mean=rep(means,times=Ns), sd=sqrt(rep(vars,times=Ns)))
        fit <- t.contrast(dv,group,contrast)
        ihat <- fit$ihat
        
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
        
        current.sim <- list(contrast, ihat, se.classic, t.classic, df.classic, p.classic, se.welch, t.welch, df.welch, p.welch) # add matrix(c(dv,group), ncol=2, dimnames=list(c(),c('dv','group'))) to save the data
        names(current.sim) <- c('contrast', 'ihat', 'se.classic', 't.classic', 'df.classic', 'p.classic', 'se.welch', 't.welch', 'df.welch', 'p.welch') # add data if saving the data
        return(current.sim)
    })
    
    names(sims)[1:length(sims)] <- paste('sim',1:length(sims),sep='')
    
    # store proportion of rejected null hypotheses
    classic.reject <- sum(lapply(sims, '[[', 'p.classic')[5:length(sims)]<=.05)/length(5:length(sims))
    welch.reject <- sum(lapply(sims, '[[', 'p.welch')[5:length(sims)]<=.05)/length(5:length(sims))
    
    #store df ratio
    df.classic.vector <- unlist(lapply(sims, '[[', 'df.classic'))
    df.welch.vector <- unlist(lapply(sims, '[[', 'df.welch'))
    df.ratio <- df.welch.vector/df.classic.vector
    df.ratio.avg <- mean(df.ratio)
    
    # compute coverage rate
    true.ihat <- contrast %*% means
    
    obs.ihat <- data.frame(lapply(sims, '[[', 'ihat'))
    classic.df <- data.frame(lapply(sims, '[[', 'df.classic'))
    welch.df <- data.frame(lapply(sims, '[[', 'df.welch'))
    classic.ses <- data.frame(lapply(sims, '[[', 'se.classic'))
    welch.ses <- data.frame(lapply(sims, '[[', 'se.welch'))
    
    t.classic <- apply(classic.df, 1, function(x){qt(.025, df=x)})        
    classic.lb <- obs.ihat-t.classic*classic.ses
    classic.ub <- obs.ihat+t.classic*classic.ses
    classic.coverage.logical <- (classic.ub-true.ihat)*(true.ihat-classic.lb)>0
    
    t.welch <- apply(welch.df, 1, function(x){qt(.025, df=x)})        
    welch.lb <- obs.ihat-t.welch*welch.ses
    welch.ub <- obs.ihat+t.welch*welch.ses
    welch.coverage.logical <- (welch.ub-true.ihat)*(true.ihat-welch.lb)>0       
    
    classic.coverage <- prop.table(table(classic.coverage.logical))[[2]]
    welch.coverage <- prop.table(table(welch.coverage.logical))[[2]]
    
    # return data 
    return(list=c(classic.reject=classic.reject, welch.reject=welch.reject, df.ratio.avg=df.ratio.avg, classic.coverage=classic.coverage, welch.coverage=welch.coverage, df.ratio, sims))
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
ve.ns20.ne.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ve.ns50.ne.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ve.ns100.ne.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(2))
set.seed(2184)
ve.ns20.ne.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ve.ns50.ne.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ve.ns100.ne.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(2))
set.seed(2184)
ve.ns20.ne.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ve.ns50.ne.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ve.ns100.ne.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))


## Save simulations
save(ve.ns20.ne.nullT,ve.ns50.ne.nullT,ve.ns100.ne.nullT,ve.ns20.ne.smalld,ve.ns50.ne.smalld,ve.ns100.ne.smalld,ve.ns20.ne.midd,ve.ns50.ne.midd,ve.ns100.ne.midd,ve.ns20.ne.bigd,ve.ns50.ne.bigd,ve.ns100.ne.bigd, file='/users/joshwondra/R-projects/Welch rule/veNeSeed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/veNeSeed2184.Rdata')


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


##### Save tables from equal Ns, equal variances simulations #####
save(obs.coverage.ve.ne, reject.null.ve.ne, file='/users/joshwondra/R-projects/Welch rule/veNeSeed2184Tables.Rdata')
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
v2.ns20.ne.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.346), vars=c(2,4), contrast=c(-1,1))
v2.ns50.ne.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.346), vars=c(2,4), contrast=c(-1,1))
v2.ns100.ne.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.346), vars=c(2,4), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(4))
set.seed(2184)
v2.ns20.ne.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.866), vars=c(2,4), contrast=c(-1,1))
v2.ns50.ne.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.866), vars=c(2,4), contrast=c(-1,1))
v2.ns100.ne.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.866), vars=c(2,4), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(4))
set.seed(2184)
v2.ns20.ne.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.386), vars=c(2,4), contrast=c(-1,1))
v2.ns50.ne.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.386), vars=c(2,4), contrast=c(-1,1))
v2.ns100.ne.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.386), vars=c(2,4), contrast=c(-1,1))


## Save simulations
save(v2.ns20.ne.nullT,v2.ns50.ne.nullT,v2.ns100.ne.nullT,v2.ns20.ne.smalld,v2.ns50.ne.smalld,v2.ns100.ne.smalld,v2.ns20.ne.midd,v2.ns50.ne.midd,v2.ns100.ne.midd, v2.ns20.ne.bigd,v2.ns50.ne.bigd,v2.ns100.ne.bigd, file='/users/joshwondra/R-projects/Welch rule/v2NeSeed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v2NeSeed2184.Rdata')


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


##### Save tables from equal Ns, double variances simulations #####
save(obs.coverage.v2.ne, reject.null.v2.ne, file='/users/joshwondra/R-projects/Welch rule/v2NeSeed2184Tables.Rdata')
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
v5.ns20.ne.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.490), vars=c(2,10), contrast=c(-1,1))
v5.ns50.ne.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.490), vars=c(2,10), contrast=c(-1,1))
v5.ns100.ne.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.490), vars=c(2,10), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(10))
set.seed(2184)
v5.ns20.ne.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.225), vars=c(2,10), contrast=c(-1,1))
v5.ns50.ne.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.225), vars=c(2,10), contrast=c(-1,1))
v5.ns100.ne.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.225), vars=c(2,10), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(10))
set.seed(2184)
v5.ns20.ne.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.960), vars=c(2,10), contrast=c(-1,1))
v5.ns50.ne.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.960), vars=c(2,10), contrast=c(-1,1))
v5.ns100.ne.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.960), vars=c(2,10), contrast=c(-1,1))


## Save simulations
save(v5.ns20.ne.nullT,v5.ns50.ne.nullT,v5.ns100.ne.nullT,v5.ns20.ne.smalld,v5.ns50.ne.smalld,v5.ns100.ne.smalld,v5.ns20.ne.midd,v5.ns50.ne.midd,v5.ns100.ne.midd, v5.ns20.ne.bigd,v5.ns50.ne.bigd,v5.ns100.ne.bigd, file='/users/joshwondra/R-projects/Welch rule/v5NeSeed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v5NeSeed2184.Rdata')


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


##### Save tables from equal Ns, 5x variances simulations #####
save(obs.coverage.v5.ne, reject.null.v5.ne, file='/users/joshwondra/R-projects/Welch rule/v5NeSeed2184Tables.Rdata')
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


##### Save tables from 1.5x Ns, equal variances simulations #####
save(obs.coverage.ve.1.5n, reject.null.ve.1.5n, file='/users/joshwondra/R-projects/Welch rule/ve1andhalfnSeed2184Tables.Rdata')
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


##### Save tables from double Ns, equal variances simulations #####
save(obs.coverage.ve.2n, reject.null.ve.2n, file='/users/joshwondra/R-projects/Welch rule/ve2nSeed2184Tables.Rdata')
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
v2.ns20.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.36), vars=c(2,4), contrast=c(-1,1))
v2.ns50.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.36), vars=c(2,4), contrast=c(-1,1))
v2.ns100.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.36), vars=c(2,4), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),30,sqrt(4))
cohen.diff(.5,50,sqrt(2),75,sqrt(4))
cohen.diff(.5,100,sqrt(2),150,sqrt(4))
set.seed(2184)
v2.ns20.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.90), vars=c(2,4), contrast=c(-1,1))
v2.ns50.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.89), vars=c(2,4), contrast=c(-1,1))
v2.ns100.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.89), vars=c(2,4), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),30,sqrt(4))
cohen.diff(.8,50,sqrt(2),75,sqrt(4))
cohen.diff(.8,100,sqrt(2),150,sqrt(4))
set.seed(2184)
v2.ns20.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.43), vars=c(2,4), contrast=c(-1,1))
v2.ns50.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.43), vars=c(2,4), contrast=c(-1,1))
v2.ns100.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.43), vars=c(2,4), contrast=c(-1,1))



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


##### Save tables from 1.5x Ns, double variances simulations #####
save(obs.coverage.v2.1.5n.ssv, reject.null.v2.1.5n.ssv, file='/users/joshwondra/R-projects/Welch rule/v21andhalfnSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v21andhalfnSeed2184Tables.Rdata')











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
v5.ns20.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.52), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.52), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.52), vars=c(2,10), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),30,sqrt(10))
cohen.diff(.5,50,sqrt(2),75,sqrt(10))
cohen.diff(.5,100,sqrt(2),150,sqrt(10))
set.seed(2184)
v5.ns20.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.31), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.31), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.30), vars=c(2,10), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),30,sqrt(10))
cohen.diff(.8,50,sqrt(2),75,sqrt(10))
cohen.diff(.8,100,sqrt(2),150,sqrt(10))
set.seed(2184)
v5.ns20.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,8.09), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,8.09), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,8.09), vars=c(2,10), contrast=c(-1,1))



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


##### Save tables from 1.5x Ns, 5x variances simulations #####
save(obs.coverage.v5.1.5n.ssv, reject.null.v5.1.5n.ssv, file='/users/joshwondra/R-projects/Welch rule/v51andhalfnSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v51andhalfnSeed2184Tables.Rdata')





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
v2.ns20.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.37), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.37), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.37), vars=c(2,4), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),40,sqrt(4))
cohen.diff(.5,50,sqrt(2),100,sqrt(4))
cohen.diff(.5,100,sqrt(2),200,sqrt(4))
set.seed(2184)
v2.ns20.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.91), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.91), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.91), vars=c(2,4), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),40,sqrt(4))
cohen.diff(.8,50,sqrt(2),100,sqrt(4))
cohen.diff(.8,100,sqrt(2),200,sqrt(4))
set.seed(2184)
v2.ns20.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.46), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.46), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.46), vars=c(2,4), contrast=c(-1,1))



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


##### Save tables from double Ns, double variances simulations #####
save(obs.coverage.v2.2n.ssv, reject.null.v2.2n.ssv, file='/users/joshwondra/R-projects/Welch rule/v22nSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v22nSeed2184Tables.Rdata')











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
v5.ns20.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.54), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.54), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.54), vars=c(2,10), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),40,sqrt(10))
cohen.diff(.5,50,sqrt(2),100,sqrt(10))
cohen.diff(.5,100,sqrt(2),200,sqrt(10))
set.seed(2184)
v5.ns20.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.36), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.36), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.35), vars=c(2,10), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),40,sqrt(10))
cohen.diff(.8,50,sqrt(2),100,sqrt(10))
cohen.diff(.8,100,sqrt(2),200,sqrt(10))
set.seed(2184)
v5.ns20.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,8.17), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,8.17), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,8.17), vars=c(2,10), contrast=c(-1,1))



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


##### Save tables from double Ns, 5x variances simulations #####
save(obs.coverage.v5.2n.ssv, reject.null.v5.2n.ssv, file='/users/joshwondra/R-projects/Welch rule/v52nSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v52nSeed2184Tables.Rdata')







##### Different Ns and variances, small group small variances #####

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
v2.ns20.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.33), vars=c(4,2), contrast=c(-1,1))
v2.ns50.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.33), vars=c(4,2), contrast=c(-1,1))
v2.ns100.1.5n.bsv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.33), vars=c(4,2), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(4),30,sqrt(2))
cohen.diff(.5,50,sqrt(4),75,sqrt(2))
cohen.diff(.5,100,sqrt(4),150,sqrt(2))
set.seed(2184)
v2.ns20.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.84), vars=c(4,2), contrast=c(-1,1))
v2.ns50.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.84), vars=c(4,2), contrast=c(-1,1))
v2.ns100.1.5n.bsv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.84), vars=c(4,2), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(4),30,sqrt(2))
cohen.diff(.8,50,sqrt(4),75,sqrt(2))
cohen.diff(.8,100,sqrt(4),150,sqrt(2))
set.seed(2184)
v2.ns20.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.34), vars=c(4,2), contrast=c(-1,1))
v2.ns50.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.34), vars=c(4,2), contrast=c(-1,1))
v2.ns100.1.5n.bsv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.34), vars=c(4,2), contrast=c(-1,1))



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


##### Save tables from 1.5x Ns, double variances simulations #####
save(obs.coverage.v2.1.5n.bsv, reject.null.v2.1.5n.bsv, file='/users/joshwondra/R-projects/Welch rule/v21andhalfnBSVSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v21andhalfnBSVSeed2184Tables.Rdata')











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
v5.ns20.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6.52), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6.52), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6.52), vars=c(2,10), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),30,sqrt(10))
cohen.diff(.5,50,sqrt(2),75,sqrt(10))
cohen.diff(.5,100,sqrt(2),150,sqrt(10))
set.seed(2184)
v5.ns20.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,7.31), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,7.31), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,7.30), vars=c(2,10), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),30,sqrt(10))
cohen.diff(.8,50,sqrt(2),75,sqrt(10))
cohen.diff(.8,100,sqrt(2),150,sqrt(10))
set.seed(2184)
v5.ns20.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,8.09), vars=c(2,10), contrast=c(-1,1))
v5.ns50.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,8.09), vars=c(2,10), contrast=c(-1,1))
v5.ns100.1.5n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,8.09), vars=c(2,10), contrast=c(-1,1))



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


##### Save tables from 1.5x Ns, 5x variances simulations #####
save(obs.coverage.v5.1.5n.ssv, reject.null.v5.1.5n.ssv, file='/users/joshwondra/R-projects/Welch rule/v51andhalfnSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v51andhalfnSeed2184Tables.Rdata')





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
v2.ns20.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.37), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.37), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.37), vars=c(2,4), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),40,sqrt(4))
cohen.diff(.5,50,sqrt(2),100,sqrt(4))
cohen.diff(.5,100,sqrt(2),200,sqrt(4))
set.seed(2184)
v2.ns20.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.91), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.91), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.91), vars=c(2,4), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),40,sqrt(4))
cohen.diff(.8,50,sqrt(2),100,sqrt(4))
cohen.diff(.8,100,sqrt(2),200,sqrt(4))
set.seed(2184)
v2.ns20.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.46), vars=c(2,4), contrast=c(-1,1))
v2.ns50.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.46), vars=c(2,4), contrast=c(-1,1))
v2.ns100.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.46), vars=c(2,4), contrast=c(-1,1))



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


##### Save tables from double Ns, double variances simulations #####
save(obs.coverage.v2.2n.ssv, reject.null.v2.2n.ssv, file='/users/joshwondra/R-projects/Welch rule/v22nSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v22nSeed2184Tables.Rdata')











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
v5.ns20.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,6.54), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,6.54), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.smalld <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,6.54), vars=c(2,10), contrast=c(-1,1))

# mid d
cohen.diff(.5,20,sqrt(2),40,sqrt(10))
cohen.diff(.5,50,sqrt(2),100,sqrt(10))
cohen.diff(.5,100,sqrt(2),200,sqrt(10))
set.seed(2184)
v5.ns20.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,7.36), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,7.36), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.midd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,7.35), vars=c(2,10), contrast=c(-1,1))

# big d
cohen.diff(.8,20,sqrt(2),40,sqrt(10))
cohen.diff(.8,50,sqrt(2),100,sqrt(10))
cohen.diff(.8,100,sqrt(2),200,sqrt(10))
set.seed(2184)
v5.ns20.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(20,40), means=c(6,8.17), vars=c(2,10), contrast=c(-1,1))
v5.ns50.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(50,100), means=c(6,8.17), vars=c(2,10), contrast=c(-1,1))
v5.ns100.2n.ssv.bigd <- t.compare(nsims=10000, Ns=c(100,200), means=c(6,8.17), vars=c(2,10), contrast=c(-1,1))



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


##### Save tables from double Ns, 5x variances simulations #####
save(obs.coverage.v5.2n.ssv, reject.null.v5.2n.ssv, file='/users/joshwondra/R-projects/Welch rule/v52nSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v52nSeed2184Tables.Rdata')







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
exp.rejects.2n[,2,] <- c(.111,.209,.339)
exp.rejects.2n[,3,] <- c(.435, .818, .971)
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
