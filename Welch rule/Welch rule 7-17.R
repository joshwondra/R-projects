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
t.compare <- function(nsims, Ns, means, vars, contrast) 
{    
    sims <- vector('list', nsims)
    group <- rep(1:length(Ns), Ns)   #vector of length N with group codes
    
    ### Create simulated data
    for(i in 1:nsims) 
    {   
        dv <- NULL
        for(j in 1:length(Ns))
        {
            dv <- c(dv, rnorm(n=Ns[j], mean=means[j], sd=sqrt(vars[j])))
        }
        
        fit <- t.contrast(dv, group, contrast)
        
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
        
        current.sim <- list(data.frame(dv, group), contrast, ihat, se.classic, t.classic, df.classic, p.classic, se.welch, t.welch, df.welch, p.welch)
        names(current.sim) <- c('data', 'contrast', 'ihat', 'se.classic', 't.classic', 'df.classic', 'p.classic', 'se.welch', 't.welch', 'df.welch', 'p.welch')
        
        sims[[i]] <- current.sim
    }
    
    names(sims)[1:length(sims)] <- paste('sim',1:length(sims),sep='')
    
    # store proportion of rejected null hypotheses
    classic.reject <- prop.table(table(lapply(sims, '[[', 7)<=.05))[[2]]
    welch.reject <- prop.table(table(lapply(sims, '[[', 11)<=.05))[[2]]
    
    # compute coverage rate
    true.ihat <- contrast %*% means
          
    obs.ihat <- data.frame(lapply(sims, '[[', 3))
    classic.df <- data.frame(lapply(sims, '[[', 6))
    welch.df <- data.frame(lapply(sims, '[[', 10))
    classic.ses <- data.frame(lapply(sims, '[[', 4))
    welch.ses <- data.frame(lapply(sims, '[[', 8))
        
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
    return(list=c(classic.reject=classic.reject, welch.reject=welch.reject, classic.coverage=classic.coverage, welch.coverage=welch.coverage, sims))
}



##### SIMULATIONS #####

#Compute mean differences to obtain specific effect sizes
cohen.diff <- function(d, N1, SD1, N2, SD2) #insert the mean, sample size, and standard deviation for each group
{
    poolSD <- sqrt(((N1-1)*SD1^2+(N2-1)*SD2^2)/(N1+N2-2)) #computes the pooled standard deviation
    diff <- d*poolSD #computes difference needed
    print(diff) #displays difference needed
}


##### Equal variances simulations #####

set.seed(2184)

# true null
ve.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ve.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,2), contrast=c(-1,1))

# small effect
cohen.diff(.2, 20, sqrt(2), 20, sqrt(2))
ve.ns20.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ve.ns50.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ve.ns100.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ve.ns200.smalld <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(2))
ve.ns20.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ve.ns50.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ve.ns100.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ve.ns200.midd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(2))
ve.ns20.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ve.ns50.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ve.ns100.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ve.ns200.bigd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))


## Save simulations
save(ve.ns20.nullT,ve.ns50.nullT,ve.ns100.nullT,ve.ns200.nullT,ve.ns20.smalld,ve.ns50.smalld,ve.ns100.smalld,ve.ns200.smalld,ve.ns20.midd,ve.ns50.midd,ve.ns100.midd,ve.ns200.midd, ve.ns20.bigd,ve.ns50.bigd,ve.ns100.bigd,ve.ns200.bigd, file='/users/joshwondra/R-projects/Welch rule/veSeed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/veSeed2184.Rdata')


##### Equal variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.ve <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.ve[,1,1] <- c(ve.ns20.nullT$classic.reject, ve.ns50.nullT$classic.reject, ve.ns100.nullT$classic.reject, ve.ns200.nullT$classic.reject)
reject.null.ve[,1,2] <- c(ve.ns20.nullT$welch.reject, ve.ns50.nullT$welch.reject, ve.ns100.nullT$welch.reject, ve.ns200.nullT$welch.reject)

# small effect
reject.null.ve[,2,1] <- c(ve.ns20.smalld$classic.reject, ve.ns50.smalld$classic.reject, ve.ns100.smalld$classic.reject, ve.ns200.smalld$classic.reject)
reject.null.ve[,2,2] <- c(ve.ns20.smalld$welch.reject, ve.ns50.smalld$welch.reject, ve.ns100.smalld$welch.reject, ve.ns200.smalld$welch.reject)

# medium effect
reject.null.ve[,3,1] <- c(ve.ns20.midd$classic.reject, ve.ns50.midd$classic.reject, ve.ns100.midd$classic.reject, ve.ns200.midd$classic.reject)
reject.null.ve[,3,2] <- c(ve.ns20.midd$welch.reject, ve.ns50.midd$welch.reject, ve.ns100.midd$welch.reject, ve.ns200.midd$welch.reject)

# large effect
reject.null.ve[,4,1] <- c(ve.ns20.bigd$classic.reject, ve.ns50.bigd$classic.reject, ve.ns100.bigd$classic.reject, ve.ns200.bigd$classic.reject)
reject.null.ve[,4,2] <- c(ve.ns20.bigd$welch.reject, ve.ns50.bigd$welch.reject, ve.ns100.bigd$welch.reject, ve.ns200.bigd$welch.reject)


##### Equal variances coverage rate #####

## store observed coverage rates
obs.coverage.ve <- array(dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.ve[,1,1] <- c(ve.ns20.nullT$classic.coverage, ve.ns50.nullT$classic.coverage, ve.ns100.nullT$classic.coverage, ve.ns200.nullT$classic.coverage)
obs.coverage.ve[,1,2] <- c(ve.ns20.nullT$welch.coverage, ve.ns50.nullT$welch.coverage, ve.ns100.nullT$welch.coverage, ve.ns200.nullT$welch.coverage)

# small effect
obs.coverage.ve[,2,1] <- c(ve.ns20.smalld$classic.coverage, ve.ns50.smalld$classic.coverage, ve.ns100.smalld$classic.coverage, ve.ns200.smalld$classic.coverage)
obs.coverage.ve[,2,2] <- c(ve.ns20.smalld$welch.coverage, ve.ns50.smalld$welch.coverage, ve.ns100.smalld$welch.coverage, ve.ns200.smalld$welch.coverage)

# medium effect
obs.coverage.ve[,3,1] <- c(ve.ns20.midd$classic.coverage, ve.ns50.midd$classic.coverage, ve.ns100.midd$classic.coverage, ve.ns200.midd$classic.coverage)
obs.coverage.ve[,3,2] <- c(ve.ns20.midd$welch.coverage, ve.ns50.midd$welch.coverage, ve.ns100.midd$welch.coverage, ve.ns200.midd$welch.coverage)

# big effect
obs.coverage.ve[,4,1] <- c(ve.ns20.bigd$classic.coverage, ve.ns50.bigd$classic.coverage, ve.ns100.bigd$classic.coverage, ve.ns200.bigd$classic.coverage)
obs.coverage.ve[,4,2] <- c(ve.ns20.bigd$welch.coverage, ve.ns50.bigd$welch.coverage, ve.ns100.bigd$welch.coverage, ve.ns200.bigd$welch.coverage)


##### Save tables from equal variances simulations #####
save(obs.coverage.ve, reject.null.ve, file='/users/joshwondra/R-projects/Welch rule/veSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/veSeed2184Tables.Rdata')



##### 1.5 times different variances simulations #####

set.seed(2184)

# true null
v1.5.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,3), contrast=c(-1,1))

# small effect
cohen.diff(.2, 20, sqrt(2), 20, sqrt(3))
v1.5.ns20.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.316), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.316), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.316), vars=c(2,3), contrast=c(-1,1))
v1.5.ns200.smalld <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.316), vars=c(2,3), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(3))
v1.5.ns20.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.791), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.791), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.791), vars=c(2,3), contrast=c(-1,1))
v1.5.ns200.midd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.791), vars=c(2,3), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(3))
v1.5.ns20.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.265), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.265), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.265), vars=c(2,3), contrast=c(-1,1))
v1.5.ns200.bigd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,7.265), vars=c(2,3), contrast=c(-1,1))


## Save simulations
save(v1.5.ns20.nullT,v1.5.ns50.nullT,v1.5.ns100.nullT,v1.5.ns200.nullT,v1.5.ns20.smalld,v1.5.ns50.smalld,v1.5.ns100.smalld,v1.5.ns200.smalld,v1.5.ns20.midd,v1.5.ns50.midd,v1.5.ns100.midd,v1.5.ns200.midd, v1.5.ns20.bigd,v1.5.ns50.bigd,v1.5.ns100.bigd,v1.5.ns200.bigd, file='/users/joshwondra/R-projects/Welch rule/v1andhalfSeed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v1andhalfSeed2184.Rdata')



##### 1.5 variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v1.5 <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v1.5[,1,1] <- c(v1.5.ns20.nullT$classic.reject, v1.5.ns50.nullT$classic.reject, v1.5.ns100.nullT$classic.reject, v1.5.ns200.nullT$classic.reject)
reject.null.v1.5[,1,2] <- c(v1.5.ns20.nullT$welch.reject, v1.5.ns50.nullT$welch.reject, v1.5.ns100.nullT$welch.reject, v1.5.ns200.nullT$welch.reject)

# small effect
reject.null.v1.5[,2,1] <- c(v1.5.ns20.smalld$classic.reject, v1.5.ns50.smalld$classic.reject, v1.5.ns100.smalld$classic.reject, v1.5.ns200.smalld$classic.reject)
reject.null.v1.5[,2,2] <- c(v1.5.ns20.smalld$welch.reject, v1.5.ns50.smalld$welch.reject, v1.5.ns100.smalld$welch.reject, v1.5.ns200.smalld$welch.reject)

# medium effect
reject.null.v1.5[,3,1] <- c(v1.5.ns20.midd$classic.reject, v1.5.ns50.midd$classic.reject, v1.5.ns100.midd$classic.reject, v1.5.ns200.midd$classic.reject)
reject.null.v1.5[,3,2] <- c(v1.5.ns20.midd$welch.reject, v1.5.ns50.midd$welch.reject, v1.5.ns100.midd$welch.reject, v1.5.ns200.midd$welch.reject)

# large effect
reject.null.v1.5[,4,1] <- c(v1.5.ns20.bigd$classic.reject, v1.5.ns50.bigd$classic.reject, v1.5.ns100.bigd$classic.reject, v1.5.ns200.bigd$classic.reject)
reject.null.v1.5[,4,2] <- c(v1.5.ns20.bigd$welch.reject, v1.5.ns50.bigd$welch.reject, v1.5.ns100.bigd$welch.reject, v1.5.ns200.bigd$welch.reject)


##### 1.5 variances coverage rate #####

## store observed coverage rates
obs.coverage.v1.5 <- array(dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v1.5[,1,1] <- c(v1.5.ns20.nullT$classic.coverage, v1.5.ns50.nullT$classic.coverage, v1.5.ns100.nullT$classic.coverage, v1.5.ns200.nullT$classic.coverage)
obs.coverage.v1.5[,1,2] <- c(v1.5.ns20.nullT$welch.coverage, v1.5.ns50.nullT$welch.coverage, v1.5.ns100.nullT$welch.coverage, v1.5.ns200.nullT$welch.coverage)

# small effect
obs.coverage.v1.5[,2,1] <- c(v1.5.ns20.smalld$classic.coverage, v1.5.ns50.smalld$classic.coverage, v1.5.ns100.smalld$classic.coverage, v1.5.ns200.smalld$classic.coverage)
obs.coverage.v1.5[,2,2] <- c(v1.5.ns20.smalld$welch.coverage, v1.5.ns50.smalld$welch.coverage, v1.5.ns100.smalld$welch.coverage, v1.5.ns200.smalld$welch.coverage)

# medium effect
obs.coverage.v1.5[,3,1] <- c(v1.5.ns20.midd$classic.coverage, v1.5.ns50.midd$classic.coverage, v1.5.ns100.midd$classic.coverage, v1.5.ns200.midd$classic.coverage)
obs.coverage.v1.5[,3,2] <- c(v1.5.ns20.midd$welch.coverage, v1.5.ns50.midd$welch.coverage, v1.5.ns100.midd$welch.coverage, v1.5.ns200.midd$welch.coverage)

# big effect
obs.coverage.v1.5[,4,1] <- c(v1.5.ns20.bigd$classic.coverage, v1.5.ns50.bigd$classic.coverage, v1.5.ns100.bigd$classic.coverage, v1.5.ns200.bigd$classic.coverage)
obs.coverage.v1.5[,4,2] <- c(v1.5.ns20.bigd$welch.coverage, v1.5.ns50.bigd$welch.coverage, v1.5.ns100.bigd$welch.coverage, v1.5.ns200.bigd$welch.coverage)


##### Save tables from 1.5 variances simulations #####
save(obs.coverage.v1.5, reject.null.v1.5, file='/users/joshwondra/R-projects/Welch rule/v1halfSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v1halfSeed2184Tables.Rdata')




##### Double variance simulations #####

set.seed(2184)

# true null
v2.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,4), contrast=c(-1,1))

# small effect
cohen.diff(.2, 20, sqrt(2), 20, sqrt(4))
v2.ns20.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.346), vars=c(2,4), contrast=c(-1,1))
v2.ns50.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.346), vars=c(2,4), contrast=c(-1,1))
v2.ns100.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.346), vars=c(2,4), contrast=c(-1,1))
v2.ns200.smalld <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.346), vars=c(2,4), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(4))
v2.ns20.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.866), vars=c(2,4), contrast=c(-1,1))
v2.ns50.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.866), vars=c(2,4), contrast=c(-1,1))
v2.ns100.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.866), vars=c(2,4), contrast=c(-1,1))
v2.ns200.midd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.866), vars=c(2,4), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(4))
v2.ns20.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.386), vars=c(2,4), contrast=c(-1,1))
v2.ns50.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.386), vars=c(2,4), contrast=c(-1,1))
v2.ns100.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.386), vars=c(2,4), contrast=c(-1,1))
v2.ns200.bigd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,7.386), vars=c(2,4), contrast=c(-1,1))


## Save simulations
save(v2.ns20.nullT,v2.ns50.nullT,v2.ns100.nullT,v2.ns200.nullT,v2.ns20.smalld,v2.ns50.smalld,v2.ns100.smalld,v2.ns200.smalld,v2.ns20.midd,v2.ns50.midd,v2.ns100.midd,v2.ns200.midd, v2.ns20.bigd,v2.ns50.bigd,v2.ns100.bigd,v2.ns200.bigd, file='/users/joshwondra/R-projects/Welch rule/v2Seed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v2Seed2184.Rdata')


##### Double variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v2 <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v2[,1,1] <- c(v2.ns20.nullT$classic.reject, v2.ns50.nullT$classic.reject, v2.ns100.nullT$classic.reject, v2.ns200.nullT$classic.reject)
reject.null.v2[,1,2] <- c(v2.ns20.nullT$welch.reject, v2.ns50.nullT$welch.reject, v2.ns100.nullT$welch.reject, v2.ns200.nullT$welch.reject)

# small effect
reject.null.v2[,2,1] <- c(v2.ns20.smalld$classic.reject, v2.ns50.smalld$classic.reject, v2.ns100.smalld$classic.reject, v2.ns200.smalld$classic.reject)
reject.null.v2[,2,2] <- c(v2.ns20.smalld$welch.reject, v2.ns50.smalld$welch.reject, v2.ns100.smalld$welch.reject, v2.ns200.smalld$welch.reject)

# medium effect
reject.null.v2[,3,1] <- c(v2.ns20.midd$classic.reject, v2.ns50.midd$classic.reject, v2.ns100.midd$classic.reject, v2.ns200.midd$classic.reject)
reject.null.v2[,3,2] <- c(v2.ns20.midd$welch.reject, v2.ns50.midd$welch.reject, v2.ns100.midd$welch.reject, v2.ns200.midd$welch.reject)

# large effect
reject.null.v2[,4,1] <- c(v2.ns20.bigd$classic.reject, v2.ns50.bigd$classic.reject, v2.ns100.bigd$classic.reject, v2.ns200.bigd$classic.reject)
reject.null.v2[,4,2] <- c(v2.ns20.bigd$welch.reject, v2.ns50.bigd$welch.reject, v2.ns100.bigd$welch.reject, v2.ns200.bigd$welch.reject)


##### Double variances coverage rate #####

## store observed coverage rates
obs.coverage.v2 <- array(dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v2[,1,1] <- c(v2.ns20.nullT$classic.coverage, v2.ns50.nullT$classic.coverage, v2.ns100.nullT$classic.coverage, v2.ns200.nullT$classic.coverage)
obs.coverage.v2[,1,2] <- c(v2.ns20.nullT$welch.coverage, v2.ns50.nullT$welch.coverage, v2.ns100.nullT$welch.coverage, v2.ns200.nullT$welch.coverage)

# small effect
obs.coverage.v2[,2,1] <- c(v2.ns20.smalld$classic.coverage, v2.ns50.smalld$classic.coverage, v2.ns100.smalld$classic.coverage, v2.ns200.smalld$classic.coverage)
obs.coverage.v2[,2,2] <- c(v2.ns20.smalld$welch.coverage, v2.ns50.smalld$welch.coverage, v2.ns100.smalld$welch.coverage, v2.ns200.smalld$welch.coverage)

# medium effect
obs.coverage.v2[,3,1] <- c(v2.ns20.midd$classic.coverage, v2.ns50.midd$classic.coverage, v2.ns100.midd$classic.coverage, v2.ns200.midd$classic.coverage)
obs.coverage.v2[,3,2] <- c(v2.ns20.midd$welch.coverage, v2.ns50.midd$welch.coverage, v2.ns100.midd$welch.coverage, v2.ns200.midd$welch.coverage)

# big effect
obs.coverage.v2[,4,1] <- c(v2.ns20.bigd$classic.coverage, v2.ns50.bigd$classic.coverage, v2.ns100.bigd$classic.coverage, v2.ns200.bigd$classic.coverage)
obs.coverage.v2[,4,2] <- c(v2.ns20.bigd$welch.coverage, v2.ns50.bigd$welch.coverage, v2.ns100.bigd$welch.coverage, v2.ns200.bigd$welch.coverage)


##### Save tables from double variances simulations #####
save(obs.coverage.v2, reject.null.v2, file='/users/joshwondra/R-projects/Welch rule/v2Seed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v2Seed2184Tables.Rdata')



##### Five times different variances simulations #####

set.seed(2184)

# true null
v5.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,10), contrast=c(-1,1))

# small effect
cohen.diff(.2, 20, sqrt(2), 20, sqrt(10))
v5.ns20.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.490), vars=c(2,10), contrast=c(-1,1))
v5.ns50.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.490), vars=c(2,10), contrast=c(-1,1))
v5.ns100.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.490), vars=c(2,10), contrast=c(-1,1))
v5.ns200.smalld <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.490), vars=c(2,10), contrast=c(-1,1))

# medium effect
cohen.diff(.5, 20, sqrt(2), 20, sqrt(10))
v5.ns20.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.225), vars=c(2,10), contrast=c(-1,1))
v5.ns50.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.225), vars=c(2,10), contrast=c(-1,1))
v5.ns100.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.225), vars=c(2,10), contrast=c(-1,1))
v5.ns200.midd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,7.225), vars=c(2,10), contrast=c(-1,1))

# large effect
cohen.diff(.8, 20, sqrt(2), 20, sqrt(10))
v5.ns20.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.960), vars=c(2,10), contrast=c(-1,1))
v5.ns50.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.960), vars=c(2,10), contrast=c(-1,1))
v5.ns100.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.960), vars=c(2,10), contrast=c(-1,1))
v5.ns200.bigd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,7.960), vars=c(2,10), contrast=c(-1,1))


## Save simulations
save(v5.ns20.nullT,v5.ns50.nullT,v5.ns100.nullT,v5.ns200.nullT,v5.ns20.smalld,v5.ns50.smalld,v5.ns100.smalld,v5.ns200.smalld,v5.ns20.midd,v5.ns50.midd,v5.ns100.midd,v5.ns200.midd, v5.ns20.bigd,v5.ns50.bigd,v5.ns100.bigd,v5.ns200.bigd, file='/users/joshwondra/R-projects/Welch rule/v5Seed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v5Seed2184.Rdata')


##### Five times variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v5 <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
reject.null.v5[,1,1] <- c(v5.ns20.nullT$classic.reject, v5.ns50.nullT$classic.reject, v5.ns100.nullT$classic.reject, v5.ns200.nullT$classic.reject)
reject.null.v5[,1,2] <- c(v5.ns20.nullT$welch.reject, v5.ns50.nullT$welch.reject, v5.ns100.nullT$welch.reject, v5.ns200.nullT$welch.reject)

# small effect
reject.null.v5[,2,1] <- c(v5.ns20.smalld$classic.reject, v5.ns50.smalld$classic.reject, v5.ns100.smalld$classic.reject, v5.ns200.smalld$classic.reject)
reject.null.v5[,2,2] <- c(v5.ns20.smalld$welch.reject, v5.ns50.smalld$welch.reject, v5.ns100.smalld$welch.reject, v5.ns200.smalld$welch.reject)

# medium effect
reject.null.v5[,3,1] <- c(v5.ns20.midd$classic.reject, v5.ns50.midd$classic.reject, v5.ns100.midd$classic.reject, v5.ns200.midd$classic.reject)
reject.null.v5[,3,2] <- c(v5.ns20.midd$welch.reject, v5.ns50.midd$welch.reject, v5.ns100.midd$welch.reject, v5.ns200.midd$welch.reject)

# large effect
reject.null.v5[,4,1] <- c(v5.ns20.bigd$classic.reject, v5.ns50.bigd$classic.reject, v5.ns100.bigd$classic.reject, v5.ns200.bigd$classic.reject)
reject.null.v5[,4,2] <- c(v5.ns20.bigd$welch.reject, v5.ns50.bigd$welch.reject, v5.ns100.bigd$welch.reject, v5.ns200.bigd$welch.reject)


##### Five times variances coverage rate #####

## store observed coverage rates
obs.coverage.v5 <- array(dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
obs.coverage.v5[,1,1] <- c(v5.ns20.nullT$classic.coverage, v5.ns50.nullT$classic.coverage, v5.ns100.nullT$classic.coverage, v5.ns200.nullT$classic.coverage)
obs.coverage.v5[,1,2] <- c(v5.ns20.nullT$welch.coverage, v5.ns50.nullT$welch.coverage, v5.ns100.nullT$welch.coverage, v5.ns200.nullT$welch.coverage)

# small effect
obs.coverage.v5[,2,1] <- c(v5.ns20.smalld$classic.coverage, v5.ns50.smalld$classic.coverage, v5.ns100.smalld$classic.coverage, v5.ns200.smalld$classic.coverage)
obs.coverage.v5[,2,2] <- c(v5.ns20.smalld$welch.coverage, v5.ns50.smalld$welch.coverage, v5.ns100.smalld$welch.coverage, v5.ns200.smalld$welch.coverage)

# medium effect
obs.coverage.v5[,3,1] <- c(v5.ns20.midd$classic.coverage, v5.ns50.midd$classic.coverage, v5.ns100.midd$classic.coverage, v5.ns200.midd$classic.coverage)
obs.coverage.v5[,3,2] <- c(v5.ns20.midd$welch.coverage, v5.ns50.midd$welch.coverage, v5.ns100.midd$welch.coverage, v5.ns200.midd$welch.coverage)

# big effect
obs.coverage.v5[,4,1] <- c(v5.ns20.bigd$classic.coverage, v5.ns50.bigd$classic.coverage, v5.ns100.bigd$classic.coverage, v5.ns200.bigd$classic.coverage)
obs.coverage.v5[,4,2] <- c(v5.ns20.bigd$welch.coverage, v5.ns50.bigd$welch.coverage, v5.ns100.bigd$welch.coverage, v5.ns200.bigd$welch.coverage)


##### Save tables from five times variances simulations #####
save(obs.coverage.v5, reject.null.v5, file='/users/joshwondra/R-projects/Welch rule/v1halfSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v1halfSeed2184Tables.Rdata')





## expected false positive rate and power

exp.rejects <- array(rep(NA, 32), dim=c(4,4,1), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8")))
exp.rejects[,1,] <- rep(.05,4)
exp.rejects[,2,] <- c(.095, .168, .291, .514)
exp.rejects[,3,] <- c(.338, .697, .940, 1)
exp.rejects[,4,] <- c(.693, .977, 1, 1)

## expected coverage rates

exp.coverage <- array(rep(.95, 64), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c('classic', 'welch')))



















##### True nulls, different Ns #####

# 1.5 variances, 5% different Ns, small group small variance
cohen.diff(.2,20,sqrt(2),21,sqrt(3))
cohen.diff(.5,20,sqrt(2),21,sqrt(3))
cohen.diff(.8,20,sqrt(2),21,sqrt(3))
set.seed(2184)
v1.5.ns20.5psmv.nullT <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns20.5psmv.smalld <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,6.32), vars=c(2,3), contrast=c(-1,1))
v1.5.ns20.5psmv.midd <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,6.79), vars=c(2,3), contrast=c(-1,1))
v1.5.ns20.5psmv.bigd <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,7.27), vars=c(2,3), contrast=c(-1,1))

cohen.diff(.2,50,sqrt(2),53,sqrt(3))
cohen.diff(.5,50,sqrt(2),53,sqrt(3))
cohen.diff(.8,50,sqrt(2),53,sqrt(3))
set.seed(2184)
v1.5.ns50.5psmv.nullT <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.5psmv.smalld <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,6.32), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.5psmv.midd <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,6.79), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.5psmv.bigd <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,7.27), vars=c(2,3), contrast=c(-1,1))

cohen.diff(.2,100,sqrt(2),105,sqrt(3))
cohen.diff(.5,100,sqrt(2),105,sqrt(3))
cohen.diff(.8,100,sqrt(2),105,sqrt(3))
set.seed(2184)
v1.5.ns100.5psmv.nullT <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.5psmv.smalld <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,6.32), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.5psmv.midd <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,6.79), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.5psmv.bigd <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,7.27), vars=c(2,3), contrast=c(-1,1))

## array to store proportion of rejected null hypotheses
v1.5.5psmv.rejects <- array(rep(NA, 32), dim=c(3,4,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

# true null
v1.5.5psmv.rejects[,1,1] <- c(v1.5.ns20.5psmv.nullT$classic.reject, v1.5.ns50.5psmv.nullT$classic.reject, v1.5.ns100.5psmv.nullT$classic.reject)
v1.5.5psmv.rejects[,1,2] <- c(v1.5.ns20.5psmv.nullT$welch.reject, v1.5.ns50.5psmv.nullT$welch.reject, v1.5.ns100.5psmv.nullT$welch.reject)

v1.5.5psmv.rejects[,2,1] <- c(v1.5.ns20.5psmv.smalld$classic.reject, v1.5.ns50.5psmv.smalld$classic.reject, v1.5.ns100.5psmv.smalld$classic.reject)
v1.5.5psmv.rejects[,2,2] <- c(v1.5.ns20.5psmv.smalld$welch.reject, v1.5.ns50.5psmv.smalld$welch.reject, v1.5.ns100.5psmv.smalld$welch.reject)

v1.5.5psmv.rejects[,1,1] <- c(v1.5.ns20.5psmv.midd$classic.reject, v1.5.ns50.5psmv.midd$classic.reject, v1.5.ns100.5psmv.midd$classic.reject)
v1.5.5psmv.rejects[,1,2] <- c(v1.5.ns20.5psmv.midd$welch.reject, v1.5.ns50.5psmv.midd$welch.reject, v1.5.ns100.5psmv.midd$welch.reject)

v1.5.5psmv.rejects[,1,1] <- c(v1.5.ns20.5psmv.bigd$classic.reject, v1.5.ns50.5psmv.bigd$classic.reject, v1.5.ns100.5psmv.bigd$classic.reject)
v1.5.5psmv.rejects[,1,2] <- c(v1.5.ns20.5psmv.bigd$welch.reject, v1.5.ns50.5psmv.bigd$welch.reject, v1.5.ns100.5psmv.bigd$welch.reject)

v1.5.5psmv.rejects


# 1.5 variances, 5% different Ns, big group small variance
cohen.diff(.2,20,sqrt(2),21,sqrt(3))
cohen.diff(.5,20,sqrt(2),21,sqrt(3))
cohen.diff(.8,20,sqrt(2),21,sqrt(3))
set.seed(2184)
v1.5.ns20.5pbv.nullT <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns20.5pbv.smalld <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,6.32), vars=c(2,3), contrast=c(-1,1))
v1.5.ns20.5pbv.midd <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,6.79), vars=c(2,3), contrast=c(-1,1))
v1.5.ns20.5pbv.bigd <- t.compare(nsims=10000, Ns=c(20,21), means=c(6,7.27), vars=c(2,3), contrast=c(-1,1))

cohen.diff(.2,50,sqrt(2),53,sqrt(3))
cohen.diff(.5,50,sqrt(2),53,sqrt(3))
cohen.diff(.8,50,sqrt(2),53,sqrt(3))
set.seed(2184)
v1.5.ns50.5pbv.nullT <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.5pbv.smalld <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,6.32), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.5pbv.midd <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,6.79), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50.5pbv.bigd <- t.compare(nsims=10000, Ns=c(50,53), means=c(6,7.27), vars=c(2,3), contrast=c(-1,1))

cohen.diff(.2,100,sqrt(2),105,sqrt(3))
cohen.diff(.5,100,sqrt(2),105,sqrt(3))
cohen.diff(.8,100,sqrt(2),105,sqrt(3))
set.seed(2184)
v1.5.ns100.5pbv.nullT <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.5pbv.smalld <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,6.32), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.5pbv.midd <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,6.79), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100.5pbv.bigd <- t.compare(nsims=10000, Ns=c(100,105), means=c(6,7.27), vars=c(2,3), contrast=c(-1,1))

# 1.5 variances, 10% different Ns
set.seed(2184)
v1.5.ns20n10p.nullT <- t.compare(nsims=10000, Ns=c(20,22), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50n10p.nullT <- t.compare(nsims=10000, Ns=c(50,55), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100n10p.nullT <- t.compare(nsims=10000, Ns=c(100,110), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns200n10p.nullT <- t.compare(nsims=10000, Ns=c(200,220), means=c(6,6), vars=c(2,3), contrast=c(-1,1))

# 1.5 variances, 20% different Ns
set.seed(2184)
v1.5.ns20n20p.nullT <- t.compare(nsims=10000, Ns=c(20,24), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50n20p.nullT <- t.compare(nsims=10000, Ns=c(50,60), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100n20p.nullT <- t.compare(nsims=10000, Ns=c(100,120), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns200n20p.nullT <- t.compare(nsims=10000, Ns=c(200,240), means=c(6,6), vars=c(2,3), contrast=c(-1,1))

# 1.5 variances, 50% different Ns
set.seed(2184)
v1.5.ns20n50p.nullT <- t.compare(nsims=10000, Ns=c(20,30), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns50n50p.nullT <- t.compare(nsims=10000, Ns=c(50,75), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns100n50p.nullT <- t.compare(nsims=10000, Ns=c(100,150), means=c(6,6), vars=c(2,3), contrast=c(-1,1))
v1.5.ns200n50p.nullT <- t.compare(nsims=10000, Ns=c(200,300), means=c(6,6), vars=c(2,3), contrast=c(-1,1))

##### Tables from 1.5 vars, true null, different Ns #####
reject.null.v1.5.diffNs.trueNull <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("5%", "10%", "20%", "50%"), c("classic", "welch")))

#classic
reject.null.v1.5.diffNs.trueNull[1,1,1] <- prop.table(table(lapply(v1.5.ns20n5p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[1,2,1] <- prop.table(table(lapply(v1.5.ns20n10p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[1,3,1] <- prop.table(table(lapply(v1.5.ns20n20p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[1,4,1] <- prop.table(table(lapply(v1.5.ns20n50p.nullT, '[[', 7)<=.05))[[2]]

reject.null.v1.5.diffNs.trueNull[2,1,1] <- prop.table(table(lapply(v1.5.ns50n5p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[2,2,1] <- prop.table(table(lapply(v1.5.ns50n10p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[2,3,1] <- prop.table(table(lapply(v1.5.ns50n20p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[2,4,1] <- prop.table(table(lapply(v1.5.ns50n50p.nullT, '[[', 7)<=.05))[[2]]

reject.null.v1.5.diffNs.trueNull[3,1,1] <- prop.table(table(lapply(v1.5.ns100n5p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[3,2,1] <- prop.table(table(lapply(v1.5.ns100n10p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[3,3,1] <- prop.table(table(lapply(v1.5.ns100n20p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[3,4,1] <- prop.table(table(lapply(v1.5.ns100n50p.nullT, '[[', 7)<=.05))[[2]]

reject.null.v1.5.diffNs.trueNull[4,1,1] <- prop.table(table(lapply(v1.5.ns200n5p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[4,2,1] <- prop.table(table(lapply(v1.5.ns200n10p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[4,3,1] <- prop.table(table(lapply(v1.5.ns200n20p.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[4,4,1] <- prop.table(table(lapply(v1.5.ns200n50p.nullT, '[[', 7)<=.05))[[2]]

#welch
reject.null.v1.5.diffNs.trueNull[1,1,2] <- prop.table(table(lapply(v1.5.ns20n5p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[1,2,2] <- prop.table(table(lapply(v1.5.ns20n10p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[1,3,2] <- prop.table(table(lapply(v1.5.ns20n20p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[1,4,2] <- prop.table(table(lapply(v1.5.ns20n50p.nullT, '[[', 11)<=.05))[[2]]

reject.null.v1.5.diffNs.trueNull[2,1,2] <- prop.table(table(lapply(v1.5.ns50n5p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[2,2,2] <- prop.table(table(lapply(v1.5.ns50n10p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[2,3,2] <- prop.table(table(lapply(v1.5.ns50n20p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[2,4,2] <- prop.table(table(lapply(v1.5.ns50n50p.nullT, '[[', 11)<=.05))[[2]]

reject.null.v1.5.diffNs.trueNull[3,1,2] <- prop.table(table(lapply(v1.5.ns100n5p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[3,2,2] <- prop.table(table(lapply(v1.5.ns100n10p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[3,3,2] <- prop.table(table(lapply(v1.5.ns100n20p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[3,4,2] <- prop.table(table(lapply(v1.5.ns100n50p.nullT, '[[', 11)<=.05))[[2]]

reject.null.v1.5.diffNs.trueNull[4,1,2] <- prop.table(table(lapply(v1.5.ns200n5p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[4,2,2] <- prop.table(table(lapply(v1.5.ns200n10p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[4,3,2] <- prop.table(table(lapply(v1.5.ns200n20p.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5.diffNs.trueNull[4,4,2] <- prop.table(table(lapply(v1.5.ns200n50p.nullT, '[[', 11)<=.05))[[2]]



# 2x variances
set.seed(2184)
v2.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,4), contrast=c(-1,1))
v2.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,4), contrast=c(-1,1))

# 5x variances
set.seed(2184)
v5.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,10), contrast=c(-1,1))
v5.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,10), contrast=c(-1,1))








## Save simulations
save(v1.5.ns20.nullT,v1.5.ns50.nullT,v1.5.ns100.nullT,v1.5.ns200.nullT,v1.5.ns20.smalld,v1.5.ns50.smalld,v1.5.ns100.smalld,v1.5.ns200.smalld,v1.5.ns20.midd,v1.5.ns50.midd,v1.5.ns100.midd,v1.5.ns200.midd, v1.5.ns20.bigd,v1.5.ns50.bigd,v1.5.ns100.bigd,v1.5.ns200.bigd, file='/users/joshwondra/R-projects/Welch rule/v1andhalfSeed2184.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v1andhalfSeed2184.Rdata')


##### 1.5 times different variances false positives and power #####

## array to store proportion of rejected null hypotheses
reject.null.v1.5 <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))
rownames(reject.null.v1.5) <- c("N=20", "N=50", "N=100", "N=200")
colnames(reject.null.v1.5) <- c("d=0", "d=.2", "d=.5", "d=.8")

# true null
reject.null.v1.5[1,1,1] <- prop.table(table(lapply(v1.5.ns20.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5[2,1,1] <- prop.table(table(lapply(v1.5.ns50.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5[3,1,1] <- prop.table(table(lapply(v1.5.ns100.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5[4,1,1] <- prop.table(table(lapply(v1.5.ns200.nullT, '[[', 7)<=.05))[[2]]
reject.null.v1.5[1,1,2] <- prop.table(table(lapply(v1.5.ns20.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5[2,1,2] <- prop.table(table(lapply(v1.5.ns50.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5[3,1,2] <- prop.table(table(lapply(v1.5.ns100.nullT, '[[', 11)<=.05))[[2]]
reject.null.v1.5[4,1,2] <- prop.table(table(lapply(v1.5.ns200.nullT, '[[', 11)<=.05))[[2]]

# small effect
reject.null.v1.5[1,2,1] <- prop.table(table(lapply(v1.5.ns20.smalld, '[[', 7)<=.05))[[2]]
reject.null.v1.5[2,2,1] <- prop.table(table(lapply(v1.5.ns50.smalld, '[[', 7)<=.05))[[2]]
reject.null.v1.5[3,2,1] <- prop.table(table(lapply(v1.5.ns100.smalld, '[[', 7)<=.05))[[2]]
reject.null.v1.5[4,2,1] <- prop.table(table(lapply(v1.5.ns200.smalld, '[[', 7)<=.05))[[2]]
reject.null.v1.5[1,2,2] <- prop.table(table(lapply(v1.5.ns20.smalld, '[[', 11)<=.05))[[2]]
reject.null.v1.5[2,2,2] <- prop.table(table(lapply(v1.5.ns50.smalld, '[[', 11)<=.05))[[2]]
reject.null.v1.5[3,2,2] <- prop.table(table(lapply(v1.5.ns100.smalld, '[[', 11)<=.05))[[2]]
reject.null.v1.5[4,2,2] <- prop.table(table(lapply(v1.5.ns200.smalld, '[[', 11)<=.05))[[2]]

# medium effect
reject.null.v1.5[1,3,1] <- prop.table(table(lapply(v1.5.ns20.midd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[2,3,1] <- prop.table(table(lapply(v1.5.ns50.midd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[3,3,1] <- prop.table(table(lapply(v1.5.ns100.midd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[4,3,1] <- prop.table(table(lapply(v1.5.ns200.midd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[1,3,2] <- prop.table(table(lapply(v1.5.ns20.midd, '[[', 11)<=.05))[[2]]
reject.null.v1.5[2,3,2] <- prop.table(table(lapply(v1.5.ns50.midd, '[[', 11)<=.05))[[2]]
reject.null.v1.5[3,3,2] <- prop.table(table(lapply(v1.5.ns100.midd, '[[', 11)<=.05))[[2]]
reject.null.v1.5[4,3,2] <- prop.table(table(lapply(v1.5.ns200.midd, '[[', 11)<=.05))[[2]]

# large effect
reject.null.v1.5[1,4,1] <- prop.table(table(lapply(v1.5.ns20.bigd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[2,4,1] <- prop.table(table(lapply(v1.5.ns50.bigd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[3,4,1] <- prop.table(table(lapply(v1.5.ns100.bigd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[4,4,1] <- 1 #prop.table(table(lapply(v1.5.ns200.bigd, '[[', 7)<=.05))[[2]]
reject.null.v1.5[1,4,2] <- prop.table(table(lapply(v1.5.ns20.bigd, '[[', 11)<=.05))[[2]]
reject.null.v1.5[2,4,2] <- prop.table(table(lapply(v1.5.ns50.bigd, '[[', 11)<=.05))[[2]]
reject.null.v1.5[3,4,2] <- prop.table(table(lapply(v1.5.ns100.bigd, '[[', 11)<=.05))[[2]]
reject.null.v1.5[4,4,2] <- 1 #prop.table(table(lapply(v1.5.ns200.bigd, '[[', 11)<=.05))[[2]]


##### 1.5 times different variances coverage rate #####

## null true
cov.ns20.nullT.v1.5 <- coverage(0,v1.5.ns20.nullT)
cov.ns50.nullT.v1.5 <- coverage(0,v1.5.ns50.nullT)
cov.ns100.nullT.v1.5 <- coverage(0,v1.5.ns100.nullT)
cov.ns200.nullT.v1.5 <- coverage(0,v1.5.ns200.nullT)

## small effect
t.ihat.smalld.v1.5 <- true.ihat(means=c(6,6.316), contrast=c(-1,1))
cov.ns20.smalld.v1.5 <- coverage(t.ihat.smalld.v1.5,v1.5.ns20.smalld)
cov.ns50.smalld.v1.5 <- coverage(t.ihat.smalld.v1.5,v1.5.ns50.smalld)
cov.ns100.smalld.v1.5 <- coverage(t.ihat.smalld.v1.5,v1.5.ns100.smalld)
cov.ns200.smalld.v1.5 <- coverage(t.ihat.smalld.v1.5,v1.5.ns200.smalld)

# medium effect
t.ihat.midd.v1.5 <- true.ihat(means=c(6,6.791), contrast=c(-1,1))
cov.ns20.midd.v1.5 <- coverage(t.ihat.midd.v1.5,v1.5.ns20.midd)
cov.ns50.midd.v1.5 <- coverage(t.ihat.midd.v1.5,v1.5.ns50.midd)
cov.ns100.midd.v1.5 <- coverage(t.ihat.midd.v1.5,v1.5.ns100.midd)
cov.ns200.midd.v1.5 <- coverage(t.ihat.midd.v1.5,v1.5.ns200.midd)

# big effect
t.ihat.bigd.v1.5 <- true.ihat(means=c(6,7.265), contrast=c(-1,1))
cov.ns20.bigd.v1.5 <- coverage(t.ihat.bigd.v1.5,v1.5.ns20.bigd)
cov.ns50.bigd.v1.5 <- coverage(t.ihat.bigd.v1.5,v1.5.ns50.bigd)
cov.ns100.bigd.v1.5 <- coverage(t.ihat.bigd.v1.5,v1.5.ns100.bigd)
cov.ns200.bigd.v1.5 <- coverage(t.ihat.bigd.v1.5,v1.5.ns200.bigd)


## store observ1.5d coverage rates
obs.coverage.v1.5 <- array(rep(NA,32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

obs.coverage.v1.5[,1,1] <- c(cov.ns20.nullT.v1.5[[1]],cov.ns50.nullT.v1.5[[1]],cov.ns100.nullT.v1.5[[1]],cov.ns200.nullT.v1.5[[1]])
obs.coverage.v1.5[,2,1] <- c(cov.ns20.smalld.v1.5[[1]],cov.ns50.smalld.v1.5[[1]],cov.ns100.smalld.v1.5[[1]],cov.ns200.smalld.v1.5[[1]])
obs.coverage.v1.5[,3,1] <- c(cov.ns20.midd.v1.5[[1]],cov.ns50.midd.v1.5[[1]],cov.ns100.midd.v1.5[[1]],cov.ns200.midd.v1.5[[1]])
obs.coverage.v1.5[,4,1] <- c(cov.ns20.bigd.v1.5[[1]],cov.ns50.bigd.v1.5[[1]],cov.ns100.bigd.v1.5[[1]],cov.ns200.bigd.v1.5[[1]])
obs.coverage.v1.5[,1,2] <- c(cov.ns20.nullT.v1.5[[2]],cov.ns50.nullT.v1.5[[2]],cov.ns100.nullT.v1.5[[2]],cov.ns200.nullT.v1.5[[2]])
obs.coverage.v1.5[,2,2] <- c(cov.ns20.smalld.v1.5[[2]],cov.ns50.smalld.v1.5[[2]],cov.ns100.smalld.v1.5[[2]],cov.ns200.smalld.v1.5[[2]])
obs.coverage.v1.5[,3,2] <- c(cov.ns20.midd.v1.5[[2]],cov.ns50.midd.v1.5[[2]],cov.ns100.midd.v1.5[[2]],cov.ns200.midd.v1.5[[2]])
obs.coverage.v1.5[,4,2] <- c(cov.ns20.bigd.v1.5[[2]],cov.ns50.bigd.v1.5[[2]],cov.ns100.bigd.v1.5[[2]],cov.ns200.bigd.v1.5[[2]])


##### Save tables from 1.5 different variances simulations #####
save(obs.coverage.v1.5, reject.null.v1.5, file='/users/joshwondra/R-projects/Welch rule/v1andhalfSeed2184Tables.Rdata')
load('/users/joshwondra/R-projects/Welch rule/v1andhalfSeed2184Tables.Rdata')