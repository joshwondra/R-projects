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
    sims <- list()
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
    
    return(sims)
}



##### Simulations #####

set.seed(2184)

## Simulations with equal variances, varying Ns

# true null
ev.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ev.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ev.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ev.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,2), contrast=c(-1,1))

# small effect
ev.ns20.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ev.ns50.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ev.ns100.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ev.ns200.smalld <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))

# medium effect
ev.ns20.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ev.ns50.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ev.ns100.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ev.ns200.midd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))

# large effect
ev.ns20.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ev.ns50.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ev.ns100.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ev.ns200.bigd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))


## Save simulations
save(ev.ns20.nullT,ev.ns50.nullT,ev.ns100.nullT,ev.ns200.nullT,ev.ns20.smalld,ev.ns50.smalld,ev.ns100.smalld,ev.ns200.smalld,ev.ns20.midd,ev.ns50.midd,ev.ns100.midd,ev.ns200.midd, ev.ns20.bigd,ev.ns50.bigd,ev.ns100.bigd,ev.ns200.bigd, file='/users/joshwondra/R-projects/equalvarsSeed2184.Rdata')
load('/users/joshwondra/R-projects/equalvarsSeed2184.Rdata')


##### FALSE POSITIVES AND POWER #####

## array to store proportion of rejected null hypotheses
reject.null <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))
rownames(reject.null) <- c("N=20", "N=50", "N=100", "N=200")
colnames(reject.null) <- c("d=0", "d=.2", "d=.5", "d=.8")

# true null
reject.null[1,1,1] <- prop.table(table(lapply(ev.ns20.nullT, '[[', 7)<=.05))[[2]]
reject.null[2,1,1] <- prop.table(table(lapply(ev.ns50.nullT, '[[', 7)<=.05))[[2]]
reject.null[3,1,1] <- prop.table(table(lapply(ev.ns100.nullT, '[[', 7)<=.05))[[2]]
reject.null[4,1,1] <- prop.table(table(lapply(ev.ns200.nullT, '[[', 7)<=.05))[[2]]
reject.null[1,1,2] <- prop.table(table(lapply(ev.ns20.nullT, '[[', 11)<=.05))[[2]]
reject.null[2,1,2] <- prop.table(table(lapply(ev.ns50.nullT, '[[', 11)<=.05))[[2]]
reject.null[3,1,2] <- prop.table(table(lapply(ev.ns100.nullT, '[[', 11)<=.05))[[2]]
reject.null[4,1,2] <- prop.table(table(lapply(ev.ns200.nullT, '[[', 11)<=.05))[[2]]

# small effect
reject.null[1,2,1] <- prop.table(table(lapply(ev.ns20.smalld, '[[', 7)<=.05))[[2]]
reject.null[2,2,1] <- prop.table(table(lapply(ev.ns50.smalld, '[[', 7)<=.05))[[2]]
reject.null[3,2,1] <- prop.table(table(lapply(ev.ns100.smalld, '[[', 7)<=.05))[[2]]
reject.null[4,2,1] <- prop.table(table(lapply(ev.ns200.smalld, '[[', 7)<=.05))[[2]]
reject.null[1,2,2] <- prop.table(table(lapply(ev.ns20.smalld, '[[', 11)<=.05))[[2]]
reject.null[2,2,2] <- prop.table(table(lapply(ev.ns50.smalld, '[[', 11)<=.05))[[2]]
reject.null[3,2,2] <- prop.table(table(lapply(ev.ns100.smalld, '[[', 11)<=.05))[[2]]
reject.null[4,2,2] <- prop.table(table(lapply(ev.ns200.smalld, '[[', 11)<=.05))[[2]]

# medium effect
reject.null[1,3,1] <- prop.table(table(lapply(ev.ns20.midd, '[[', 7)<=.05))[[2]]
reject.null[2,3,1] <- prop.table(table(lapply(ev.ns50.midd, '[[', 7)<=.05))[[2]]
reject.null[3,3,1] <- prop.table(table(lapply(ev.ns100.midd, '[[', 7)<=.05))[[2]]
reject.null[4,3,1] <- prop.table(table(lapply(ev.ns200.midd, '[[', 7)<=.05))[[2]]
reject.null[1,3,2] <- prop.table(table(lapply(ev.ns20.midd, '[[', 11)<=.05))[[2]]
reject.null[2,3,2] <- prop.table(table(lapply(ev.ns50.midd, '[[', 11)<=.05))[[2]]
reject.null[3,3,2] <- prop.table(table(lapply(ev.ns100.midd, '[[', 11)<=.05))[[2]]
reject.null[4,3,2] <- prop.table(table(lapply(ev.ns200.midd, '[[', 11)<=.05))[[2]]

# large effect
reject.null[1,4,1] <- prop.table(table(lapply(ev.ns20.bigd, '[[', 7)<=.05))[[2]]
reject.null[2,4,1] <- prop.table(table(lapply(ev.ns50.bigd, '[[', 7)<=.05))[[2]]
reject.null[3,4,1] <- prop.table(table(lapply(ev.ns100.bigd, '[[', 7)<=.05))[[2]]
reject.null[4,4,1] <- 1 #prop.table(table(lapply(ev.ns200.bigd, '[[', 7)<=.05))[[2]]
reject.null[1,4,2] <- prop.table(table(lapply(ev.ns20.bigd, '[[', 11)<=.05))[[2]]
reject.null[2,4,2] <- prop.table(table(lapply(ev.ns50.bigd, '[[', 11)<=.05))[[2]]
reject.null[3,4,2] <- prop.table(table(lapply(ev.ns100.bigd, '[[', 11)<=.05))[[2]]
reject.null[4,4,2] <- 1 #prop.table(table(lapply(ev.ns200.bigd, '[[', 11)<=.05))[[2]]

diff.eqvar <- reject.null[,,1]-reject.null[,,2]

## expected false positive rate and power

exp.rejects <- array(rep(NA, 32), dim=c(4,4,1), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8")))
exp.rejects[,1,] <- rep(.05,4)
exp.rejects[,2,] <- c(.095, .168, .291, .514)
exp.rejects[,3,] <- c(.338, .697, .940, 1)
exp.rejects[,4,] <- c(.693, .977, 1, 1)

reject.null-array(rep(exp.rejects,2),dim=c(4,4,2))

##### COVERAGE RATE #####

true.ihat <- function(means, contrast){
    ihat <- contrast %*% means
    return(ihat)
}

coverage <- function(true.ihat, obs.ihat, df, st.error){
    
    df <- as.matrix(df)
    t <- apply(df, 1, function(x){qt(.025, df=x)})               
    lb <- obs.ihat-t*st.error
    ub <- obs.ihat+t*st.error
    coverage <- (ub-true.ihat)*(true.ihat-lb)>0
    result <- prop.table(table(coverage))[[2]]
    return(result)
}


## null true
t.ihat.nullT <- 0

# Ns = 20
obs.ihat.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 3))
classic.df.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 6))
welch.df.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 10))
classic.ses.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 4))
welch.ses.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 8))
cov.classic.ns20.nullT <- coverage(t.ihat.nullT, obs.ihat.ns20.nullT, classic.df.ns20.nullT, classic.ses.ns20.nullT)
cov.welch.ns20.nullT <- coverage(t.ihat.nullT, obs.ihat.ns20.nullT, welch.df.ns20.nullT, welch.ses.ns20.nullT)

# Ns = 50
obs.ihat.ns50.nullT <- data.frame(lapply(ev.ns50.nullT, '[[', 3))
classic.df.ns50.nullT <- data.frame(lapply(ev.ns50.nullT, '[[', 6))
welch.df.ns50.nullT <- data.frame(lapply(ev.ns50.nullT, '[[', 10))
classic.ses.ns50.nullT <- data.frame(lapply(ev.ns50.nullT, '[[', 4))
welch.ses.ns50.nullT <- data.frame(lapply(ev.ns50.nullT, '[[', 8))
cov.classic.ns50.nullT <- coverage(t.ihat.nullT, obs.ihat.ns50.nullT, classic.df.ns50.nullT, classic.ses.ns50.nullT)
cov.welch.ns50.nullT <- coverage(t.ihat.nullT, obs.ihat.ns50.nullT, welch.df.ns50.nullT, welch.ses.ns50.nullT)

# Ns = 100
obs.ihat.ns100.nullT <- data.frame(lapply(ev.ns100.nullT, '[[', 3))
classic.df.ns100.nullT <- data.frame(lapply(ev.ns100.nullT, '[[', 6))
welch.df.ns100.nullT <- data.frame(lapply(ev.ns100.nullT, '[[', 10))
classic.ses.ns100.nullT <- data.frame(lapply(ev.ns100.nullT, '[[', 4))
welch.ses.ns100.nullT <- data.frame(lapply(ev.ns100.nullT, '[[', 8))
cov.classic.ns100.nullT <- coverage(t.ihat.nullT, obs.ihat.ns100.nullT, classic.df.ns100.nullT, classic.ses.ns100.nullT)
cov.welch.ns100.nullT <- coverage(t.ihat.nullT, obs.ihat.ns100.nullT, welch.df.ns100.nullT, welch.ses.ns100.nullT)

# Ns = 200
obs.ihat.ns200.nullT <- data.frame(lapply(ev.ns200.nullT, '[[', 3))
classic.df.ns200.nullT <- data.frame(lapply(ev.ns200.nullT, '[[', 6))
welch.df.ns200.nullT <- data.frame(lapply(ev.ns200.nullT, '[[', 10))
classic.ses.ns200.nullT <- data.frame(lapply(ev.ns200.nullT, '[[', 4))
welch.ses.ns200.nullT <- data.frame(lapply(ev.ns200.nullT, '[[', 8))
cov.classic.ns200.nullT <- coverage(t.ihat.nullT, obs.ihat.ns200.nullT, classic.df.ns200.nullT, classic.ses.ns200.nullT)
cov.welch.ns200.nullT <- coverage(t.ihat.nullT, obs.ihat.ns200.nullT, welch.df.ns200.nullT, welch.ses.ns200.nullT)


## small effect
t.ihat.smalld <- true.ihat(means=c(6,6.283), contrast=c(-1,1))

# Ns = 20
obs.ihat.ns20.smalld <- data.frame(lapply(ev.ns20.smalld, '[[', 3))
classic.df.ns20.smalld <- data.frame(lapply(ev.ns20.smalld, '[[', 6))
welch.df.ns20.smalld <- data.frame(lapply(ev.ns20.smalld, '[[', 10))
classic.ses.ns20.smalld <- data.frame(lapply(ev.ns20.smalld, '[[', 4))
welch.ses.ns20.smalld <- data.frame(lapply(ev.ns20.smalld, '[[', 8))
cov.classic.ns20.smalld <- coverage(t.ihat.smalld, obs.ihat.ns20.smalld, classic.df.ns20.smalld, classic.ses.ns20.smalld)
cov.welch.ns20.smalld <- coverage(t.ihat.smalld, obs.ihat.ns20.smalld, welch.df.ns20.smalld, welch.ses.ns20.smalld)

# Ns = 50
obs.ihat.ns50.smalld <- data.frame(lapply(ev.ns50.smalld, '[[', 3))
classic.df.ns50.smalld <- data.frame(lapply(ev.ns50.smalld, '[[', 6))
welch.df.ns50.smalld <- data.frame(lapply(ev.ns50.smalld, '[[', 10))
classic.ses.ns50.smalld <- data.frame(lapply(ev.ns50.smalld, '[[', 4))
welch.ses.ns50.smalld <- data.frame(lapply(ev.ns50.smalld, '[[', 8))
cov.classic.ns50.smalld <- coverage(t.ihat.smalld, obs.ihat.ns50.smalld, classic.df.ns50.smalld, classic.ses.ns50.smalld)
cov.welch.ns50.smalld <- coverage(t.ihat.smalld, obs.ihat.ns50.smalld, welch.df.ns50.smalld, welch.ses.ns50.smalld)

# Ns = 100
obs.ihat.ns100.smalld <- data.frame(lapply(ev.ns100.smalld, '[[', 3))
classic.df.ns100.smalld <- data.frame(lapply(ev.ns100.smalld, '[[', 6))
welch.df.ns100.smalld <- data.frame(lapply(ev.ns100.smalld, '[[', 10))
classic.ses.ns100.smalld <- data.frame(lapply(ev.ns100.smalld, '[[', 4))
welch.ses.ns100.smalld <- data.frame(lapply(ev.ns100.smalld, '[[', 8))
cov.classic.ns100.smalld <- coverage(t.ihat.smalld, obs.ihat.ns100.smalld, classic.df.ns100.smalld, classic.ses.ns100.smalld)
cov.welch.ns100.smalld <- coverage(t.ihat.smalld, obs.ihat.ns100.smalld, welch.df.ns100.smalld, welch.ses.ns100.smalld)

# Ns = 200
obs.ihat.ns200.smalld <- data.frame(lapply(ev.ns200.smalld, '[[', 3))
classic.df.ns200.smalld <- data.frame(lapply(ev.ns200.smalld, '[[', 6))
welch.df.ns200.smalld <- data.frame(lapply(ev.ns200.smalld, '[[', 10))
classic.ses.ns200.smalld <- data.frame(lapply(ev.ns200.smalld, '[[', 4))
welch.ses.ns200.smalld <- data.frame(lapply(ev.ns200.smalld, '[[', 8))
cov.classic.ns200.smalld <- coverage(t.ihat.smalld, obs.ihat.ns200.smalld, classic.df.ns200.smalld, classic.ses.ns200.smalld)
cov.welch.ns200.smalld <- coverage(t.ihat.smalld, obs.ihat.ns200.smalld, welch.df.ns200.smalld, welch.ses.ns200.smalld)


# medium effect
t.ihat.midd <- true.ihat(means=c(6,6.707), contrast=c(-1,1))

# Ns = 20
obs.ihat.ns20.midd <- data.frame(lapply(ev.ns20.midd, '[[', 3))
classic.df.ns20.midd <- data.frame(lapply(ev.ns20.midd, '[[', 6))
welch.df.ns20.midd <- data.frame(lapply(ev.ns20.midd, '[[', 10))
classic.ses.ns20.midd <- data.frame(lapply(ev.ns20.midd, '[[', 4))
welch.ses.ns20.midd <- data.frame(lapply(ev.ns20.midd, '[[', 8))
cov.classic.ns20.midd <- coverage(t.ihat.midd, obs.ihat.ns20.midd, classic.df.ns20.midd, classic.ses.ns20.midd)
cov.welch.ns20.midd <- coverage(t.ihat.midd, obs.ihat.ns20.midd, welch.df.ns20.midd, welch.ses.ns20.midd)

# Ns = 50
obs.ihat.ns50.midd <- data.frame(lapply(ev.ns50.midd, '[[', 3))
classic.df.ns50.midd <- data.frame(lapply(ev.ns50.midd, '[[', 6))
welch.df.ns50.midd <- data.frame(lapply(ev.ns50.midd, '[[', 10))
classic.ses.ns50.midd <- data.frame(lapply(ev.ns50.midd, '[[', 4))
welch.ses.ns50.midd <- data.frame(lapply(ev.ns50.midd, '[[', 8))
cov.classic.ns50.midd <- coverage(t.ihat.midd, obs.ihat.ns50.midd, classic.df.ns50.midd, classic.ses.ns50.midd)
cov.welch.ns50.midd <- coverage(t.ihat.midd, obs.ihat.ns50.midd, welch.df.ns50.midd, welch.ses.ns50.midd)

# Ns = 100
obs.ihat.ns100.midd <- data.frame(lapply(ev.ns100.midd, '[[', 3))
classic.df.ns100.midd <- data.frame(lapply(ev.ns100.midd, '[[', 6))
welch.df.ns100.midd <- data.frame(lapply(ev.ns100.midd, '[[', 10))
classic.ses.ns100.midd <- data.frame(lapply(ev.ns100.midd, '[[', 4))
welch.ses.ns100.midd <- data.frame(lapply(ev.ns100.midd, '[[', 8))
cov.classic.ns100.midd <- coverage(t.ihat.midd, obs.ihat.ns100.midd, classic.df.ns100.midd, classic.ses.ns100.midd)
cov.welch.ns100.midd <- coverage(t.ihat.midd, obs.ihat.ns100.midd, welch.df.ns100.midd, welch.ses.ns100.midd)

# Ns = 200
obs.ihat.ns200.midd <- data.frame(lapply(ev.ns200.midd, '[[', 3))
classic.df.ns200.midd <- data.frame(lapply(ev.ns200.midd, '[[', 6))
welch.df.ns200.midd <- data.frame(lapply(ev.ns200.midd, '[[', 10))
classic.ses.ns200.midd <- data.frame(lapply(ev.ns200.midd, '[[', 4))
welch.ses.ns200.midd <- data.frame(lapply(ev.ns200.midd, '[[', 8))
cov.classic.ns200.midd <- coverage(t.ihat.midd, obs.ihat.ns200.midd, classic.df.ns200.midd, classic.ses.ns200.midd)
cov.welch.ns200.midd <- coverage(t.ihat.midd, obs.ihat.ns200.midd, welch.df.ns200.midd, welch.ses.ns200.midd)


# big effect
t.ihat.bigd <- true.ihat(means=c(6,7.131), contrast=c(-1,1))

# Ns = 20
obs.ihat.ns20.bigd <- data.frame(lapply(ev.ns20.bigd, '[[', 3))
classic.df.ns20.bigd <- data.frame(lapply(ev.ns20.bigd, '[[', 6))
welch.df.ns20.bigd <- data.frame(lapply(ev.ns20.bigd, '[[', 10))
classic.ses.ns20.bigd <- data.frame(lapply(ev.ns20.bigd, '[[', 4))
welch.ses.ns20.bigd <- data.frame(lapply(ev.ns20.bigd, '[[', 8))
cov.classic.ns20.bigd <- coverage(t.ihat.bigd, obs.ihat.ns20.bigd, classic.df.ns20.bigd, classic.ses.ns20.bigd)
cov.welch.ns20.bigd <- coverage(t.ihat.bigd, obs.ihat.ns20.bigd, welch.df.ns20.bigd, welch.ses.ns20.bigd)

# Ns = 50
obs.ihat.ns50.bigd <- data.frame(lapply(ev.ns50.bigd, '[[', 3))
classic.df.ns50.bigd <- data.frame(lapply(ev.ns50.bigd, '[[', 6))
welch.df.ns50.bigd <- data.frame(lapply(ev.ns50.bigd, '[[', 10))
classic.ses.ns50.bigd <- data.frame(lapply(ev.ns50.bigd, '[[', 4))
welch.ses.ns50.bigd <- data.frame(lapply(ev.ns50.bigd, '[[', 8))
cov.classic.ns50.bigd <- coverage(t.ihat.bigd, obs.ihat.ns50.bigd, classic.df.ns50.bigd, classic.ses.ns50.bigd)
cov.welch.ns50.bigd <- coverage(t.ihat.bigd, obs.ihat.ns50.bigd, welch.df.ns50.bigd, welch.ses.ns50.bigd)

# Ns = 100
obs.ihat.ns100.bigd <- data.frame(lapply(ev.ns100.bigd, '[[', 3))
classic.df.ns100.bigd <- data.frame(lapply(ev.ns100.bigd, '[[', 6))
welch.df.ns100.bigd <- data.frame(lapply(ev.ns100.bigd, '[[', 10))
classic.ses.ns100.bigd <- data.frame(lapply(ev.ns100.bigd, '[[', 4))
welch.ses.ns100.bigd <- data.frame(lapply(ev.ns100.bigd, '[[', 8))
cov.classic.ns100.bigd <- coverage(t.ihat.bigd, obs.ihat.ns100.bigd, classic.df.ns100.bigd, classic.ses.ns100.bigd)
cov.welch.ns100.bigd <- coverage(t.ihat.bigd, obs.ihat.ns100.bigd, welch.df.ns100.bigd, welch.ses.ns100.bigd)

# Ns = 200
obs.ihat.ns200.bigd <- data.frame(lapply(ev.ns200.bigd, '[[', 3))
classic.df.ns200.bigd <- data.frame(lapply(ev.ns200.bigd, '[[', 6))
welch.df.ns200.bigd <- data.frame(lapply(ev.ns200.bigd, '[[', 10))
classic.ses.ns200.bigd <- data.frame(lapply(ev.ns200.bigd, '[[', 4))
welch.ses.ns200.bigd <- data.frame(lapply(ev.ns200.bigd, '[[', 8))
cov.classic.ns200.bigd <- coverage(t.ihat.bigd, obs.ihat.ns200.bigd, classic.df.ns200.bigd, classic.ses.ns200.bigd)
cov.welch.ns200.bigd <- coverage(t.ihat.bigd, obs.ihat.ns200.bigd, welch.df.ns200.bigd, welch.ses.ns200.bigd)


## store observed coverage rates
obs.coverage <- exp.coverage <- array(dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))

obs.coverage[,1,1] <- c(cov.classic.ns20.nullT,cov.classic.ns50.nullT,cov.classic.ns100.nullT,cov.classic.ns200.nullT)
obs.coverage[,2,1] <- c(cov.classic.ns20.smalld,cov.classic.ns50.smalld,cov.classic.ns100.smalld,cov.classic.ns200.smalld)
obs.coverage[,3,1] <- c(cov.classic.ns20.midd,cov.classic.ns50.midd,cov.classic.ns100.midd,cov.classic.ns200.midd)
obs.coverage[,4,1] <- c(cov.classic.ns20.bigd,cov.classic.ns50.bigd,cov.classic.ns100.bigd,cov.classic.ns200.bigd)
obs.coverage[,1,2] <- c(cov.welch.ns20.nullT,cov.welch.ns50.nullT,cov.welch.ns100.nullT,cov.welch.ns200.nullT)
obs.coverage[,2,2] <- c(cov.welch.ns20.smalld,cov.welch.ns50.smalld,cov.welch.ns100.smalld,cov.welch.ns200.smalld)
obs.coverage[,3,2] <- c(cov.welch.ns20.midd,cov.welch.ns50.midd,cov.welch.ns100.midd,cov.welch.ns200.midd)
obs.coverage[,4,2] <- c(cov.welch.ns20.bigd,cov.welch.ns50.bigd,cov.welch.ns100.bigd,cov.welch.ns200.bigd)

## expected coverage rates

exp.coverage <- array(rep(.95, 64), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c('classic', 'welch')))


obs.coverage
exp.coverage
reject.null
exp.rejects


##### Save tables from simulations #####
save(obs.coverage, exp.coverage, reject.null, exp.rejects, file='/users/joshwondra/R-projects/Seed2184EqVarTables.Rdata')
load('/users/joshwondra/R-projects/Seed2184EqVarTables.Rdata')



##### Plot data in LaTeX using Sweave #####
Sweave('Welch manuscript.Rnw')
