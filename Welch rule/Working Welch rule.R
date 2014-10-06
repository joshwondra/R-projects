# Next step: extend to more than two groups


#### Step 1: Create simulations that examine Type I error and power under ideal conditions with two groups ####

t.compare.old <- function(nsims, Ns, means, vars, contrast) 
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

        fit.classic <- t.test(dv ~ group, var.equal=TRUE)
        t.classic <- fit.classic$statistic
        df.classic <- fit.classic$parameter
        p.classic <- fit.classic$p.value

        fit.welch <- t.test(dv ~ group, var.equal=FALSE)
        t.welch <- fit.welch$statistic
        df.welch <- fit.welch$parameter
        p.welch <- fit.welch$p.value

        current.sim <- list(data.frame(dv, group), fit.classic, t.classic, df.classic, p.classic, fit.welch, t.welch, df.welch, p.welch)
        names(current.sim) <- c('data', 'fit.classic', 't.classic', 'df.classic', 'p.classic', 'fit.welch', 't.welch', 'df.welch', 'p.welch')
        
        sims[[i]] <- current.sim
    }
        
    names(sims)[1:length(sims)] <- paste('sim',1:length(sims),sep='')
    
    return(sims)
}

#### Step 2: Create t test for any contrast ####

y <- rnorm(60)
x.diff <- factor(c(rep(1,30), rep(2,10), rep(3,20)))
x.eq <- factor(c(rep(1,20), rep(2,20), rep(3,20)))
contrast <- c(-1, 0, 1)

m1 <- aov(y~x.diff)
m1
sqrt(contrast^2 %*% (by(y,x.diff,var)/by(y,x.diff,length)))
sqrt(1.0374*(sum(contrast^2/by(y,x.diff,length))))
vars <- by(y,x.diff,var)
Ns <- by(y,x.diff,length)
df.classic <- sum(Ns)-length(Ns)
mse <- sum(vars*(Ns-1))/df.classic

m2 <- aov(y~x.eq)
m2
sqrt(contrast^2 %*% (by(y,x.eq,var)/by(y,x.eq,length)))
sqrt(.9954*(sum(contrast^2/by(y,x.eq,length))))
vars <- by(y,x.eq,var)
Ns <- by(y,x.eq,length)
df.classic <- sum(Ns)-length(Ns)
mse <- sum(var*(Ns-1))/df.classic
fit.contrast(m2, 'x.eq', c(-1,0,1))

t.denominator <- sqrt(contrast^2 %*% (vars/ns))
test <- sqrt()

#Need to compute different t values for welch and classic
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



##### DEBUG STEP 2 #####

test <- t.compare(nsims=2, Ns=c(10,10), means=c(2,4), vars=c(1,1))

#debug t.contrast with two groups, equal sample sizes, equal variances
group1 <- rnorm(n=20, mean=3, sd=1)
group2 <- rnorm(n=20, mean=5, sd=1)
y <- c(group1, group2)
x <- c(rep(1,20), rep(2,20))
args(t.contrast)
debug1 <- t.contrast(y,x,c(-1,1))
debug1
#classic t=4.73, df=38, p=3.08e-05, se=.36
#welch t=4.73, df=36.98, p=3.25e-05, se=.34, ihat=1.60
t.test(y~x, var.equal=TRUE) #t=4.73, df=38, p=3.08e-05
t.test(y~x, var.equal=FALSE) #t=4.73, df=36.99, p=3.25e-05
#OKAY!

#debug t.contrast with two groups, unequal sample sizes, equal variances
group1 <- rnorm(n=40, mean=3, sd=1)
group2 <- rnorm(n=20, mean=5, sd=1)
y <- c(group1, group2)
x <- c(rep(1,40), rep(2,20))
args(t.contrast)
debug2 <- t.contrast(y,x,c(-1,1))
debug2
#classic t=6.20, df=58, p=6.37e-08, se=.29
#welch t=6.04, df=35.55, p=6.49e-07, se=.30, ihat=1.78
t.test(y~x, var.equal=TRUE) #t=6.20, df=58, p=6.37e-08
t.test(y~x, var.equal=FALSE) #t=6.04, df=35.56, p=6.49e-07
#OKAY!

#debug t.contrast with two groups, unequal sample sizes, unequal variances
group1 <- rnorm(n=40, mean=3, sd=1)
group2 <- rnorm(n=20, mean=5, sd=5)
y <- c(group1, group2)
x <- c(rep(1,40), rep(2,20))
args(t.contrast)
debug3 <- t.contrast(y,x,c(-1,1))
debug3
#classic t=3.24, df=58, p=.002, se=.68
#welch t=2.36, df=20.01, p=.028, se=.94, ihat=2.22
t.test(y~x, var.equal=TRUE) #t=3.24, df=58, p=.002
t.test(y~x, var.equal=FALSE) #t=2.36, df=20.01, p=.028
#OKAY!

#debug t.contrast with three groups, unequal sample sizes, equal variances
group1 <- rnorm(20, 20, 1)
group2 <- rnorm(15, 20, 1)
group3 <- rnorm(30, 21, 1)
y <- c(group1, group2, group3)
x <- c(rep(1,20), rep(2,15), rep(3,30))
contrast1 <- c(-1,-1,2)
contrast2 <- c(-1,1,0)

library(car)
x.cont1 <- recode(x, "1:2=-1;3=2")
x.cont2 <- recode(x, "1=-1;2=1;3=0")
fit <- lm(y~x.cont1+x.cont2)
summary(fit) #df=62; contrast 1: t=6.15, p=6.11e-08; contrast 2: t=-.011, p=.991

args(t.contrast)
debug4 <- t.contrast(y, x, contrast1)
debug4 #se.classic=.49, t.classic=6.15, df.classic=62, p.classic=6.11e-08
debug5 <- t.contrast(y, x, contrast2)
debug5 #se.classic=.34, t.classic=-.011, df.classic=62, p.classic=.991
#OKAY!



#### Step 3: Integrate t test contrasts into simulation function ####

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


#### DEBUG STEP 3 ####
args(t.compare.full)

#one simulation, equal Ns, equal means, equal vars, two groups
debug5 <- t.compare.full(nsims=1, Ns=c(10,10), means=c(10,10), vars=c(2,2), contrast=c(-1,1))
debug5 
#se.classic = .71, t.classic=-.07, df.classic=18, p.classic=.94
#se.welch = .71, t.welch=-.07, df.welch=12.46, p.welch=.94
debug5.data <- debug5$sim1$data
t.test(dv ~ group, data=debug5.data, var.equal=TRUE) #t = .07, df=18, p=.94
t.test(dv ~ group, data=debug5.data, var.equal=FALSE) #t = .07, df=12.46, p=.94

#one simulation, different Ns, different means, different vars, two groups
debug6 <- t.compare.full(nsims=1, Ns=c(30,10), means=c(11,10), vars=c(3,2), contrast=c(-1,1))
debug6 
#se.classic = .44, t.classic=-1.68, df.classic=38, p.classic=.10
#se.welch = .42, t.welch=-1.76, df.welch=16.85, p.welch=.10
debug6.data <- debug6$sim1$data
t.test(dv ~ group, data=debug6.data, var.equal=TRUE) #t = 1.68, df=38, p=.10
t.test(dv ~ group, data=debug6.data, var.equal=FALSE) #t = 1.76, df=16.85, p=.10


#one simulation, different Ns, different means, different vars, three groups
debug7.1 <- t.compare.full(nsims=1, Ns=c(30,10,20), means=c(11,10,12), vars=c(3,5,2), contrast=c(-1,-1, 2))
debug7.2 <- t.compare.full(nsims=1, Ns=c(30,10,20), means=c(11,10,12), vars=c(3,5,2), contrast=c(-1, 1, 0))

debug7.1data <- debug7.1$sim1$data
debug7.2data <- debug7.2$sim1$data
debug7.1data$cont1 <- recode(debug7.1data$group, "1:2=-1;3=2")
debug7.1data$cont2 <- recode(debug7.1data$group, "1=-1;2=1;3=0")
debug7.2data$cont1 <- recode(debug7.1data$group, "1:2=-1;3=2")
debug7.2data$cont2 <- recode(debug7.1data$group, "1=-1;2=1;3=0")

debug7.1 #t.classic=4.75, df.classic=57, p.classic=1.40e-05
debug7.2 #t.classic=1.03, df.classic=57, p.classic=.31
summary(lm(dv ~ cont1+cont2, data=debug7.1data)) #cont1: t = 4.75, p = 1.4e-05, df=57
summary(lm(dv ~ cont1+cont2, data=debug7.2data)) #cont1: t = 1.03, p = .31, df=57


#test to see what happens if we have multiple contrasts
args(t.compare.full)
t.compare.full(nsims=1, Ns=c(20,10,30), means=c(11,12,10), vars=c(3,2,1), contrast=matrix(data=c(-1,-1,2,-1,1,0), nrow=3, ncol=2))

args(t.contrast)
dv <- c(rnorm(10,10), rnorm(20,9), rnorm(15,11))
groups <- c(rep(1,10),rep(2,20),rep(3,15))
cont.vec <- c(-1,-1,2)
cont.mat <- matrix(c(-1,-1,2,-1,1,0), nrow=2)

t.contrast(dv,groups,cont.vec)
t.contrast(dv,groups,cont.mat)

debug9 <- t.compare.full(nsims=1, Ns=c(20,10,30), means=c(11,12,10), vars=c(3,2,1), contrast=matrix(data=c(-1,-1,2,-1,1,0), nrow=2, byrow=TRUE, dimnames=list(c("contrast1", "contrast2"))))

debug9.data <- debug9$sim1$data
debug9.data$cont1 <- recode(debug9.data$group, "1:2=-1;3=2")
debug9.data$cont2 <- recode(debug9.data$group, "1=-1;2=1;3=0")
debug9
#contrast1 t = -2.74, df=57, p=8.17e-03
#contrast2 t = 4.21, p=9.18e-05
summary(lm(dv ~ cont1+cont2, data=debug9.data))
#cont1 t=-2.74, p=.00817
#cont2 t=4.21, p=9.18e-05

##### Simulations #####

#Compute mean differences to obtain specific effect sizes
cohen.diff <- function(d, N1, SD1, N2, SD2) #insert the mean, sample size, and standard deviation for each group
{
    poolSD <- sqrt(((N1-1)*SD1^2+(N2-1)*SD2^2)/(N1+N2-2)) #computes the pooled standard deviation
    diff <- d*poolSD #computes difference needed
    print(diff) #displays difference needed
}

set.seed(2184)
set.seed(5614)

## store equal variances, varying Ns, with null true and false

#old versions
#reject.null <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=1k"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))
#rownames(reject.null) <- c("N=20", "N=50", "N=100", "N=1k")
#colnames(reject.null) <- c("d=0", "d=.2", "d=.5", "d=.8")

reject.null <- array(rep(NA, 32), dim=c(4,4,2), dimnames=list(c("N=20", "N=50", "N=100", "N=200"), c("d=0", "d=.2", "d=.5", "d=.8"), c("classic", "welch")))
rownames(reject.null) <- c("N=20", "N=50", "N=100", "N=200")
colnames(reject.null) <- c("d=0", "d=.2", "d=.5", "d=.8")

set.seed(2184)

# Equal variances, varying Ns, null true
ev.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ev.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ev.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
ev.ns200.nullT <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6), vars=c(2,2), contrast=c(-1,1))
#store this
reject.null[1,1,1] <- prop.table(table(lapply(ev.ns20.nullT, '[[', 7)<=.05))[[2]]
reject.null[2,1,1] <- prop.table(table(lapply(ev.ns50.nullT, '[[', 7)<=.05))[[2]]
reject.null[3,1,1] <- prop.table(table(lapply(ev.ns100.nullT, '[[', 7)<=.05))[[2]]
reject.null[4,1,1] <- prop.table(table(lapply(ev.ns200.nullT, '[[', 7)<=.05))[[2]]
reject.null[1,1,2] <- prop.table(table(lapply(ev.ns20.nullT, '[[', 11)<=.05))[[2]]
reject.null[2,1,2] <- prop.table(table(lapply(ev.ns50.nullT, '[[', 11)<=.05))[[2]]
reject.null[3,1,2] <- prop.table(table(lapply(ev.ns100.nullT, '[[', 11)<=.05))[[2]]
reject.null[4,1,2] <- prop.table(table(lapply(ev.ns200.nullT, '[[', 11)<=.05))[[2]]


# Equal variances, varying Ns, small effect
args(cohen.diff)
cohen.diff(.2, 20, sqrt(2), 20, sqrt(2))
ev.ns20.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ev.ns50.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ev.ns100.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
ev.ns200.smalld <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.283), vars=c(2,2), contrast=c(-1,1))
#store this
reject.null[1,2,1] <- prop.table(table(lapply(ev.ns20.smalld, '[[', 7)<=.05))[[2]]
reject.null[2,2,1] <- prop.table(table(lapply(ev.ns50.smalld, '[[', 7)<=.05))[[2]]
reject.null[3,2,1] <- prop.table(table(lapply(ev.ns100.smalld, '[[', 7)<=.05))[[2]]
reject.null[4,2,1] <- prop.table(table(lapply(ev.ns200.smalld, '[[', 7)<=.05))[[2]]
reject.null[1,2,2] <- prop.table(table(lapply(ev.ns20.smalld, '[[', 11)<=.05))[[2]]
reject.null[2,2,2] <- prop.table(table(lapply(ev.ns50.smalld, '[[', 11)<=.05))[[2]]
reject.null[3,2,2] <- prop.table(table(lapply(ev.ns100.smalld, '[[', 11)<=.05))[[2]]
reject.null[4,2,2] <- prop.table(table(lapply(ev.ns200.smalld, '[[', 11)<=.05))[[2]]


# Equal variances, varying Ns, medium effect
args(cohen.diff)
cohen.diff(.5, 20, sqrt(2), 20, sqrt(2))
ev.ns20.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ev.ns50.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ev.ns100.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
ev.ns200.midd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,6.707), vars=c(2,2), contrast=c(-1,1))
#store this
reject.null[1,3,1] <- prop.table(table(lapply(ev.ns20.midd, '[[', 7)<=.05))[[2]]
reject.null[2,3,1] <- prop.table(table(lapply(ev.ns50.midd, '[[', 7)<=.05))[[2]]
reject.null[3,3,1] <- prop.table(table(lapply(ev.ns100.midd, '[[', 7)<=.05))[[2]]
reject.null[4,3,1] <- prop.table(table(lapply(ev.ns200.midd, '[[', 7)<=.05))[[2]]
reject.null[1,3,2] <- prop.table(table(lapply(ev.ns20.midd, '[[', 11)<=.05))[[2]]
reject.null[2,3,2] <- prop.table(table(lapply(ev.ns50.midd, '[[', 11)<=.05))[[2]]
reject.null[3,3,2] <- prop.table(table(lapply(ev.ns100.midd, '[[', 11)<=.05))[[2]]
reject.null[4,3,2] <- prop.table(table(lapply(ev.ns200.midd, '[[', 11)<=.05))[[2]]

# Equal variances, varying Ns, large effect
args(cohen.diff)
cohen.diff(.8, 20, sqrt(2), 20, sqrt(2))
ev.ns20.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ev.ns50.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ev.ns100.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
ev.ns200.bigd <- t.compare(nsims=10000, Ns=c(200,200), means=c(6,7.131), vars=c(2,2), contrast=c(-1,1))
#store this
reject.null[1,4,1] <- prop.table(table(lapply(ev.ns20.bigd, '[[', 7)<=.05))[[2]]
reject.null[2,4,1] <- prop.table(table(lapply(ev.ns50.bigd, '[[', 7)<=.05))[[2]]
reject.null[3,4,1] <- prop.table(table(lapply(ev.ns100.bigd, '[[', 7)<=.05))[[2]]
reject.null[4,4,1] <- prop.table(table(lapply(ev.ns200.bigd, '[[', 7)<=.05))[[2]]
reject.null[1,4,2] <- prop.table(table(lapply(ev.ns20.bigd, '[[', 11)<=.05))[[2]]
reject.null[2,4,2] <- prop.table(table(lapply(ev.ns50.bigd, '[[', 11)<=.05))[[2]]
reject.null[3,4,2] <- prop.table(table(lapply(ev.ns100.bigd, '[[', 11)<=.05))[[2]]
reject.null[4,4,2] <- prop.table(table(lapply(ev.ns200.bigd, '[[', 11)<=.05))[[2]]

diff.eqvar <- reject.null[,,1]-reject.null[,,2]

## expected false positive rate and power

#old
#exp.rejects <- array(rep(NA, 32), dim=c(4,4,1), dimnames=list(c("N=20", "N=50", "N=100", "N=1k"), c("d=0", "d=.2", "d=.5", "d=.8")))

exp.rejects <- array(rep(NA, 32), dim=c(3,3,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5")))
exp.rejects[,1,] <- rep(.05,3)
exp.rejects[,2,] <- c(.095, .168, .291)
exp.rejects[,3,] <- c(.338, .697, .94)
#exp.rejects[,4,] <- c(.693, .977, 1, 1)


##### Work on coverage rate ####

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

t.ihat.smalld <- true.ihat(means=c(6,6.283), contrast=c(-1,1))
t.ihat.midd <- true.ihat(means=c(6,6.707), contrast=c(-1,1))
t.ihat.bigd <- true.ihat(means=c(6,7.131), contrast=c(-1,1))

# null true
t.ihat.nullT <- 0
obs.ihat.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 3))
classic.df.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 6))
welch.df.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 10))
classic.ses.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 4))
welch.ses.ns20.nullT <- data.frame(lapply(ev.ns20.nullT, '[[', 8))
cov.classic.ns20.nullT <- coverage(t.ihat.nullT, obs.ihat.ns20.nullT, classic.df.ns20.nullT, classic.ses.ns20.nullT)

df <- as.matrix(classic.df.ns20.nullT)

apply(df, 1, function(x){qt(.025, df=x)})
