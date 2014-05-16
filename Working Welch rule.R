# Next step: extend to more than two groups


#### Step 1: Create simulations that examine Type I error and power under ideal conditions with two groups ####

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
    df.classic <- sum(Ns)-length(Ns)
    mse <- sum(vars*(Ns+1))/df.classic
    se.ihat <- sqrt(mse*(sum(contrast^2/Ns)))
    df.welch <- 
}



##### DEBUG #####

test <- t.compare(nsims=2, Ns=c(10,10), means=c(2,4), vars=c(1,1))


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

reject.null <- array(rep(NA, 32), dim=c(3,3,2), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5"), c("classic", "welch")))
rownames(reject.null) <- c("N=20", "N=50", "N=100")
colnames(reject.null) <- c("d=0", "d=.2", "d=.5")


# Equal variances, varying Ns, null true
ev.ns20.nullT <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6), vars=c(2,2))
ev.ns50.nullT <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6), vars=c(2,2))
ev.ns100.nullT <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6), vars=c(2,2))
#ev.ns1000.nullT <- t.compare(nsims=10000, Ns=c(1000,1000), means=c(6,6), vars=c(2,2))
#store this
reject.null[1,1,1] <- prop.table(table(lapply(ev.ns20.nullT, '[[', 5)<=.05))[[2]]
reject.null[2,1,1] <- prop.table(table(lapply(ev.ns50.nullT, '[[', 5)<=.05))[[2]]
reject.null[3,1,1] <- prop.table(table(lapply(ev.ns100.nullT, '[[', 5)<=.05))[[2]]
#reject.null[4,1,1] <- prop.table(table(lapply(ev.ns1000.nullT, '[[', 5)<=.05))[[2]]
reject.null[1,1,2] <- prop.table(table(lapply(ev.ns20.nullT, '[[', 9)<=.05))[[2]]
reject.null[2,1,2] <- prop.table(table(lapply(ev.ns50.nullT, '[[', 9)<=.05))[[2]]
reject.null[3,1,2] <- prop.table(table(lapply(ev.ns100.nullT, '[[', 9)<=.05))[[2]]
#reject.null[4,1,2] <- prop.table(table(lapply(ev.ns1000.nullT, '[[', 9)<=.05))[[2]]


# Equal variances, varying Ns, small effect
args(cohen.diff)
cohen.diff(.2, 20, sqrt(2), 20, sqrt(2))
ev.ns20.smalld <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.283), vars=c(2,2))
ev.ns50.smalld <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.283), vars=c(2,2))
ev.ns100.smalld <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.283), vars=c(2,2))
#ev.ns1000.smalld <- t.compare(nsims=10000, Ns=c(1000,1000), means=c(6,6.283), vars=c(2,2))
#store this
reject.null[1,2,1] <- prop.table(table(lapply(ev.ns20.smalld, '[[', 5)<=.05))[[2]]
reject.null[2,2,1] <- prop.table(table(lapply(ev.ns50.smalld, '[[', 5)<=.05))[[2]]
reject.null[3,2,1] <- prop.table(table(lapply(ev.ns100.smalld, '[[', 5)<=.05))[[2]]
#reject.null[4,2,1] <- 1 #prop.table(table(lapply(ev.ns1000.smalld, '[[', 5)<=.05))[[2]]
reject.null[1,2,2] <- prop.table(table(lapply(ev.ns20.smalld, '[[', 9)<=.05))[[2]]
reject.null[2,2,2] <- prop.table(table(lapply(ev.ns50.smalld, '[[', 9)<=.05))[[2]]
reject.null[3,2,2] <- prop.table(table(lapply(ev.ns100.smalld, '[[', 9)<=.05))[[2]]
#reject.null[4,2,2] <- 1 #prop.table(table(lapply(ev.ns1000.smalld, '[[', 9)<=.05))[[2]]


# Equal variances, varying Ns, medium effect
args(cohen.diff)
cohen.diff(.5, 20, sqrt(2), 20, sqrt(2))
ev.ns20.midd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,6.707), vars=c(2,2))
ev.ns50.midd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,6.707), vars=c(2,2))
ev.ns100.midd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,6.707), vars=c(2,2))
#ev.ns1000.midd <- t.compare(nsims=10000, Ns=c(1000,1000), means=c(6,6.707), vars=c(2,2))
#store this
reject.null[1,3,1] <- prop.table(table(lapply(ev.ns20.midd, '[[', 5)<=.05))[[2]]
reject.null[2,3,1] <- prop.table(table(lapply(ev.ns50.midd, '[[', 5)<=.05))[[2]]
reject.null[3,3,1] <- prop.table(table(lapply(ev.ns100.midd, '[[', 5)<=.05))[[2]]
#reject.null[4,3,1] <- 1 #prop.table(table(lapply(ev.ns1000.midd, '[[', 5)<=.05))[[2]]
reject.null[1,3,2] <- prop.table(table(lapply(ev.ns20.midd, '[[', 9)<=.05))[[2]]
reject.null[2,3,2] <- prop.table(table(lapply(ev.ns50.midd, '[[', 9)<=.05))[[2]]
reject.null[3,3,2] <- prop.table(table(lapply(ev.ns100.midd, '[[', 9)<=.05))[[2]]
#reject.null[4,3,2] <- 1 #prop.table(table(lapply(ev.ns1000.midd, '[[', 9)<=.05))[[2]]

# Equal variances, varying Ns, large effect
args(cohen.diff)
cohen.diff(.8, 20, sqrt(2), 20, sqrt(2))
ev.ns20.bigd <- t.compare(nsims=10000, Ns=c(20,20), means=c(6,7.131), vars=c(2,2))
ev.ns50.bigd <- t.compare(nsims=10000, Ns=c(50,50), means=c(6,7.131), vars=c(2,2))
ev.ns100.bigd <- t.compare(nsims=10000, Ns=c(100,100), means=c(6,7.131), vars=c(2,2))
#ev.ns1000.bigd <- t.compare(nsims=10000, Ns=c(1000,1000), means=c(6,7.131), vars=c(2,2))
#store this
reject.null[1,4,1] <- prop.table(table(lapply(ev.ns20.bigd, '[[', 5)<=.05))[[2]]
reject.null[2,4,1] <- 1 #prop.table(table(lapply(ev.ns50.bigd, '[[', 5)<=.05))[[2]]
reject.null[3,4,1] <- 1 #prop.table(table(lapply(ev.ns100.bigd, '[[', 5)<=.05))[[2]]
#reject.null[4,4,1] <- 1 #prop.table(table(lapply(ev.ns1000.bigd, '[[', 5)<=.05))[[2]]
reject.null[1,4,2] <- prop.table(table(lapply(ev.ns20.bigd, '[[', 9)<=.05))[[2]]
reject.null[2,4,2] <- 1 #prop.table(table(lapply(ev.ns50.bigd, '[[', 9)<=.05))[[2]]
reject.null[3,4,2] <- 1 #prop.table(table(lapply(ev.ns100.bigd, '[[', 9)<=.05))[[2]]
#reject.null[4,4,2] <- 1 #prop.table(table(lapply(ev.ns1000.bigd, '[[', 9)<=.05))[[2]]

diff.eqvar <- reject.null[,,1]-reject.null[,,2]

## expected false positive rate and power

#old
#exp.rejects <- array(rep(NA, 32), dim=c(4,4,1), dimnames=list(c("N=20", "N=50", "N=100", "N=1k"), c("d=0", "d=.2", "d=.5", "d=.8")))

exp.rejects <- array(rep(NA, 32), dim=c(3,3,1), dimnames=list(c("N=20", "N=50", "N=100"), c("d=0", "d=.2", "d=.5")))
exp.rejects[,1,] <- rep(.05,3)
exp.rejects[,2,] <- c(.095, .168, .291)
exp.rejects[,3,] <- c(.338, .697, .94)
#exp.rejects[,4,] <- c(.693, .977, 1, 1)
