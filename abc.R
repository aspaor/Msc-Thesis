rm(list=ls())
library(abc)
## number of abc replications
args <- commandArgs(trailingOnly = TRUE)

out <- "5"

doSims <- FALSE
header <- FALSE
mig1 <- 0
bot1 <- 1
index <- "XXX"

if (length(args)==1) {
    doSims <- as.logical(args[1])
}
if (length(args)==2) {
    doSims <- as.logical(args[1])
    header <- as.logical(args[2])
}
if (length(args)==5) {
    doSims <- as.logical(args[1])
    header <- as.logical(args[2])
    mig1 <- as.numeric(args[3])
    bot1 <- as.numeric(args[4])
    index <- args[5]
}

plotNam=c("PCA", "PCA Summary Statistics", "ABC", expression("F_ST"), expression("f_3"), "Prevosti Distance", "Number of shared polymorphisms", "migration rate", "bottleneck strength")


write(c(doSims, header, mig1, bot1, index), stderr())

##print(doSims)
##print(header)
##print(mig1)
##print(bot1)
##stop("check")

## first generate 100 datasets and use them as pseudoobserved data
## this part you have it already from ms-prime
samplesize  <- 30
pops <- 3
reps <- 100

num.config <- rep(samplesize/pops, pops)
config <- paste(num.config, collapse=" ")
cmd <- paste("ms ", samplesize, reps, "-I ", pops, paste(rep(samplesize/pops, pops), collapse=" "), 0, "-em 0 1 3 ", mig1, " -em 0 3 1 ",mig1, "-ej 0.1 2 3 -ej 1 3 1 -t 100 -r 1000 100 -en 0.001 2 ", bot1, " -en 0.1 2 1 > ", paste("ms.obs.", index, sep=""))
system(cmd)
## I convert the dataset to summary statistics using the msABC
cmd <- paste("CoMuStats -input ", paste("ms.obs.", index, sep=""), " -ms -npop ", pops, paste( rep(samplesize/pops, pops), collapse=" "), " -sepPops -pairwiseFst > ", paste("ms.obs.stats.", index, sep=""))
system(cmd)
obs.stats <- read.table( paste("ms.obs.stats.",index, sep=""), h=TRUE)

## get ms Data
read.ms.dat <- function(filepath, n=samplesize){
    out.list <- list()
    con = file(filepath, "r")
    cnt <- 1
    while ( TRUE ) {
        line = readLines(con, n = 1)
        if ( length(line) == 0 ) {
            break
        }
        if( grepl("[01]", x=line, perl=TRUE) &
            !grepl("[^01]", x=line,  perl=TRUE)){
            int.cnt <- ((cnt-1)%/%n)+1
            int.cnt.2 <- (cnt-1) %% n
            ##print(c(int.cnt, int.cnt.2, cnt))
            if(int.cnt.2 == 0){
                out.list[[int.cnt]] <- list()
            }
            ##out.list[[int.cnt]] <-
            out.list[[int.cnt]] <- append(out.list[[int.cnt]], strsplit(line, split=""))
            cnt <- cnt+1
        }
    }
    out.mat <- lapply(out.list, function(x){
        tmp.mat <- matrix(as.numeric(unlist(x)), nrow=n, byrow=TRUE)
    })
    close(con)
    return(out.mat)
}
ms.obs.data <- read.ms.dat(paste("ms.obs.", index, sep=""))

mig.names <- list()
ii <- 1
for(i in 1:pops){
    for(j in 1:pops){
        mig.names[[ii]] <- paste(i, j, sep=".")
        ii <- ii+1
    }
}

abcreps <- 200000
unifMin <- 0
unifMax <- 10
migMat.prior <- c()
pars <- c()
if( doSims == TRUE){
    ## Now let's create a large dataset with simulations to use it as a base dataset for the abc
    migMat.prior <- matrix(runif(abcreps*pops*pops, unifMin, unifMax), nrow=abcreps, ncol=pops*pops)
    pars <- as.data.frame(migMat.prior)
    colnames(pars) <- unlist(mig.names)
    write.table(x=pars, file="migmat.prior", quote=F, sep="\t")
    mig.cmd <- apply(migMat.prior, 1, paste, collapse=" ")
    msABCcmds <- paste("ms ", samplesize, 1, "-I ", pops, config, 0, "-ma ", mig.cmd, "-t 100")
    system(paste(msABCcmds[1], ">ms.sims"))
    sapply(2:length(msABCcmds), function(x){
        ##if(x %% 500 == 0){ print(x)}
        system(paste(msABCcmds[x], "|sed '1,2d' >> ms.sims"))
    })
    cmd.sims <- paste("CoMuStats -input ms.sims -ms -npop ", pops, paste( rep(samplesize/pops, pops), collapse=" "), " -sepPops -pairwiseFst > ms.sims.stats")
    system(cmd.sims)
}

sims.stats <- read.table("ms.sims.stats", h=TRUE)
if(doSims == FALSE){
    pars <- read.table("migmat.prior", h=TRUE)
}



mig.names <- list()
ii <- 1
for(i in 1:pops){
    for(j in 1:pops){
        mig.names[[ii]] <- paste(i, j, sep=".")
        ii <- ii+1
    }
}


f3 <- function(ms, config, o, a, b){
    inds <- rep(1:length(config), config)
    freqs <- lapply(1:ncol(ms), function(x){
        col <- ms[,x]
        sapply(1:length(config), function(i){
            myinds <- inds == i
            f <- sum(ms[myinds, x])/sum(myinds)
            return(f)
        })
    })
    mean(unlist(lapply(freqs, function(x){ (x[o]-x[a])*(x[o]-x[b]) })))
}

f3all <- function(ms, config, o, a, b){
    inds <- rep(1:length(config), config)
    freqs <- lapply(1:ncol(ms), function(x){
        col <- ms[,x]
        sapply(1:length(config), function(i){
            myinds <- inds == i
            f <- sum(ms[myinds, x])/sum(myinds)
            return(f)
        })
    })
    v1 <- mean(unlist(lapply(freqs, function(x){ (x[o]-x[a])*(x[o]-x[b]) })))
    v2 <- mean(unlist(lapply(freqs, function(x){ (x[a]-x[o])*(x[a]-x[b]) })))
    v3 <- mean(unlist(lapply(freqs, function(x){ (x[b]-x[o])*(x[b]-x[a]) })))
    return(c(v1, v2, v3))    
}


freq.stat <- function(ms, config, o, a, b){
    inds <- rep(1:length(config), config)
    freqs <- lapply(1:ncol(ms), function(x){
        col <- ms[,x]
        sapply(1:length(config), function(i){
            myinds <- inds == i
            f <- sum(ms[myinds, x])/sum(myinds)
            return(f)
        })
    })
    v1 <- mean(unlist(lapply(freqs, function(x){1 - abs(x[a] - x[b])})))
    v2 <- mean(unlist(lapply(freqs, function(x){1 - abs(x[b] - x[o])})))
    v3 <- mean(unlist(lapply(freqs, function(x){1 - abs(x[o] - x[a])})))
    return(c(v1, v2, v3))    
}

shared.stat <- function(ms, config, o, a, b){
    inds <- rep(1:length(config), config)
    freqs <- lapply(1:ncol(ms), function(x){
        col <- ms[,x]
        sapply(1:length(config), function(i){
            myinds <- inds == i
            f <- sum(ms[myinds, x])/sum(myinds)
            return(f)
        })
    })
    combs <- combn(length(config), 2)
    res <- apply(combs, 2, function(x){
        mean(unlist(lapply(freqs, function(f){ f[x[1]] * f[x[2]] > 0})))
    })
    return(res)
}

##pars <- as.data.frame(migMat.prior)
##colnames(pars) <- unlist(mig.names)
fsts <- grep(pattern="Fst", x=names(sims.stats), ignore.case=F)

res <- list()
for(i in 1:reps){
    myabc <- abc(target=as.data.frame(obs.stats[i,fsts]), param=pars, sumstat=sims.stats[,fsts], tol=0.01, method='loclinear', hcorr=TRUE)
    tmp.mat <- matrix(summary(myabc)[3,], nrow=pops, byrow=TRUE)
    uind <- upper.tri(tmp.mat)
    tmp.mat <- (tmp.mat + t(tmp.mat))/2
    res[[i]] <- tmp.mat[uind]
}
res.mat <- matrix(unlist(res), ncol=pops*(pops-1)/2, byrow=TRUE)
abc.cor <- sum(apply(res.mat, 1, which.max)  == 3)/nrow(res.mat)

getSampleDist <- function(data, config){
    nclus <- length(config)
    combs <- combn(nclus,2)
    res <- matrix(0, nrow=nclus, ncol=nclus)
    cols <- rep(1:length(config), config)
    means <- sapply(1:length(config), function(x){
        m <- data[cols==x,]
        if(class(m) == "numeric"){
            m <- matrix(m, nrow=1)
        }
        apply(m, 2, mean)
    })
    means <- t(means)
    as.numeric(dist(means, method="euclidean", diag=FALSE))
}

## Now I start with the PCA process for each one of the observed datasets
pcas <- list()
pca.dist <- list()
for(i in 1:length(ms.obs.data)){
    pcas[[i]] <- prcomp(ms.obs.data[[i]], center=TRUE, scale.=TRUE)
    pca.dist[[i]] <- getSampleDist( pcas[[i]]$x[,c(1,2)], num.config)
}
pca.dist.mat <- matrix(unlist(pca.dist), nrow=reps, byrow=TRUE)
pca.cor <- sum(apply(pca.dist.mat, 1, which.min)  == 3)/nrow(pca.dist.mat)



## Now I start with the PCA process for each one of the observed datasets
pop.stats <- intersect(grep(".*\\d$", names(obs.stats), value=TRUE), grep("^Fst", x=names(obs.stats), perl=TRUE, invert=TRUE, value=TRUE))
pcas.2 <- list()
pca.dist.2 <- list()
for(i in 1:length(ms.obs.data)){
    tmp.mat <- matrix(as.numeric(obs.stats[i,pop.stats]), nrow=length(num.config), byrow=TRUE)
    good.inds <- which(apply(tmp.mat, 2, var) > 0)
    pcas.2[[i]] <- prcomp(tmp.mat[, good.inds], center=TRUE, scale.=TRUE)
    pca.dist.2[[i]] <- getSampleDist( pcas.2[[i]]$x[,c(1,2)], rep(1, length(num.config)))
}
pca.dist.mat.2 <- matrix(unlist(pca.dist.2), nrow=reps, byrow=TRUE)
pca.cor.2 <- sum(apply(pca.dist.mat.2, 1, which.min)  == 3)/nrow(pca.dist.mat.2)

## shared.stat
shared.stat.values <- lapply(ms.obs.data, shared.stat, num.config, 1,2,3)
shared.cor.tmp <- lapply(shared.stat.values, function(x){ which.max(x) == 3} )
shared.cor <- sum(unlist(shared.cor.tmp))/length(unlist(shared.cor.tmp))
## freq stat
freq.stat.values <- lapply(ms.obs.data, freq.stat, num.config, 1,2,3)
freq.cor.tmp <- lapply(freq.stat.values, function(x){ which.max(x) == 1} )
freq.cor <- sum(unlist(freq.cor.tmp))/length(unlist(freq.cor.tmp))
## f3values
f3values <- lapply( ms.obs.data, f3all, num.config, 1,2,3)
## the greatest values should correspond to the outgroup
f3.cor.tmp <- lapply(f3values, function(x){ which.max(x) == 1} )
f3.cor <- sum(unlist(f3.cor.tmp))/length(unlist(f3.cor.tmp))
### fst based
fsts <- grep(pattern="Fst_", x=names(sims.stats), ignore.case=F)
fst.values <- list()
for(i in 1:length(ms.obs.data)){
    fst.values[[i]] <- obs.stats[i,fsts]
}
fst.mat <- matrix(unlist(fst.values), nrow=reps, byrow=TRUE)
fst.cor <- sum(apply(fst.mat, 1, which.min)  == 3)/nrow(fst.mat)

if(header==TRUE){
    mat=matrix(c(PCA=pca.cor, pca.cor.2, ABC=abc.cor, fst.cor, f3.cor, freq.cor, shared.cor, mig=mig1, bot=bot1),nrow=1)
    colnames(mat)=plotNam
    write.table(x=mat, file=paste("abc_res_m", mig1, "_b", bot1, ".res",out, sep=""), row.names=F, col.names=T, quote=F, sep="\t", append=FALSE)
}else{
    mat=matrix(c(PCA=pca.cor, pca.cor.2, ABC=abc.cor, fst.cor, f3.cor, freq.cor, shared.cor, mig=mig1, bot=bot1),nrow=1)
    colnames(mat)=plotNam
    write.table(x=mat, file=paste("abc_res_m", mig1, "_b", bot1, ".res",out, sep=""), row.names=F, col.names=F, quote=F, sep="\t", append=TRUE)
}

## now let's play a bit with ms data
## F3





## pdf("allplots.pdf")
## for(i in 1:length(ms.obs.data)){
##     cols <- rep(1:length(num.config), num.config)
##     plot(pcas[[i]]$x[,1], pcas[[i]]$x[,2], col=cols)
## }
## dev.off()

