rm(list=ls())
library(foreach)
library(doMC)
registerDoMC(24) 

N <- 2000

migs <- c(runif(N, 0, 4))
bots.log <- c(runif(N, -4, 0))

system("Rscript abc.R FALSE TRUE 0 1 0 > log")
i <- 1

foreach(i=1:(N)) %dopar% {
    cmd <- paste("Rscript abc.R FALSE FALSE ", migs[i], 10^bots.log[i], i, " > log5" )
    system(cmd)    
}

