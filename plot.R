i <- 5
system( paste("cat abc_res_m0_b1.res", i, " > tmp_all.tmpres",i, sep=""))
system( paste("cat abc_res_m*.*_b*.*.res", i, " >> tmp_all.tmpres", i, sep="") )
a <- read.table( paste("tmp_all.tmpres",i, sep=""), h=TRUE, sep="\t")
dim(a)
names(a)
combs <- combn(x=ncol(a) - 2, 2)
maxColorValue <- 1000



plotcol <- function(cols, pdfname){
    pdf(pdfname, height=5*ceiling(ncol(combs)/3), width=5*3)
    layout(mat=matrix(1:ncol(combs), byrow=TRUE, nrow=ceiling(ncol(combs)/3)))
    for(i in 1:ncol(combs)){
        plot(a[,combs[1,i]], a[,combs[2,i]], col=cols, xlab=colnames(a)[combs[1,i]], ylab=colnames(a)[combs[2,i]], xlim=c(0,1), ylim=c(0,1), pch=19, cex=1.2)
        abline(a=0, b=1)
    }
    dev.off()
}

palette <- colorRampPalette(c("yellow", "green", "orange", "red", "purple", "black"))(maxColorValue)
cols <- palette[cut(a$mig, maxColorValue)]
plotcol(cols, paste("plot_mig", i, ".pdf", sep=""))

cols <- palette[cut(a$bot, maxColorValue)]
plotcol(cols, paste("plot_bot", i, ".pdf", sep=""))

### density plot
pdf(paste("density_plots", i, ".pdf", sep=""))
for(i in 1:(ncol(a)-2)){
    if(i > 1)par(new=TRUE)
    den <- density(a[,i])
    plot(den$x, den$y, xlab="", type='l', ylab="", xlim=c(0,1), col=i, axes=F, ylim=c(0, max(den$y*1.2)), lwd=2)
}
axis(side=1, tick=TRUE)
axis(side=2)
mtext(text="Success rate", side=1, 3, cex=1.4)
mtext(text="Density", side=2, 3, cex=1.4)
legend("topright", legend=colnames(a)[1:(ncol(a)-2)], col=1:(ncol(a)-2), pch=19)
dev.off()


## 2D density plot
require(lattice)
names(a)
x <- a$migration.rate
y <- a$bottleneck.strength
r <- a$F_ST
grid <- expand.grid(x=x, y=y)
grid$z <- r
pdf("test.pdf")
levelplot(F_ST~bottleneck.strength*migration.rate, a, cuts=50, scales=list(log="e"), xlab="",
          ylab="", main="Weird Function", sub="with log scales",
          colorkey = T, region = TRUE) # try cm.colors() or terrain.colors()
dev.off()
x <- seq(pi/4, 5 * pi, length = 100)
y <- seq(pi/4, 5 * pi, length = 100)
r <- as.vector(sqrt(outer(x^2, y^2, "+")))
grid <- expand.grid(x=x, y=y)
grid$z <- cos(r^2) * exp(-r/(pi^3))
levelplot(z~x*y, grid, cuts = 50, scales=list(log="e"), xlab="",
          ylab="", main="Weird Function", sub="with log scales",
          colorkey = FALSE, region = TRUE)


a$mig
a$bot
