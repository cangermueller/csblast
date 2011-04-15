opt.bin_size <- 0.5
opt.res <- 300
opt.size <- c(8, 8)
opt.legend.cex <- 1.0
opt.linewidth <- 1.5
opt.pch <- 18
opt.linetype <- 1
opt.args <- commandArgs(TRUE)
if (length(opt.args) < 3) {
  print("Missing arguments!")
  print("plot.R OUT-FILE TITLE IN-FILE+")
  quit()
}
opt.outfile <- opt.args[1]
opt.title <- opt.args[2]

data.values <- list()
data.mean <- vector()
data.labels <- vector()
data.min <- 100
data.max <- 0
for (i in 1:(length(opt.args)-2)) {
  d <- read.csv(opt.args[i + 2], head = F)
  data.values[[i]] <- d[,1]
  data.mean[i] <- mean(d[,1])
  data.labels[i] <- sprintf("%s: %3.2f", basename(opt.args[i + 2]), data.mean[i])
  data.min <- min(data.min, data.values[[i]])
  data.max <- max(data.max, data.values[[i]])
}
data.n <- length(data.values)
data.breaks <- seq(floor(data.min) - 0.5 * opt.bin_size, round(data.max + 0.5) + 0.5 * opt.bin_size, by=opt.bin_size)
data.freq <- list()
data.freq.max <- 0
for (i in 1:data.n) {
  h <- hist(data.values[[i]], data.breaks, plot=F)
  data.freq[[i]] <- h$counts
  sum <- sum(data.freq[[i]])
  for (j in 1:length(data.freq[[i]])) data.freq[[i]][j] <- data.freq[[i]][j] / sum
  data.freq.max <- max(data.freq.max, data.freq[[i]])
}
data.x <- data.breaks[1:length(data.breaks)-1] + 0.5 * opt.bin_size

pdf(opt.outfile, title=opt.title, opt.size[1], opt.size[2])

xrange <- range(data.breaks)
yrange <- range(0, data.freq.max + 0.1)
plot(xrange, yrange, type="n", xlab="Neff", ylab="Frequency", xaxt="n")
plot.colors <- rainbow(data.n)
plot.linetype <- rep(opt.linetype, data.n)
plot.pch <- rep(opt.pch, data.n)
for (i in 1:data.n) {
  lines(data.x, data.freq[[i]], type="l", lwd=opt.linewidth,
      lty=plot.linetype[i], col=plot.colors[i], pch=plot.pch) 
}
axis(1, at=data.x, labels=data.x)

legend("topright", hor = FALSE, cex = opt.legend.cex,
	data.labels, col=plot.colors, pch=plot.pch)
