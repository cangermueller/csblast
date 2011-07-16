opt.args <- commandArgs(TRUE)
if (length(opt.args) < 5) {
  print("plot_lines.R OUTFILE TITLE XLAB YLAB FILE+")
  quit(status=1)
}
opt.outfile <- opt.args[1]
opt.title <- opt.args[2]
opt.xlab <- opt.args[3]
opt.ylab <- opt.args[4]
opt.files <- opt.args[5:length(opt.args)]
opt.files.n <- length(opt.files)

opt.size <- c(8, 8)
opt.legend.cex <- 1.0
opt.linewidth <- 1.8
opt.colors <- rainbow(opt.files.n + 1)
opt.colors <- opt.colors[-2]
opt.linetypes <- seq(1, opt.files.n + 1)
opt.linetypes <- opt.linetypes[-3]
opt.pch <- rep(18, opt.files.n)

data <- list(read.table(opt.files[1]))
data.n <- opt.files.n
data.x <- c()
data.y <- c()
data.labels <- c()
for (i in 1:opt.files.n) {
  data[[i]] <- read.table(opt.files[i])
  data.x <- c(data.x, data[[i]][,1])
  data.y <- c(data.y, data[[i]][,2])
  data.labels[i] <- basename(opt.files[i])
}
data.x <- unique(data.x)
data.y <- unique(data.y)

pdf(opt.outfile, title=opt.title, opt.size[1], opt.size[2])

plot(range(data.x), range(data.y), type="n", xaxt="n", yaxt="n", xlab=opt.xlab, ylab=opt.ylab)
for (i in 1:data.n) {
  lines(data[[i]][,1], data[[i]][,2], type="b", lty=opt.linetypes[i],
    lwd=opt.linewidth, col=opt.colors[i], pch=opt.pch[i])
}
axis(1, data.x)
axis(2, pretty(data.y))

legend("topright", hor=F, cex=opt.legend.cex,
	data.labels, col=opt.colors, lty=opt.linetypes, lwd=opt.linewidth)
