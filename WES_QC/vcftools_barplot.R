#install.packages("ggplot2")
library(ggplot2)

wd = getwd()

split = strsplit(wd,"/")[[1]][[5]]

batch = read.table(paste0(split,".","het"), header=T)

png("test.png",width= 7,
  height    = 4,
  units     = "in",
  res       = 1200,
  pointsize = 2)

ggplot(data=batch, aes(x=INDV, y=F)) + geom_bar(aes(color=F>0.1),stat="identity",fill="steelblue") + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=3)) + ylim(-0.5,0.5)

dev.off()
