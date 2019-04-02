#install.packages("ggplot2")
library(ggplot2)

batch = read.table("BATCH12.het", header=T)

png("Inbreeding_Coefficient_mqc.png",res =600)

ggplot(data=batch, aes(x=INDV, y=F)) + geom_bar(aes(color=F>0.1),stat="identity",fill="steelblue") + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=3)) + ylim(-0.5,0.5)

dev.off()
