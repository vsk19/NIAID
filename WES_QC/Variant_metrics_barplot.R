library(ggplot2)
wd = getwd()

split = strsplit(wd,"/")[[1]][[5]]

png(paste0(split,"_Het2Hom_mqc.png"),width = 12,
	height=6,units="in",res=1200)
metrics<-read.table(paste0(split,".","variant_calling_detail_metrics"),header = TRUE,skip = 6)
ggplot(data=metrics, aes(x=SAMPLE_ALIAS, y=HET_HOMVAR_RATIO)) + geom_bar(aes(color=HET_HOMVAR_RATIO<1.5),stat="identity",fill="steelblue") + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=5))
dev.off()
