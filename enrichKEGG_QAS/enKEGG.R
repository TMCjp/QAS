#安装需要的R包
pkgs <- rownames(installed.packages())
#本程序所需包
packages_1=c('clusterProfiler','ggplot2','org.Hs.eg.db')
#将要安装的R包
packages_2=setdiff(packages_1,pkgs) 
#安装所需R包 
if(length(packages_2)>0){
  install.packages(packages_2, repos="http://cran.rstudio.com/")  
  source("http://bioconductor.org/biocLite.R")
  biocLite(packages_2)
}

args <- commandArgs(TRUE)
#suppressMessages(library(RSvgDevice))

suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))

data <- read.table(args[1], header=TRUE, sep=",", na.strings=c('NA'))
sample_unique <- unique(data[,1])
mat <- match(sample_unique,data[,1])
data <- data[mat,]
data <- na.omit(data)

de<-data[,1]
de <- as.character(de)

symbol <- select(org.Hs.eg.db,  #数据库
                 keys=de, #我们提供的ID
                 columns="ENTREZID", #需要在数据库里取出哪些ID信息
                 keytype="SYMBOL" #我们提供的ID属于哪种类型
				 )

de <- symbol$ENTREZID
#GO enrichment reuqires  clusterProfiler package
kk <- enrichKEGG(de, organism=args[2], pvalueCutoff=as.numeric(args[3]), pAdjustMethod=args[4], qvalueCutoff=as.numeric(args[5]), readable=TRUE)

kegg_analysis<-as.data.frame(summary(kk))
kegg_analysis$Description <- gsub(pattern=",", replacement=" ", kegg_analysis$Description)
write.table(kegg_analysis,args[6],sep=",",row.names=FALSE)

#devSVG(args[8])
png(args[7],width = 850, height = 550, pointsize = 12)
barplot(kk, showCategory=10)
dev.off()


