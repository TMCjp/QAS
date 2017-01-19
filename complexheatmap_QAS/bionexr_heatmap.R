args <- commandArgs(TRUE)

#######
#读入并处理下载的突变数据
mutation <- read.table(args[1],header = TRUE,sep=",")
mutation <- mutation[,c(1,7)]

#合并标示符，以便进行计数
merge_data <- paste(mutation$Hugo_Symbol,mutation$Tumor_Sample_Barcode,sep="_")
#利用table函数将邻接表转变为计数邻接矩阵
count_data <- as.matrix(table(merge_data))
rm(merge_data)

#计数完成后进行拆分，并合并成matrix
symbol <- do.call('rbind',strsplit(rownames(count_data),"_"))
count_data <- cbind(symbol,count_data)
rownames(count_data) <- NULL
rm(symbol)

#构建最终关系矩阵
dat <- matrix(as.numeric(0),
              nrow=length(unique(mutation$Hugo_Symbol)),
              ncol=length(unique(mutation$Tumor_Sample_Barcode))
)
rownames(dat) <- unique(mutation$Hugo_Symbol)
colnames(dat) <- unique(mutation$Tumor_Sample_Barcode)
dat[as.matrix(count_data[,c(1,2)])] <- as.numeric(count_data[,3])

# dat <- cbind(symbol=rownames(dat),dat)
# rownames(dat) <- NULL

rm(count_data)
#######


#######
#读入下载的表达量数据，并进行处理
express <- read.table(args[2],header = TRUE,sep=",")

colnames(express) <- gsub(".", "-", colnames(express),fixed = TRUE)
# express <- cbind(symbol=rownames(express),express)
# rownames(express) <- NULL
#######


#######
#读入bionexr结果数据，对结果进行处理，提取出驱动基因
bionexr_result <- read.table(args[3],header=TRUE,sep=",")

genename <- agrep("name", colnames(bionexr_result), max = list(sub = 0))
symbol <- unique(bionexr_result[,genename])

if(length(symbol) > 50){
  symbol <- symbol[1:50]
}


#symbol_match <- match(symbol,mutation$Hugo_Symbol)
#result <- mutation[symbol_match,]$Hugo_Symbol

dat_match <- match(symbol,rownames(dat))
dat <- dat[dat_match,]
dat <- as.matrix(unique(dat))


express_match <- match(symbol,rownames(express))
express <- express[express_match,]
express <- as.matrix(unique(express))


#write.table(dat,args[2],sep=",",row.names=FALSE)

library(ComplexHeatmap)
#只取30个样本进行展示
ht1 <- Heatmap(dat[,1:30], name = "mutation", column_title = "mutation Heatmap")
ht2 <- Heatmap(express[,1:30], name = "expression", column_title = "expression Heatmap")

png(args[4], width = 960, height = 850, units = "px", pointsize = 12)
ht1+ht2

dat <- cbind(symbol=rownames(dat),dat)
rownames(dat) <- NULL
write.table(dat,args[5],sep=",",row.names=FALSE)

express <- cbind(symbol=rownames(express),express)
rownames(express) <- NULL
write.table(express,args[6],sep=",",row.names=FALSE)
