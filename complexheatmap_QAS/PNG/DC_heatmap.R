#安装需要的R包
pkgs <- rownames(installed.packages())
#本程序所需包
packages_1=c('data.table','ComplexHeatmap','circlize','RColorBrewer')
#将要安装的R包
packages_2=setdiff(packages_1,pkgs) 
#安装所需R包 
if(length(packages_2)>0){
  install.packages(packages_2, repos="http://cran.rstudio.com/")	
  source("http://bioconductor.org/biocLite.R")
  biocLite(packages_2)
}

args <- commandArgs(TRUE)

library(data.table)
mut_exp <- fread(args[1], sep=',')

mut_exp <- as.data.frame(mut_exp)   #fread读入后是data.table格式需要先转换


#通过列名匹配，来区分不同组学数据
samplenames <- colnames(mut_exp)
mutation <- agrep("mutation", samplenames, max = list(sub = 0))
mRNAseq_RSEM <- agrep("uncv2_mrnaseq_rsem",samplenames, max = list(sub = 0))
rm(samplenames)

#单独提取symbol，并与不同组学数据进行结合，拆分数据为不同组学数据表格形式
genesymbol <- mut_exp[,1]

mutation <- mut_exp[,mutation]
mutation <- cbind(genesymbol,mutation)
mRNAseq_RSEM <- mut_exp[,mRNAseq_RSEM]
mRNAseq_RSEM <- cbind(genesymbol,mRNAseq_RSEM)
rm(genesymbol,mut_exp)


#规范行名
rownames(mutation) <- mutation[,1]
mutation <- mutation[,-1]
rownames(mRNAseq_RSEM) <- mRNAseq_RSEM[,1]
mRNAseq_RSEM <- mRNAseq_RSEM[,-1]


#转为数值矩阵并去缺失，拆分为两个表
mutation <- as.matrix(mutation)
mRNAseq_RSEM <- as.matrix(mRNAseq_RSEM)

suppressWarnings(storage.mode(mutation) <- "numeric")
mutation <- na.omit(mutation)
suppressWarnings(storage.mode(mRNAseq_RSEM) <- "numeric")
mRNAseq_RSEM <- na.omit(mRNAseq_RSEM)


###读入驱动基因分析结果并进行相应处理
bionexr_result <- fread(args[2], sep=",")

genename <- agrep("genename", colnames(bionexr_result), max = list(sub = 0))
symbol <- unique(bionexr_result[,genename])

if(length(symbol) > 50){
  symbol <- symbol[1:50]
}


#通过symbol匹配原始数据表格
# mutation匹配部分
dat_match <- match(symbol,rownames(mutation))
dat_match <- na.omit(dat_match)
mutation <- mutation[dat_match,]
mutation <- as.matrix(unique(mutation))

# mRNAseq_RSEM匹配部分
dat_match <- match(symbol,rownames(mRNAseq_RSEM))
dat_match <- na.omit(dat_match)
mRNAseq_RSEM <- mRNAseq_RSEM[dat_match,]
mRNAseq_RSEM <- as.matrix(unique(mRNAseq_RSEM))

rm(dat_match,genename,symbol,bionexr_result)



###读入分型结果并进行处理
icluster_result <- fread(args[3])

#通过samplenames匹配原始数据表格
samplenames <- as.matrix(icluster_result)[,1]

# mutation匹配部分
mut_col <- colnames(mutation)
mut_col <- substr(mut_col, 1, 15)
colnames(mutation) <- mut_col
mut_col <- match(samplenames, mut_col)
mut_col <- na.omit(mut_col)
mutation <- mutation[,mut_col]

# mRNAseq_RSEM匹配部分
exp_col <- colnames(mRNAseq_RSEM)
exp_col <- substr(exp_col, 1, 15)
colnames(mRNAseq_RSEM) <- exp_col
exp_col <- match(samplenames, exp_col)
exp_col <- na.omit(exp_col)
mRNAseq_RSEM <- mRNAseq_RSEM[,exp_col]

rm(exp_col,mut_col,samplenames)


#确保两个热图行数一致
intersection <- intersect(rownames(mutation),rownames(mRNAseq_RSEM))

math <- match(intersection,rownames(mutation))
math <- na.omit(math)
mutation <- mutation[math,]

math <- match(intersection,rownames(mRNAseq_RSEM))
math <- na.omit(math)
mRNAseq_RSEM <- mRNAseq_RSEM[math,]


#对样本进行排序，根据Icluster的分类结果
colnames(icluster_result)[2] <- 'class'
ord1 <- order(icluster_result$class)
mutation <- mutation[,ord1]
mRNAseq_RSEM <- mRNAseq_RSEM[,ord1]
rm(ord1)


#热图绘制模块
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#增添bar信息
ord1 <- order(icluster_result$class)
icluster_result <- icluster_result[ord1,]
icluster_result <- as.data.frame(icluster_result)
rownames(icluster_result) <- icluster_result[,1]
icluster_result <- icluster_result[,-1]
color_num <- max(icluster_result)   #最大聚类
icluster_result <- as.character(icluster_result)
icluster_result <- as.data.frame(icluster_result)
colnames(icluster_result) <- 'class'

#对颜色进行定义
if(color_num == 2){
  colors <- list(class = c("1" =  "red", "2" = "yellow"))
}
if(color_num == 3){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green"))
}
if(color_num == 4){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green", "4"="LightGrey"))
}
if(color_num == 5){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green", "4"="LightGrey","5"="SandyBrown"))
}
if(color_num == 6){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green", "4"="LightGrey","5"="SandyBrown","6"="DarkGreen"))
}
if(color_num == 7){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green", "4"="LightGrey","5"="SandyBrown","6"="DarkGreen","7"="DeepSkyBlue"))
}
if(color_num == 8){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green", "4"="LightGrey","5"="SandyBrown","6"="DarkGreen","7"="DeepSkyBlue","8"="Magenta"))
}
if(color_num == 9){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green", "4"="LightGrey","5"="SandyBrown","6"="DarkGreen","7"="DeepSkyBlue","8"="Magenta","9"="Azure"))
}
if(color_num == 10){
  colors <- list(class = c("1" =  "red", "2" = "yellow", "3"= "green", "4"="LightGrey","5"="SandyBrown","6"="DarkGreen","7"="DeepSkyBlue","8"="Magenta","9"="Azure","10"="blue"))
}

ha1 = HeatmapAnnotation(df = icluster_result, col = colors)


#绘制多张热图
ht1 = Heatmap(mRNAseq_RSEM, cluster_columns = FALSE, 
              top_annotation = ha1,
              #col = colorRamp2(c(min(mRNAseq_RSEM),mad(mRNAseq_RSEM),max(mRNAseq_RSEM)), c("RoyalBlue", "green", "red")),
              column_names_gp = gpar(fontsize = 5),show_column_names = T, show_row_names = FALSE,
              name = "expression", column_title = "expression Heatmap")

ht2 = Heatmap(mutation,cluster_columns = FALSE,
              top_annotation = ha1,
              column_names_gp = gpar(fontsize = 5), show_column_names = T,
              name = "mutation", column_title = "mutation Heatmap")

png(args[4],width = 850, height = 550, pointsize = 12)
ht1 + ht2
dev.off()

