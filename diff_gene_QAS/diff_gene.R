#已安装的R包
pkgs <- rownames(installed.packages())
#本程序所需包
packages_1=c("edgeR","locfit")
#将要安装的R包
packages_2=setdiff(packages_1,pkgs) 
#安装所需R包 
if(length(packages_2)>0){
    install.packages(packages_2, repos = "http://mirrors.xmu.edu.cn/CRAN/")    
    source("http://bioconductor.org/biocLite.R")
	options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
    biocLite(packages_2)
}

args <- commandArgs(TRUE)
library(edgeR)
library(Cairo)
#数据读入与前期预处理
x <- read.csv(args[1])
rownames(x) <- x[,1]
x <- x[,-1]

#根据样本条码信息寻找正常组织与疾病组织
col_n <- colnames(x)
i <- nchar(col_n[1])
col_n_simple <- substr(col_n, (i-1), i)

#去除噪音信息（非癌症、非正常的样本）
if(length(which(col_n_simple != "01" & col_n_simple != "11")) != 0){
  x <- x[,-which(col_n_simple != "01" & col_n_simple != "11")]
  col_n_simple <- col_n_simple[-which(col_n_simple != "01" & col_n_simple != "11")]
} else {
  x <- x
}

#判断正常和疾病两种类型是否都有
if(length(which(col_n_simple == "11")) == 0){
  stop('The data that you selected have not normal sample to compare')
}
if(length(which(col_n_simple == "01")) == 0){
  stop('The data that you selected have not disease sample to compare')
}

#改变信息，为实验设计
col_n_simple <- gsub("01","normal",col_n_simple)
col_n_simple <- gsub("11","disease",col_n_simple)
rm(col_n,i)

#实验设计，测定差异基因
y <- DGEList(counts=x,group=col_n_simple)
y <- calcNormFactors(y)
design <- model.matrix(~col_n_simple)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
# topTags(qlf)
#将所有的差异情况生成一个矩阵
tags <- topTags(qlf, n=nrow(x), adjust.method="BH", sort.by=args[2], p.value=1)
tags_matrix <- as.data.frame(tags@.Data[1])
rm(tags)

tags_result <- tags_matrix[1:as.numeric(args[3]),]
tags_result <- cbind(symbol=rownames(tags_result),tags_result)

#存储结果
write.table(tags_result,args[4],sep=",",row.names = FALSE)

#Plot log-fold change against log-counts per million(logCPM), with DE genes(差异表达基因) highlighted
#分别标记差异和正常基因
de <- decideTestsDGE(qlf)
detags <- rownames(y)[as.logical(de)]
#绘图并添加辅助线
CairoPNG(args[5], width = 750, height = 650, units = "px", pointsize = 12)
plotSmear(qlf, de.tags=detags, main="log-fold change against log-counts per million, with DE genes highlighted")
abline(h=c(-1, 1), col="blue")
