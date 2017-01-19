#安装需要的R包
pkgs <- rownames(installed.packages())
#本程序所需包
packages_1=c('data.table','RMySQL','d3heatmap','networkD3')
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
library(RMySQL)

all_gene <- fread(args[1])
diff_gene <- fread(args[2])

symbols <- diff_gene$symbol
samples <- colnames(all_gene)[2:ncol(all_gene)]

rm(all_gene,diff_gene)



#从后台直接加载表达值数据
selected_sample <- toString(samples)
selected_sample <- gsub(', ','|',selected_sample,fixed = TRUE)

#连接数据库
conn <- dbConnect(MySQL(),
                  dbname = "questionworkflow",
                  username="root", password="lihl@1123",
                  host="111.198.139.95",port=3306)  #192.168.99.20

cancer <- args[3]
selected_matrix <- dbGetQuery(conn, paste0("SELECT * FROM ",cancer,"_uncv2_mrnaseq_rsem ", " where sample_name REGEXP '",selected_sample,"'"))

dbDisconnect(conn)
rm(conn,cancer,selected_sample)



#数据处理阶段
genename <- fread("/home/ubuntu/galaxy/tools/prodigy_test/genenames/geneSort.csv",header = FALSE)
selected_value <- as.character(selected_matrix[,2])
sample_name <- as.character(selected_matrix[,1])
rm(selected_matrix)

value_list <- alist()

for(i in 1:length(selected_value)){
  value_list[i] <- strsplit(selected_value[i],'#',fixed = TRUE)
}
rm(i)

selected_value <- as.matrix(as.data.frame(value_list))
rm(value_list)

rownames(selected_value) <- as.character(as.matrix(genename)[1,])
colnames(selected_value) <- sample_name
rm(genename,sample_name)

#转为数值型矩阵并去缺失
suppressWarnings(storage.mode(selected_value) <- "numeric")
selected_value <- na.omit(selected_value)

#symbol进行相应的匹配
math <- match(symbols,rownames(selected_value))
math <- na.omit(math)
selected_value <- selected_value[math,]


library(d3heatmap)
d3 <- d3heatmap(selected_value, scale="column", colors="Blues")

#使用networkD3的saveNetwork的函数来保存html文件
d3_random <- runif(1,0,10000000000)
d3_random <- round(d3_random)

library(networkD3)

system(paste0('mkdir ',"/var/www/QAS_D3heatmap/",d3_random,"_D3heatmap/"))
saveNetwork(d3,file=paste0("/var/www/QAS_D3heatmap/",d3_random,"_D3heatmap/","heatmap.html"),selfcontained = FALSE)
write.table(paste0("/var/www/QAS_D3heatmap/",d3_random,"_D3heatmap/","heatmap.html"),args[4],sep=',',row.names=F,col.names=F,quote=FALSE)