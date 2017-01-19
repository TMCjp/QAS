#安装需要的R包
pkgs <- rownames(installed.packages())
#本程序所需包
packages_1=c('ComplexHeatmap','circlize','colorspace','GetoptLong','RSvgDevice','data.table','RMySQL')
#将要安装的R包
packages_2=setdiff(packages_1,pkgs) 
#安装所需R包 
if(length(packages_2)>0){
  install.packages(packages_2, repos="http://cran.rstudio.com/")	
  source("http://bioconductor.org/biocLite.R")
  biocLite(packages_2)
}

suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(colorspace)))
suppressWarnings(suppressMessages(library(GetoptLong)))
suppressWarnings(suppressMessages(library(RSvgDevice)))

argv <- commandArgs(TRUE);


library(data.table)
library(RMySQL)

all_gene <- fread(argv[1])
diff_gene <- fread(argv[2])

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

cancer <- argv[3]
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



#设定各项参数（原complexheatmap代码）
cluster_rows = as.logical(argv[4])
row_hclust_side = argv[5]
row_names_side = argv[6]
cluster_columns = as.logical(argv[7])
column_names_side = argv[8]
column_hclust_side = argv[9]
km = as.numeric(argv[10])
outPath = "/var/www/heatmapoutput/"
numberRange = c(as.numeric(argv[11]),as.numeric(argv[12]),as.numeric(argv[13]))
colorRange = c(argv[14],argv[15],argv[16])
colorRange = gsub(pattern="__pd__", replacement="#", colorRange)

ha_mix_col_his_val = argv[17]
ha_mix_col_den_val = argv[18]
ha_mix_col_vio_val = argv[19]
ha_mix_col_heatmap_val = argv[20]
ha_mix_row_his_val = argv[21]
ha_mix_row_den_val = argv[22]
ha_mix_row_vio_val = argv[23]
ha_mix_row_heatmap_val = argv[24]

fileNamePrefix <- runif(1,0,1000000000)
fileNamePrefix <- round(fileNamePrefix)



#构建复杂热图
ff_1 <- hclust(dist(t(selected_value)))
selected_value <- selected_value[, ff_1$order]

kc <- kmeans(selected_value, km)
cla <- kc$cluster
clas <- unique(kc$cluster)
data <- NULL
for(i in 1:km) {
	cs <- names(cla[cla==i])
  #cs <- names(cla[cla==clas[i]])
  protein_us <- selected_value[cs, ]
  ff <- hclust(dist(protein_us))
  write.csv(as.data.frame(protein_us[ff$order, ]), paste(outPath,fileNamePrefix,"inte_", i,".csv", sep="", collapse = ""))
}


ha_mix_col_his = anno_histogram(selected_value, which = "column")                                                             

ha_mix_col_den = anno_density(selected_value, type = "line", which = "column")

ha_mix_col_vio = anno_density(selected_value, type = "violin", which = "column")

ha_mix_col_heatmap = anno_density(selected_value, type = "heatmap", which = "column")

ha_mix_row_his = HeatmapAnnotation(histogram = anno_histogram(selected_value, which = "row"),which = "row",width = unit(2, "cm"))
ha_mix_row_den = HeatmapAnnotation(density_line = anno_density(selected_value, type = "line", which = "row"),which = "row",width = unit(2, "cm"))
ha_mix_row_vio = HeatmapAnnotation(violin = anno_density(selected_value, type = "violin", which = "row"),which = "row",width = unit(2, "cm"))
ha_mix_row_heatmap = HeatmapAnnotation(heatmap = anno_density(selected_value, type = "heatmap", which = "row"), which = "row",width = unit(2, "cm"))


col_anno <- c(ha_mix_col_his_val, ha_mix_col_den_val, ha_mix_col_vio_val, ha_mix_col_heatmap_val)
row_anno <- c(ha_mix_row_his_val, ha_mix_row_den_val, ha_mix_row_vio_val, ha_mix_row_heatmap_val)


col_anno <- col_anno[which(col_anno!="NULL")]
options(expressions=50000)
devSVG(paste(outPath,fileNamePrefix,"complex.SVG",sep=""))
if(length(col_anno)==0) {
  Heatmap(selected_value, name = "foo") + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==1) {
  Heatmap(selected_value, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1])), top_annotation_height = unit(2, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==2) {
  Heatmap(selected_value, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2])), top_annotation_height = unit(4, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==3){
  Heatmap(selected_value, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2]), col_anno_3=get(col_anno[3])), top_annotation_height = unit(6, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==4){
  Heatmap(selected_value, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2]), col_anno_3=get(col_anno[3]), col_anno_4=get(col_anno[4])), top_annotation_height = unit(8, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
}
dev.off()
write.table(paste(outPath,fileNamePrefix,"complex.SVG",sep=""),argv[25],row.names=FALSE,col.names=FALSE)
