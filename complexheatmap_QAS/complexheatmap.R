suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(colorspace)))
suppressWarnings(suppressMessages(library(GetoptLong)))
suppressWarnings(suppressMessages(library(RSvgDevice)))

argv <- commandArgs(TRUE);
inputPath = argv[1]
cluster_rows = as.logical(argv[3])
row_hclust_side = argv[4]
row_names_side = argv[5]
cluster_columns = as.logical(argv[6])
column_names_side = argv[7]
column_hclust_side = argv[8]
km = as.numeric(argv[9])
outPath = "/var/www/heatmapoutput/"
numberRange = c(as.numeric(argv[10]),as.numeric(argv[11]),as.numeric(argv[12]))
colorRange = c(argv[13],argv[14],argv[15])
colorRange = gsub(pattern="__pd__", replacement="#", colorRange)

ha_mix_col_his_val = argv[16]
ha_mix_col_den_val = argv[17]
ha_mix_col_vio_val = argv[18]
ha_mix_col_heatmap_val = argv[19]
ha_mix_row_his_val = argv[20]
ha_mix_row_den_val = argv[21]
ha_mix_row_vio_val = argv[22]
ha_mix_row_heatmap_val = argv[23]

fileNamePrefix <- runif(1,0,1000000000)
fileNamePrefix <- round(fileNamePrefix)
mat <- read.table(inputPath, header=TRUE,sep=",")
mat <- as.matrix(mat)

mat <- na.omit(mat)
mat_unique <- unique(mat[,1])
math <- match(mat_unique,mat[,1])
mat <- mat[math,]
rownames(mat) <- as.factor(mat[,1])
mat <- mat[,-1]
#只取前50行的数据，进行展示
mat <- mat[1:50,]

protein <- mat
protein[is.na(protein)] <- 0

ff_1 <- hclust(dist(t(protein)))
protein <- protein[, ff_1$order]

kc <- kmeans(protein, km)
cla <- kc$cluster
clas <- unique(kc$cluster)
data <- NULL
for(i in 1:km) {
	cs <- names(cla[cla==i])
  #cs <- names(cla[cla==clas[i]])
  protein_us <- protein[cs, ]
  ff <- hclust(dist(protein_us))
  write.csv(as.data.frame(protein_us[ff$order, ]), paste(outPath,fileNamePrefix,"inte_", i,".csv", sep="", collapse = ""))
}


ha_mix_col_his = anno_histogram(mat, which = "column")                                                             

ha_mix_col_den = anno_density(mat, type = "line", which = "column")

ha_mix_col_vio = anno_density(mat, type = "violin", which = "column")

ha_mix_col_heatmap = anno_density(mat, type = "heatmap", which = "column")

ha_mix_row_his = HeatmapAnnotation(histogram = anno_histogram(mat, which = "row"),which = "row",width = unit(2, "cm"))
ha_mix_row_den = HeatmapAnnotation(density_line = anno_density(mat, type = "line", which = "row"),which = "row",width = unit(2, "cm"))
ha_mix_row_vio = HeatmapAnnotation(violin = anno_density(mat, type = "violin", which = "row"),which = "row",width = unit(2, "cm"))
ha_mix_row_heatmap = HeatmapAnnotation(heatmap = anno_density(mat, type = "heatmap", which = "row"), which = "row",width = unit(2, "cm"))


col_anno <- c(ha_mix_col_his_val, ha_mix_col_den_val, ha_mix_col_vio_val, ha_mix_col_heatmap_val)
row_anno <- c(ha_mix_row_his_val, ha_mix_row_den_val, ha_mix_row_vio_val, ha_mix_row_heatmap_val)


col_anno <- col_anno[which(col_anno!="NULL")]
options(expressions=50000)
devSVG(paste(outPath,fileNamePrefix,"complex.SVG",sep=""))
if(length(col_anno)==0) {
  Heatmap(mat, name = "foo") + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==1) {
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1])), top_annotation_height = unit(2, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==2) {
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2])), top_annotation_height = unit(4, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==3){
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2]), col_anno_3=get(col_anno[3])), top_annotation_height = unit(6, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==4){
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2]), col_anno_3=get(col_anno[3]), col_anno_4=get(col_anno[4])), top_annotation_height = unit(8, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
}
dev.off()
write.table(fileNamePrefix,argv[24],row.names=FALSE,col.names=FALSE)

png(argv[25])
if(length(col_anno)==0) {
  Heatmap(mat, name = "foo") + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==1) {
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1])), top_annotation_height = unit(2, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==2) {
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2])), top_annotation_height = unit(4, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==3){
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2]), col_anno_3=get(col_anno[3])), top_annotation_height = unit(6, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
} else if(length(col_anno)==4){
  Heatmap(mat, name = "foo", cluster_rows =cluster_rows,row_dend_side=row_hclust_side,row_names_side=row_names_side,cluster_columns=cluster_columns,column_names_side=column_names_side,column_dend_side=column_hclust_side, split = data.frame(cla), combined_name_fun = NULL,col = colorRamp2(numberRange,colorRange),top_annotation = HeatmapAnnotation(col_anno_1=get(col_anno[1]), col_anno_2=get(col_anno[2]), col_anno_3=get(col_anno[3]), col_anno_4=get(col_anno[4])), top_annotation_height = unit(8, "cm")) + get(row_anno[1]) + get(row_anno[2]) + get(row_anno[3]) +get(row_anno[4])
}
dev.off()