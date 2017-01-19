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



#仅针对此数据，统一样本名
colname_mut <- unique(mutation$Tumor_Sample_Barcode)
colname_exp <- unique(colnames(express))
mat <- match(colname_mut,colname_exp)
express <- express[,mat]


#heatmap bar class
library(ComplexHeatmap)
library(circlize)


mut_ha_mix_top = HeatmapAnnotation(violin = anno_density(dat, type = "violin", gp = gpar(fill = "red")))
mut_ha_mix_right = HeatmapAnnotation(histogram = anno_histogram(dat, gp = gpar(fill = "blue"),which = "row"), #row表明是行注释
                                     density_line = anno_density(dat, type = "line",gp = gpar(col = "red"),which = "row"),
                                     boxplot = anno_boxplot(dat,which="row"),
                                which = "row",  #用于行注释，默认是"column"
                                width = unit(3, "cm"))

exp_ha_mix_top = HeatmapAnnotation(violin = anno_density(express, type = "violin", gp = gpar(fill = "red")))
exp_ha_mix_right = HeatmapAnnotation(histogram = anno_histogram(express, gp = gpar(fill = "blue"),which = "row"),
                                which = "row",  #用于行注释，默认是"column"
                                width = unit(1, "cm"))

								
heat_ann <- read.table(args[4])

# ha1 = HeatmapAnnotation(df = colname_mut, col = list(class = c("1" =  "red", "2" = "yellow", "3"= "green")))
# ha2 = HeatmapAnnotation(points = anno_points(rnorm(10)))
# ht1 = Heatmap(dat,cluster_columns = FALSE, name = "ht1", column_title = "Heatmap 1", top_annotation = ha1)
# ht2 = Heatmap(mat, name = "ht2", column_title = "Heatmap 2", top_annotation = ha2, show_heatmap_legend = FALSE)

ha1 = HeatmapAnnotation(df = heat_ann, col = list(class = c("1" =  "red", "2" = "yellow", "3"= "green")))
ht1 = mut_ha_mix_right + Heatmap(dat,cluster_columns = FALSE, cluster_rows = FALSE,
              name = "mutation", column_title = "mutation Heatmap", 
              top_annotation = ha1,col = colorRamp2(c(0, 1, 3), c("RoyalBlue", "green", "red")),
              show_column_names = FALSE,
              bottom_annotation = mut_ha_mix_top)
ht2 = exp_ha_mix_right + Heatmap(express, cluster_columns = FALSE,
              name = "expression", column_title = "expression Heatmap", 
              top_annotation = ha1, col = colorRamp2(c(min(express),mad(express),max(express)), c("RoyalBlue", "green", "red")),
              show_column_names = FALSE,
              bottom_annotation = exp_ha_mix_top)

png(args[5], width = 960, height = 850, units = "px", pointsize = 12)
ht1+ht2

dat <- t(dat)
dat <- cbind(sample=rownames(dat),dat)
dat <- cbind(heat_ann,dat)
write.table(dat,args[6],row.names=FALSE,sep=",")