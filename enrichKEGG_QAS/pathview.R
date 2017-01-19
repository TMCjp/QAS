#安装需要的R包
pkgs <- rownames(installed.packages())
#本程序所需包
packages_1=c('pathview')
#将要安装的R包
packages_2=setdiff(packages_1,pkgs) 
#安装所需R包 
if(length(packages_2)>0){
  install.packages(packages_2, repos="http://cran.rstudio.com/")	
  source("http://bioconductor.org/biocLite.R")
  biocLite(packages_2)
}

args <- commandArgs(TRUE)

diff_gene <- read.table(args[1], header=TRUE, sep=",", na.strings=c('NA'))
enrichKEGG <- read.table(args[2], header=TRUE, sep=",", na.strings=c('NA'))

library("pathview")

KEGG_random <- runif(1,0,10000000000)
KEGG_random <- round(KEGG_random)

galaxy_dir <- getwd()

setwd("/var/www/QAS_enrichKEGG/result/")
hsa04110 <- pathview(gene.data  = as.character(diff_gene$symbol),
                     pathway.id = as.character(enrichKEGG$ID[[1]]),
                     species    = args[3],   #'hsa'
                     out.suffix = KEGG_random, 
                     gene.idtype = "SYMBOL",
                     kegg.dir = "/var/www/QAS_enrichKEGG/pathway/")

result_dir <- "/var/www/QAS_enrichKEGG/result/"

setwd(galaxy_dir)
results <- paste0(result_dir, as.character(enrichKEGG$ID[[1]]), '.', KEGG_random, '.png')

write.table(results, args[4], sep = ',', row.names = FALSE, col.names = FALSE)