d1 <- Read10X_h5("glioblastoma/GSE70630_OK/Glioma_GSE70630_expression.h5")
genedata <- data.frame(gene1 = rownames(d1), gene2 = rownames(d1),stringsAsFactors = F)
d2 <- read.table('glioblastoma/GSE70630_OK/Glioma_GSE70630_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
obj <- create_miRTalk(sc_data = d1, sc_celltype = as.character(d2$Celltype..original.), species = "Human")
obj@data$ndata <- d1
obj <- find_cci(obj, mir_info = mir_info, mir2tar = mir2tar)

percent_num <- function(x){
  return(length(x[x > 0]))
}
