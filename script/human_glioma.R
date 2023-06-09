library(miRTalk)
d1 <- Read10X_h5("Glioma_GSE84465_expression.h5")
d2 <- read.table('Glioma_GSE84465_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)

a2 <- d2[,c("Cell","Celltype..minor.lineage.")]
colnames(a2)[2] <- "Celltype"
a2[a2$Celltype == "AC-like Malignant",]$Celltype <- "Malignant"
a2[a2$Celltype == "M1",]$Celltype <- "Monocyte"

obj <- create_miRTalk(sc_data = d1, sc_celltype = a2$Celltype, species = "Human",if_normalize = F)
obj <- find_miRNA(obj, mir_info = mir_info)
obj <- find_hvtg(obj)
obj <- find_miRTalk(obj, mir2tar = mir2tar,use_n_cores = 16)

celltype_col <- ggsci::pal_d3(palette = "category10")(7)
plot_miRTalk_sankey(object = obj,celltype_color = celltype_col,node_pad = 0)
plot_miRTalk_chord(obj,celltype_color = celltype_col)
plot_miR_heatmap(obj)

############### bubble
cci_sim <- get_miRTalk_cci(obj)
cci1 <- cci_sim[cci_sim$celltype_sender == "Malignant",]
cci1 <- cci1[cci1$miRNA == "hsa-miR-125b-5p",]

y_name <-  unique(c(cci1$celltype_sender, cci1$celltype_receiver))
x_name <- unique(cci1$target_gene)

y_name <- data.frame(names = y_name, y = 7:1,stringsAsFactors = F)
x_name <- data.frame(names = x_name, x = 1:length(x_name),stringsAsFactors = F)

cci1$x <- 0
cci1$y <- 0

for (i in 1:nrow(cci1)) {
  d1 <- cci1$celltype_receiver[i]
  d1 <- y_name[y_name$names == d1,]$y
  cci1$y[i] <- d1

  d1 <- cci1$target_gene[i]
  d1 <- x_name[x_name$names == d1,]$x
  cci1$x[i] <- d1
}

cci12 <- cci1[1,]
cci12$x <- 1
cci12$y <- 1
cci12$score <- 0
cci12$prob <- min(cci1$prob)
cci1 <- rbind(cci1,cci12)
p1 <- ggplot2::ggplot(data = cci1) + ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = prob,
    size = score)) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+scale_color_viridis(option = "C")+
    ggplot2::scale_x_continuous(name = "targets", breaks = x_name$x, labels = x_name$names, limits = c(1, max(x_name$x)))+
    scale_size_continuous(range = c(2,10))
