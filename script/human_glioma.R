library(miRTalk)
d1 <- Read10X_h5("Glioma_GSE84465_expression.h5")
d2 <- read.table('Glioma_GSE84465_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
a2 <- d2[,c("Cell","Celltype..minor.lineage.")]
colnames(a2)[2] <- "Celltype"
a2[a2$Celltype == "AC-like Malignant",]$Celltype <- "Malignant"
a2[a2$Celltype == "M1",]$Celltype <- "Monocyte"
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
obj <- create_miRTalk(sc_data = d1, sc_celltype = a2$Celltype,
                      species = "Human", condition = rep("cancer", ncol(d1)), evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, use_n_cores = 16)


######################### miRNA mapping
mrna_exp <- readRDS("~/workspace/miRTalk/test_data/glioblastoma/GSE84465_OK/revise/tcga_gbm_rna.rds")
mrna_meta <- as.data.frame(mrna_exp@colData)
mrna_data_list <- mrna_exp@assays@data@listData
mrna_data1 <- mrna_data_list[[1]]
mrna_gene <- as.data.frame(mrna_exp@rowRanges)
mrna_data_list <- lapply(mrna_data_list, function(x){
  rownames(x) <- rownames(mrna_gene)
  genenames <- as.data.frame(table(mrna_gene$gene_name),stringsAsFactors = F)
  genenames1 <- genenames[genenames$Freq == 1,]
  genenames2 <- genenames[genenames$Freq > 1,]
  mrna_gene1 <- mrna_gene[mrna_gene$gene_name %in% genenames1$Var1,]
  mrna_gene2 <- mrna_gene[mrna_gene$gene_name %in% genenames2$Var1,]
  x1 <- x[rownames(mrna_gene1),]
  rownames(x1) <- mrna_gene1$gene_name
  x2 <- x[rownames(mrna_gene2),]
  x3 <- x2[grep(pattern = "PAR_Y",rownames(x2)),]
  x2 <- x2[!rownames(x2) %in% rownames(x3),]
  x2_list <- list()
  for (i in 1:nrow(genenames2)) {
    a1 <- mrna_gene2[mrna_gene2$gene_name == genenames2$Var1[i],]
    a1 <- a1[rownames(a1) %in% rownames(x2),]
    x_tmp <- x2[rownames(a1),]
    if (nrow(a1) > 1) {
      x_tmp <- colMeans(x_tmp)
    }
    x2_list[[i]] <- x_tmp
  }
  x2_list <- as.data.frame(t(as.data.frame(x2_list)))
  rownames(x2_list) <- genenames2$Var1
  x_final <- rbind(x1,x2_list)
  colnames(x_final) <- mrna_meta$sample
  x_final <- rev_gene(x_final, data_type = "count",species = "Human",geneinfo = geneinfo)
  return(x_final)
})
genes_mirna <- mrna_gene[mrna_gene$gene_type == "miRNA",]
res_all_sig <- readRDS("~/workspace/miRTalk/test_data/glioblastoma/GSE84465_OK/revise/res_survival_tcga_sig.rds")
mir_info_human <- mir_info[mir_info$species == "Human",c(1,2,3)]
mir_genes <- mir_info_human[mir_info_human$miRNA_mature %in% res_all_sig$miRNA,]
mir_genes <- unique(mir_genes$gene)
mrna_data_scoring <- lapply(mrna_data_list, function(x){
  x_mirrna <- x[rownames(x) %in% genes_mirna$gene_name,]
  x_mirrna <- CreateSeuratObject(x_mirrna)
  x_mirrna <- NormalizeData(x_mirrna)
  x_mirrna <- x_mirrna[["RNA"]]@data
  x_mirrna <- x_mirrna[rownames(x_mirrna) %in% mir_genes,]
  # x_condition
  x1_condition <- mrna_meta[,c("sample","definition")]
  x1_condition$colnames <- colnames(x_mirrna)
  colnames(x_mirrna) <- x1_condition$definition
  x1_test <- apply(x_mirrna, 1, function(x){
    x_normal <- x[names(x) == "Solid Tissue Normal"]
    x_tumor <- x[names(x) != "Solid Tissue Normal"]
    w_test <- as.numeric(wilcox.test(x_normal, x_tumor)$p.value)
    t_test1 <- as.numeric(t.test(x_normal, x_tumor)$p.value)
    t_test2 <- as.numeric(t.test(x_normal, x_tumor, var.equal = T)$p.value)
    return(c(mean(x_normal), mean(x_tumor),w_test,t_test1,t_test2))
  })
  x1_test <- t(x1_test)
  x1_test[is.nan(x1_test)] <- 1
  colnames(x1_test) <- c("mean_normal", "mean_tumor", "w_test", "t_test1", "t_test2")
  x1_test <- as.data.frame(x1_test)
})
res_all_sig_tmp <- res_all_sig[,c("miRNA","targetgene","mti_type")]
colnames(res_all_sig_tmp)[1] <- "miRNA_mature"
mrna_data_merge <- list()
for (i in 1:length(mrna_data_scoring)) {
  x1_test <- mrna_data_scoring[[i]]
  x1_test$change <- apply(x1_test[,c(1,2)], 1, function(x){
    if (x[1] > x[2]) {
      return("decrease")
    }
    if (x[1] < x[2]) {
      return("increase")
    }
    if (x[1] == x[2]) {
      return("no_change")
    }
  })
  x1_test$sig_w <- "not"
  x1_test[x1_test$w_test < 0.1,]$sig_w <- "sig"
  x1_test$sig_t1 <- "not"
  x1_test[x1_test$t_test1 < 0.1,]$sig_t1 <- "sig"
  x1_test$sig_t2 <- "not"
  x1_test[x1_test$t_test2 < 0.1,]$sig_t2 <- "sig"
  d_tmp <- data.frame()
  for (j in 1:nrow(x1_test)) {
    mirgenename <- rownames(x1_test)[j]
    d1 <- mir_info_human[mir_info_human$gene == mirgenename,]
    d2 <- x1_test[d1$gene,]
    d2$miRNA_mature <- d1$miRNA_mature
    d2$miRNA_gene <- d1$gene
    d_tmp <- rbind(d_tmp, d2)
  }
  d_tmp <- d_tmp[d_tmp$miRNA_mature %in% res_all_sig_tmp$miRNA,]
  d_tmp <- merge(res_all_sig_tmp, d_tmp)
  mrna_data_merge[[i]] <- d_tmp
}

d_tmp <- mrna_data_merge[[2]]


######################### onco-MiTI and TS-MiTI plot
res_summary <- readRDS("res_summary.rds")
res_summary <- res_summary[res_summary$summary != "unsure",]
evbiog_human <- evbiog[evbiog$species == "Human",]$gene
risc_human <- risc[risc$species == "Human",]$gene
tcga_patient <- read.table('data_clinical_patient.txt',sep = "\t",header = T)
tcga_patient$time_os <- tcga_patient$OS_MONTHS
tcga_patient$status_os <- 0
tcga_patient[tcga_patient$OS_STATUS == "1:DECEASED",]$status_os <- 1
tcga_patient$time_dfs <- tcga_patient$DFS_MONTHS
tcga_patient$status_dfs <- 0
tcga_patient[tcga_patient$DFS_STATUS == "1:Recurred",]$status_dfs <- 1

# miRNA
tcga_mirna1 <- readr::read_table('data_mirna.txt')
# mRNA
tcga_mrna2 <- readr::read_table('data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt')
score_minmaxfun <- function(x){
  x_minmax <- (x-min(x))/(max(x)-min(x))
  return(x_minmax)
}
score_rankfun <- function(x){
  x_rank <- rank(x)
  x_rank <- x_rank/max(x_rank)
  return(x_rank)
}
genenames <- as.data.frame(table(tcga_mrna2$Hugo_Symbol))
genenames <- genenames[genenames$Freq == 1,]
tcga_mrna2 <- tcga_mrna2[tcga_mrna2$Hugo_Symbol %in% as.character(genenames$Var1),]
tcga_mrna2 <- as.data.frame(tcga_mrna2)
rownames(tcga_mrna2) <- tcga_mrna2$Hugo_Symbol
tcga_mrna2 <- tcga_mrna2[,-1]
colnames(tcga_mrna2) <- c(colnames(tcga_mrna2)[-1],"other")
tcga_mrna2 <- tcga_mrna2[,1:206]
tcga_mrna2 <- as.matrix(tcga_mrna2)
tcga_mrna2[which(is.na(tcga_mrna2))] <- 0
tcga_mrna2 <- rev_gene(data = tcga_mrna2, data_type = "count",species = "Human", geneinfo = geneinfo)
# scoring
obj_seurat1 <- CreateSeuratObject(tcga_mrna2)
obj_seurat1 <- AddModuleScore(obj_seurat1, features = list(c1 = evbiog_human, c2 = risc_human))
moduelscore <- data.frame(patient = colnames(tcga_mrna2), evbiog2 = obj_seurat1$Cluster1, risc2 = obj_seurat1$Cluster2,stringsAsFactors = F)
moduelscore$evbiog2 <- score_minmaxfun(moduelscore$evbiog2)
moduelscore$risc2 <- score_minmaxfun(moduelscore$risc2)
# miRNA
tcga_mirna1 <- as.data.frame(tcga_mirna1)
rownames(tcga_mirna1) <- tcga_mirna1$Hugo_Symbol
tcga_mirna1 <- tcga_mirna1[,-1]
tcga_mirna1 <- as.matrix(tcga_mirna1)
tcga_mirna1[which(is.na(tcga_mirna1))] <- 0
miR_samples <- intersect(moduelscore$patient, colnames(tcga_mirna1))
evbiog_score2 <- moduelscore[miR_samples,]$evbiog2
risc_score2 <- moduelscore[miR_samples,]$risc2
rownames(tcga_patient) <- paste0(tcga_patient$PATIENT_ID,"-01")
tcga_patient <- tcga_patient[miR_samples,]
tcga_targene2 <- apply(tcga_mrna2, 2, function(x){
  x_rank <- rank(-x)
  x_rank <- x_rank/max(x_rank)
  return(x_rank)
})

library(survminer)
for (i in 1:nrow(res_summary)) {
  mirname_new <- res_summary$miRNA_matched_cbiop[i]
  targetgene <- res_summary$targetgene[i]
  type <- res_summary$mti_type[i]
  mirname_exp1 <- tcga_mirna1[mirname_new, miR_samples]
  targetgene_exp2 <- tcga_targene2[targetgene, miR_samples]
  sscore_mir1_gene2 <- evbiog_score2*mirname_exp1*targetgene_exp2*risc_score2
  tcga_patient1 <- tcga_patient
  tcga_patient1$score <- sscore_mir1_gene2
  tcga_patient1 <- tcga_patient1[order(-tcga_patient1$score),]
  tcga_patient1$group <- "low"
  tcga_patient1[1:99,]$group <- "high"
  # os
  if (res_summary$os_p[i] < res_summary$dfs_p[i]) {
    fit <- survfit(Surv(time_os,status_os) ~ group,  data = tcga_patient1)
    p1 <- ggsurvplot(fit,
                     data = tcga_patient1,
                     conf.int = TRUE,
                     pval = TRUE,
                     risk.table = F,
                     surv.median.line = "hv",
                     add.all = F,
                     palette = "hue")
    pdf(file = paste0(i,"_",type,"_",mirname_new, "-",targetgene,"_os.pdf"))
    print(p1)
    dev.off()
  }
  if (res_summary$os_p[i] > res_summary$dfs_p[i]) {
    fit <- survfit(Surv(time_dfs,status_dfs) ~ group,  data = tcga_patient1)
    p1 <- ggsurvplot(fit,
                     data = tcga_patient1,
                     conf.int = TRUE,
                     pval = TRUE,
                     risk.table = F,
                     surv.median.line = "hv",
                     add.all = F,
                     palette = "hue")
    pdf(file = paste0(i,"_",type,"_",mirname_new, "-",targetgene,"_dfs.pdf"))
    print(p1)
    dev.off()
  }
}

for (i in 1:nrow(res_summary)) {
  mirname_new <- res_summary$miRNA_matched_cbiop[i]
  targetgene <- res_summary$targetgene[i]
  type <- res_summary$mti_type[i]
  mirname_exp1 <- tcga_mirna1[mirname_new, miR_samples]
  mirname_exp1 <- (mirname_exp1-min(mirname_exp1))/(max(mirname_exp1)-min(mirname_exp1))
  targetgene_exp2 <- tcga_mrna2[targetgene, miR_samples]
  targetgene_exp2 <- (targetgene_exp2-min(targetgene_exp2))/(max(targetgene_exp2)-min(targetgene_exp2))
  testdata <- data.frame(mir = as.numeric(mirname_exp1), target = as.numeric(targetgene_exp2),stringsAsFactors = F)
  a1 <- correlation::correlation(data = testdata,method = "spearman")
  spear_r <- round(as.numeric(a1$r),4)
  spear_p <- as.numeric(a1$p)
  text_data <- data.frame(x = 0.5, y = 1, text = paste0(spear_r, ";", spear_p),stringsAsFactors = F)
  p1 <- ggplot()+geom_point(data = testdata, aes(x = mir, y = target))+geom_smooth(data = testdata,aes(x = mir, y = target),method = "lm")+geom_text(mapping = aes(x = x,y=y, label = text),data = text_data)
  ggsave(filename = paste0(i,"_",type,"_",mirname_new, "-",targetgene,"_cor.pdf"),plot = p1,device = "pdf",width = 6,height = 6)
}
# change
res_all_sig <- readRDS("res_survival_tcga_sig.rds")
colnames(res_all_sig)[1] <- "miRNA_mature"
mir_info_human <- mir_info[mir_info$species == "Human", c(1,2)]
mir_info_human <- merge(res_all_sig, mir_info_human)
cci_miRnames <- unique(mir_info_human$miRNA)
mir_exp <- readRDS("tcga_gbm_mirna.rds")
mir_exp_count <- mir_exp[,c(1, grep(pattern = "read_count", x = colnames(mir_exp)))]
rownames(mir_exp_count) <- mir_exp_count$miRNA_ID
mir_exp_count <- mir_exp_count[,-1]
mir_exp_count <- CreateSeuratObject(mir_exp_count)
mir_exp_count <- NormalizeData(mir_exp_count)
mir_exp_count <- mir_exp_count[["RNA"]]@data
mir_exp_meta <- data.frame(sample = colnames(mir_exp_count), condition = c("normal", "tumor",
                                                                           "normal","tumor","tumor","tumor","tumor","normal","normal",
                                                                           "tumor", "tumor", "tumor", "normal"), stringsAsFactors = F)
mirnames <- rownames(mir_exp_count)
mirnames <- stringr::str_replace(string = mirnames,pattern = "hsa-mir",replacement = "hsa-miR")
mirnames <- data.frame(mir = mirnames, mir2 = mirnames,stringsAsFactors = F)
#
cci_miRnames_map <- cci_miRnames[cci_miRnames %in% mirnames$mir]
cci_miRnames_unmap <- cci_miRnames[!cci_miRnames %in% mirnames$mir]
rownames(mir_exp_count) <- mirnames$mir
mir_exp_list <- list()
for (i in 1:length(cci_miRnames_unmap)) {
  d1_names <- mirnames$mir[grep(pattern = cci_miRnames_unmap[i], x = mirnames$mir)]
  mir_exp_list[[i]] <- colMeans(mir_exp_count[d1_names,])
}
mir_exp_list <- as.data.frame(t(as.data.frame(mir_exp_list)))
mir_exp_list <- as(object = as.matrix(mir_exp_list), Class = 'dgCMatrix')
rownames(mir_exp_list) <- cci_miRnames_unmap
colnames(mir_exp_list) <- mir_exp_meta$condition
mir_exp_count <- mir_exp_count[cci_miRnames_map,]
colnames(mir_exp_count) <- mir_exp_meta$condition
mir_exp_count <- rbind(mir_exp_count,mir_exp_list)
for (i in 1:nrow(res_summary)) {
  mirname_new <- res_summary$miRNA_matched[i]
  type <- res_summary$mti_type[i]
  x_normal <- mir_exp_count[mirname_new,]
  x_normal <- as.numeric(x_normal[names(x_normal) == "normal"])
  x_tumor <- mir_exp_count[mirname_new,]
  x_tumor <- as.numeric(x_tumor[names(x_tumor) == "tumor"])
  data_plot <- data.frame(x = c("normal", "tumor"), y = c(mean(x_normal), mean(x_tumor)),
                          ymin = c(mean(x_normal) - sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) - sd(x_tumor)/sqrt(length(x_tumor))),
                          ymax = c(mean(x_normal) + sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) + sd(x_tumor)/sqrt(length(x_tumor))),stringsAsFactors = F)

  p1 <- ggplot(data_plot)+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))
  ggsave(filename = paste0(i,"_",type,"_",mirname_new, "_change.png"),plot = p1,device = "png",width = 6,height = 6)
}


## cor
mirname_new <- "hsa-miR-10b"
targetgene <- "PABPC1"
mirname_exp1 <- tcga_mirna1[mirname_new, miR_samples]
mirname_exp1 <- (mirname_exp1-min(mirname_exp1))/(max(mirname_exp1)-min(mirname_exp1))
targetgene_exp2 <- tcga_mrna2[targetgene, miR_samples]
targetgene_exp2 <- (targetgene_exp2-min(targetgene_exp2))/(max(targetgene_exp2)-min(targetgene_exp2))
testdata <- data.frame(mir = as.numeric(mirname_exp1), target = as.numeric(targetgene_exp2),stringsAsFactors = F)
a1 <- correlation::correlation(data = testdata,method = "spearman")
spear_r <- round(as.numeric(a1$r),4)
spear_p <- as.numeric(a1$p)
text_data <- data.frame(x = 0.5, y = 1, text = paste0(spear_r, ";", spear_p),stringsAsFactors = F)
p1 <- ggplot()+geom_point(data = testdata, aes(x = mir, y = target))+geom_smooth(data = testdata,aes(x = mir, y = target),method = "lm")+geom_text(mapping = aes(x = x,y=y, label = text),data = text_data)
ggsave(filename = paste0("Other",mirname_new, "-",targetgene,"_cor.pdf"),plot = p1,device = "pdf",width = 6,height = 6)

mirname_new <- "hsa-miR-211"
targetgene <- "PABPC1"
mirname_exp1 <- tcga_mirna1[mirname_new, miR_samples]
mirname_exp1 <- (mirname_exp1-min(mirname_exp1))/(max(mirname_exp1)-min(mirname_exp1))
targetgene_exp2 <- tcga_mrna2[targetgene, miR_samples]
targetgene_exp2 <- (targetgene_exp2-min(targetgene_exp2))/(max(targetgene_exp2)-min(targetgene_exp2))
testdata <- data.frame(mir = as.numeric(mirname_exp1), target = as.numeric(targetgene_exp2),stringsAsFactors = F)
a1 <- correlation::correlation(data = testdata,method = "spearman")
spear_r <- round(as.numeric(a1$r),4)
spear_p <- as.numeric(a1$p)
text_data <- data.frame(x = 0.5, y = 1, text = paste0(spear_r, ";", spear_p),stringsAsFactors = F)
p1 <- ggplot()+geom_point(data = testdata, aes(x = mir, y = target))+geom_smooth(data = testdata,aes(x = mir, y = target),method = "lm")+geom_text(mapping = aes(x = x,y=y, label = text),data = text_data)
ggsave(filename = paste0("Other",mirname_new, "-",targetgene,"_cor.pdf"),plot = p1,device = "pdf",width = 6,height = 6)
