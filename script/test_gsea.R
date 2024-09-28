# function

library(Seurat)
library(clusterProfiler)
devtools::load_all('~/github/miRTalk/')
# bladder cancer
setwd('~/workspace/miRTalk/test_data/bladder/')
d1 <- Read10X_h5("BLCA_GSE130001_expression.h5")
d2 <- read.table('BLCA_GSE130001_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$Sample)

evbiog <- evbiog[evbiog$species == "Human",]
colnames(evbiog) <- c("term", "gene")

res_tmp <- NULL
for (i in 1:length(sampleid)) {
  sc_meta <- d2[d2$Sample == sampleid[i],]
  sc_data <- d1[,sc_meta$Cell]
  celltypes <- unique(sc_meta$Celltype..major.lineage.)
  res_tmp_sample <- data.frame(sample = sampleid[i], celltypes = celltypes, gsea_p = 1, stringsAsFactors = F)
  for (j in 1:length(celltypes)) {
      sc_meta1 <- sc_meta[sc_meta$Celltype..major.lineage. == celltypes[j],]
      sc_data1 <- sc_data[,sc_meta1$Cell]
      sc_data1 <- rowMeans(sc_data1)
      sc_data1 <- scale(sc_data1)
      genelist <- as.numeric(sc_data1)
      names(genelist) <- rownames(sc_data1)
      genelist <- genelist[order(-genelist)]
      gsea_res <- GSEA(geneList = genelist, TERM2GENE = evbiog,verbose = F)
      gsea_res <- gsea_res@result
      if (nrow(gsea_res) == 1) {
          res_tmp_sample$gsea_p[j] <- gsea_res$pvalue
      }
  }
  res_tmp <- rbind(res_tmp, res_tmp_sample)
}

# chol
setwd('~/workspace/miRTalk/test_data/chol/')
load("TLMA-SC-A044/sc_data_raw.rda")
load("TLMA-SC-A044/sc_meta.rda")
d1 <- sc_data_raw
d2 <- sc_meta
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
d2[d2$Specimen %in% c("ICC24_1", "ICC24_2"),]$Specimen <- "ICC24"
sampleid <- unique(d2$Specimen)
d1 <- CreateSeuratObject(d1)
d1 <- NormalizeData(d1)
d1 <- d1[["RNA"]]@data

res_tmp <- NULL
for (i in 1:length(sampleid)) {
  sc_meta <- d2[d2$Specimen == sampleid[i],]
  sc_data <- d1[,sc_meta$Id]
  celltypes <- unique(sc_meta$Name_Celltype)
  res_tmp_sample <- data.frame(sample = sampleid[i], celltypes = celltypes, gsea_p = 1, stringsAsFactors = F)
  for (j in 1:length(celltypes)) {
    sc_meta1 <- sc_meta[sc_meta$Name_Celltype == celltypes[j],]
    sc_data1 <- sc_data[,sc_meta1$Id]
    if (nrow(sc_meta1) > 1) {
        sc_data1 <- rowMeans(sc_data1)
        sc_data1 <- scale(sc_data1)
        genelist <- as.numeric(sc_data1)
        names(genelist) <- rownames(sc_data1)
    }
    genelist <- genelist[order(-genelist)]
    gsea_res <- GSEA(geneList = genelist, TERM2GENE = evbiog,verbose = F)
    gsea_res <- gsea_res@result
    if (nrow(gsea_res) == 1) {
      res_tmp_sample$gsea_p[j] <- gsea_res$pvalue
    }
  }
  res_tmp <- rbind(res_tmp, res_tmp_sample)
}

# ovarian
setwd('~/workspace/miRTalk/test_data/OV/')
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "treatment-naive",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$patient_id)
d1 <- CreateSeuratObject(d1)
d1 <- NormalizeData(d1)
d1 <- d1[["RNA"]]@data

res_tmp <- NULL
for (i in 1:length(sampleid)) {
  sc_meta <- d2[d2$patient_id == sampleid[i],]
  sc_data <- d1[,sc_meta$cell]
  celltypes <- unique(sc_meta$cell_subtype)
  res_tmp_sample <- data.frame(sample = sampleid[i], celltypes = celltypes, gsea_p = 1, stringsAsFactors = F)
  for (j in 1:length(celltypes)) {
    sc_meta1 <- sc_meta[sc_meta$cell_subtype == celltypes[j],]
    sc_data1 <- sc_data[,sc_meta1$cell]
    if (nrow(sc_meta1) > 1) {
      sc_data1 <- rowMeans(sc_data1)
      sc_data1 <- scale(sc_data1)
      genelist <- as.numeric(sc_data1)
      names(genelist) <- rownames(sc_data1)
    }
    genelist <- genelist[order(-genelist)]
    gsea_res <- GSEA(geneList = genelist, TERM2GENE = evbiog,verbose = F)
    gsea_res <- gsea_res@result
    if (nrow(gsea_res) == 1) {
      res_tmp_sample$gsea_p[j] <- gsea_res$pvalue
    }
  }
  res_tmp <- rbind(res_tmp, res_tmp_sample)
}


d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "post-NACT",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$patient_id)
d1 <- CreateSeuratObject(d1)
d1 <- NormalizeData(d1)
d1 <- d1[["RNA"]]@data

res_tmp <- NULL
for (i in 1:length(sampleid)) {
  sc_meta <- d2[d2$patient_id == sampleid[i],]
  sc_data <- d1[,sc_meta$cell]
  celltypes <- unique(sc_meta$cell_subtype)
  res_tmp_sample <- data.frame(sample = sampleid[i], celltypes = celltypes, gsea_p = 1, stringsAsFactors = F)
  for (j in 1:length(celltypes)) {
    sc_meta1 <- sc_meta[sc_meta$cell_subtype == celltypes[j],]
    sc_data1 <- sc_data[,sc_meta1$cell]
    if (nrow(sc_meta1) > 1) {
      sc_data1 <- rowMeans(sc_data1)
      sc_data1 <- scale(sc_data1)
      genelist <- as.numeric(sc_data1)
      names(genelist) <- rownames(sc_data1)
    }
    genelist <- genelist[order(-genelist)]
    gsea_res <- GSEA(geneList = genelist, TERM2GENE = evbiog,verbose = F)
    gsea_res <- gsea_res@result
    if (nrow(gsea_res) == 1) {
      res_tmp_sample$gsea_p[j] <- gsea_res$pvalue
    }
  }
  res_tmp <- rbind(res_tmp, res_tmp_sample)
}


# kidney
samplename <- unique(sc_meta$Time.point)
sc_data <- rev_gene(data = sc_data, data_type = "count",species = "Mouse",geneinfo = geneinfo)

d2 <- sc_meta
d1 <- sc_data

evbiog <- evbiog[evbiog$species == "Mouse",]

res_tmp <- NULL
sampleid <- samplename[1]
for (i in 1:length(sampleid)) {
  sc_meta <- d2[d2$Time.point == sampleid[i],]
  sc_data <- d1[,sc_meta$id]
  celltypes <- unique(sc_meta$Annotation.Level.2)
  res_tmp_sample <- data.frame(sample = sampleid[i], celltypes = celltypes, gsea_p = 1, stringsAsFactors = F)
  for (j in 1:length(celltypes)) {
    sc_meta1 <- sc_meta[sc_meta$Annotation.Level.2 == celltypes[j],]
    sc_data1 <- sc_data[,sc_meta1$id]
    if (nrow(sc_meta1) > 1) {
      sc_data1 <- rowMeans(sc_data1)
      sc_data1 <- scale(sc_data1)
      genelist <- as.numeric(sc_data1)
      names(genelist) <- rownames(sc_data1)
    }
    genelist <- genelist[order(-genelist)]
    gsea_res <- GSEA(geneList = genelist, TERM2GENE = evbiog,verbose = F)
    gsea_res <- gsea_res@result
    if (nrow(gsea_res) == 1) {
      res_tmp_sample$gsea_p[j] <- gsea_res$pvalue
    }
  }
  res_tmp <- rbind(res_tmp, res_tmp_sample)
}

res_tmp <- NULL
sampleid <- samplename[1]
for (i in 1:length(sampleid)) {
  sc_meta <- d2[d2$Time.point != sampleid[i],]
  sc_data <- d1[,sc_meta$id]
  celltypes <- unique(sc_meta$Annotation.Level.2)
  res_tmp_sample <- data.frame(sample = sampleid[i], celltypes = celltypes, gsea_p = 1, stringsAsFactors = F)
  for (j in 1:length(celltypes)) {
    sc_meta1 <- sc_meta[sc_meta$Annotation.Level.2 == celltypes[j],]
    sc_data1 <- sc_data[,sc_meta1$id]
    if (nrow(sc_meta1) > 1) {
      sc_data1 <- rowMeans(sc_data1)
      sc_data1 <- scale(sc_data1)
      genelist <- as.numeric(sc_data1)
      names(genelist) <- rownames(sc_data1)
    }
    genelist <- genelist[order(-genelist)]
    gsea_res <- GSEA(geneList = genelist, TERM2GENE = evbiog,verbose = F)
    gsea_res <- gsea_res@result
    if (nrow(gsea_res) == 1) {
      res_tmp_sample$gsea_p[j] <- gsea_res$pvalue
    }
  }
  res_tmp <- rbind(res_tmp, res_tmp_sample)
}  
  
  
  
# TLMA-SC-A031
setwd('~/workspace/miRTalk/test_data/liver/TLMA-SC-A031')
load('sc_data_raw.rda')
load('sc_meta.rda')
load('sc_meta_specimen.rda')
sc_data_raw <- rev_gene(sc_data_raw,data_type = 'count',species = 'Rat',geneinfo = geneinfo)

evbiog <- evbiog[evbiog$species == "Rat",]

d1 <- sc_data_raw
d2 <- sc_meta
d2$Status <- "Normal"
d2[d2$Specimen %in% c("FDL1", "FDL2", "FDL3"),]$Status <- "NAFLD (high-fat diet, 8 weeks)"
samplename <- unique(sc_meta_specimen$Status)

res_tmp <- NULL

for (i in 1:length(samplename)) {
  sc_meta <- d2[d2$Status != samplename[i],]
  sc_data <- d1[,sc_meta$Id]
  celltypes <- unique(sc_meta$Name_Celltype)
  res_tmp_sample <- data.frame(sample = samplename[i], celltypes = celltypes, gsea_p = 1, stringsAsFactors = F)
  for (j in 1:length(celltypes)) {
    sc_meta1 <- sc_meta[sc_meta$Name_Celltype == celltypes[j],]
    sc_data1 <- sc_data[,sc_meta1$Id]
    if (nrow(sc_meta1) > 1) {
      sc_data1 <- rowMeans(sc_data1)
      sc_data1 <- scale(sc_data1)
      genelist <- as.numeric(sc_data1)
      names(genelist) <- rownames(sc_data1)
    }
    genelist <- genelist[order(-genelist)]
    gsea_res <- GSEA(geneList = genelist, TERM2GENE = evbiog,verbose = F)
    gsea_res <- gsea_res@result
    if (nrow(gsea_res) == 1) {
      res_tmp_sample$gsea_p[j] <- gsea_res$pvalue
    }
  }
  res_tmp <- rbind(res_tmp, res_tmp_sample)
}  