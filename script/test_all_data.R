devtools::load_all('~/github/miRTalk/')

# bladder cancer
setwd('~/workspace/miRTalk/test_data/bladder/')
d1 <- Read10X_h5("BLCA_GSE130001_expression.h5")
d2 <- read.table('BLCA_GSE130001_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$Celltype..major.lineage.,
                      species = "Human",condition = d2$Sample, evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_bladder.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_bladder.csv")

#---samples
d1 <- Read10X_h5("BLCA_GSE130001_expression.h5")
d2 <- read.table('BLCA_GSE130001_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$Sample)
res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$Sample == sampleid[i],]
  # 1.0
  print(paste0("bladder_",i,"_1.0"))
  a3 <- a2[,c("Cell","Celltype..major.lineage.")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]),
                        species = "Human", condition = rep(sampleid[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = F)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
saveRDS(res_tmp,file = "~/workspace/miRTalk/revise/benchmark/cci_bladder_samples.rds")
write.csv(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_bladder_samples.csv")


# chol
setwd('~/workspace/miRTalk/test_data/chol/')
load("TLMA-SC-A044/sc_data_raw.rda")
load("TLMA-SC-A044/sc_meta.rda")
d1 <- sc_data_raw
d2 <- sc_meta
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
d2[d2$Specimen %in% c("ICC24_1", "ICC24_2"),]$Specimen <- "ICC24"
obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$Celltype,
                      species = "Human",condition = d2$Specimen, evbiog = evbiog, risc = risc, if_normalize = T)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_chol.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_chol.csv")
#---samples
load("TLMA-SC-A044/sc_data_raw.rda")
load("TLMA-SC-A044/sc_meta.rda")
d1 <- sc_data_raw
d2 <- sc_meta
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
d2[d2$Specimen %in% c("ICC24_1", "ICC24_2"),]$Specimen <- "ICC24"
sampleid <- unique(d2$Specimen)
res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$Specimen == sampleid[i],]
  # 1.0
  print(paste0("chol_",i,"_1.0"))
  a3 <- a2[,c("Id","Celltype")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]),
                        species = "Human", condition = rep(sampleid[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = T)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
saveRDS(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_chol_samples.rds")
write.csv(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_chol_samples.csv")


# ovarian
setwd('~/workspace/miRTalk/test_data/OV/')
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
d2$condition <- paste0(d2$treatment_phase,"_", d2$patient_id)
obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$cell_subtype,
                      species = "Human",condition = d2$condition, evbiog = evbiog, risc = risc, if_normalize = T)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_all_filter.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_all_filter.csv")

obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_all.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_all.csv")

#---before
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "treatment-naive",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$cell_subtype,
                      species = "Human",condition = d2$patient_id, evbiog = evbiog, risc = risc, if_normalize = T)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_before_filter.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_before_filter.csv")

obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_before.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_before.csv")

#---after
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "post-NACT",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$cell_subtype,
                      species = "Human",condition = d2$patient_id, evbiog = evbiog, risc = risc, if_normalize = T)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_after_filter.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_after_filter.csv")
obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_after.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_after.csv")


#---paired samples
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
samplenames <- unique(d2$patient_id)
res_all <- data.frame()
for (i in 1:length(samplenames)) {
    a2 <- d2[d2$patient_id == samplenames[i],]
    a1 <- d1[,a2$cell]
    obj <- create_miRTalk(sc_data = a1, sc_celltype = a2$cell_subtype,
                          species = "Human",condition = a2$treatment_phase, evbiog = evbiog, risc = risc, if_normalize = T)
    obj <- find_hvtg(object = obj)
    obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
    obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
    cci <- obj@cci
    cci$sample <- samplenames[i]
    res_all <- rbind(res_all, cci)
}
saveRDS(res_all, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_paired_samples_filter.rds")
write.csv(res_all, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_paired_samples_filter.csv")

res_all <- data.frame()
for (i in 1:length(samplenames)) {
  a2 <- d2[d2$patient_id == samplenames[i],]
  a1 <- d1[,a2$cell]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = a2$cell_subtype,
                        species = "Human",condition = a2$treatment_phase, evbiog = evbiog, risc = risc, if_normalize = T)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
  cci <- obj@cci
  cci$sample <- samplenames[i]
  res_all <- rbind(res_all, cci)
}
saveRDS(res_all, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_paired_samples.rds")
write.csv(res_all, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_paired_samples.csv")



#---samples-before
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "treatment-naive",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$patient_id)
res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$patient == sampleid[i],]
  # 1.0
  print(paste0("ovarian_",i,"_1.0"))
  a3 <- a2[,c("cell","cell_subtype")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]),
                        species = "Human", condition = rep(sampleid[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = T)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
saveRDS(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_before_filter.rds")
write.csv(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_before_filter.csv")

res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$patient == sampleid[i],]
  # 1.0
  print(paste0("ovarian_",i,"_1.0"))
  a3 <- a2[,c("cell","cell_subtype")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]),
                        species = "Human", condition = rep(sampleid[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = T)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
saveRDS(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_before.rds")
write.csv(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_before.csv")

#---samples-after
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "post-NACT",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$patient_id)
res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$patient == sampleid[i],]
  # 1.0
  print(paste0("ovarian_",i,"_1.0"))
  a3 <- a2[,c("cell","cell_subtype")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]),
                        species = "Human", condition = rep(sampleid[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = T)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
saveRDS(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_after_filter.rds")
write.csv(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_after_filter.csv")

res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$patient == sampleid[i],]
  # 1.0
  print(paste0("ovarian_",i,"_1.0"))
  a3 <- a2[,c("cell","cell_subtype")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]),
                        species = "Human", condition = rep(sampleid[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = T)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
saveRDS(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_after.rds")
write.csv(res_tmp, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_samples_after.csv")
