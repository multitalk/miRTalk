devtools::load_all('~/github/miRTalk/')

# GBM cancer
setwd('~/workspace/miRTalk/test_data/glioblastoma/GSE84465_OK/')
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

obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
saveRDS(obj, file = "revise/obj.rds")

obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
saveRDS(obj, file = "revise/obj_filter.rds")

# kidney
setwd('~/workspace/miRTalk/cases/kidney_fibrosis/')
load("rawdata.rda")
sc_data <- rev_gene(sc_data, data_type = 'count',species = 'Mouse',geneinfo = geneinfo)
sc_meta[sc_meta$Time.point != "Uninjured",]$Time.point <- "Injured"
obj <- create_miRTalk(sc_data = sc_data, sc_celltype = sc_meta$Annotation.Level.2,
                      species = "Mouse", condition = sc_meta$Time.point, evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
saveRDS(obj, file = "revise/obj_filter.rds")

obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
saveRDS(obj, file = "revise/obj.rds")

#---samples
a2 <- sc_meta[sc_meta$Time.point == "Uninjured",]
a1 <- sc_data[,a2$id]
obj <- create_miRTalk(sc_data = a1, sc_celltype = a2$Annotation.Level.2,
                      species = "Mouse", condition = a2$Time.point, evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
cci <- obj@cci
saveRDS(obj, file = "revise/obj_uninjured_filter.rds")

obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
cci <- obj@cci
saveRDS(obj, file = "revise/obj_uninjured.rds")

a2 <- sc_meta[sc_meta$Time.point == "Injured",]
a1 <- sc_data[,a2$id]
obj <- create_miRTalk(sc_data = a1, sc_celltype = a2$Annotation.Level.2,
                      species = "Mouse", condition = a2$Time.point, evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
saveRDS(obj, file = "revise/obj_injured_filter.rds")
obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
saveRDS(obj, file = "revise/obj_injured.rds")

# Liver
setwd('~/workspace/miRTalk/test_data/liver/TLMA-SC-A031')
load('sc_data_raw.rda')
load('sc_meta.rda')
load('sc_meta_specimen.rda')
sc_data_raw <- rev_gene(sc_data_raw,data_type = 'count',species = 'Rat',geneinfo = geneinfo)

res_condition <- list()
samplename <- unique(sc_meta$Specimen)
for (i in 1:length(samplename)) {
  print(i)
  a2 <- sc_meta[sc_meta$Specimen == samplename, ]
  a1 <- sc_data_raw[,a2$Id]
  obj <- create_miRTalk(sc_data = a1,sc_celltype = a2$Name_Celltype,
                        species = "Rat",condition = rep(samplename[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = F)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
  res_condition[[i]] <- obj
}
names(res_condition) <- samplename
saveRDS(res_condition, file = '../../../cases/rat_LT/obj_samples_filter.rds')

res_condition <- list()
samplename <- unique(sc_meta$Specimen)
for (i in 1:length(samplename)) {
  print(i)
  a2 <- sc_meta[sc_meta$Specimen == samplename, ]
  a1 <- sc_data_raw[,a2$Id]
  obj <- create_miRTalk(sc_data = a1,sc_celltype = a2$Name_Celltype,
                        species = "Rat",condition = rep(samplename[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = F)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
  res_condition[[i]] <- obj
}
names(res_condition) <- samplename
saveRDS(res_condition, file = '../../../cases/rat_LT/obj_samples.rds')

sc_meta$condition <- "FDL"
sc_meta[sc_meta$Specimen %in% c("CDL1", "CDL2","CDL3"),]$condition <- "CDL"
res_condition <- list()
samplename <- unique(sc_meta$condition)
for (i in 1:length(samplename)) {
  print(i)
  a2 <- sc_meta[sc_meta$condition == samplename[i], ]
  a1 <- sc_data_raw[,a2$Id]
  obj <- create_miRTalk(sc_data = a1,sc_celltype = a2$Name_Celltype,
                        species = "Rat",condition = rep(samplename[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = F)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
  res_condition[[i]] <- obj
}
names(res_condition) <- samplename
saveRDS(res_condition, file = '../../../cases/rat_LT/obj_condition_filter.rds')

res_condition <- list()
samplename <- unique(sc_meta$condition)
for (i in 1:length(samplename)) {
  print(i)
  a2 <- sc_meta[sc_meta$condition == samplename[i], ]
  a1 <- sc_data_raw[,a2$Id]
  obj <- create_miRTalk(sc_data = a1,sc_celltype = a2$Name_Celltype,
                        species = "Rat",condition = rep(samplename[i], ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = F)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
  res_condition[[i]] <- obj
}
names(res_condition) <- samplename
saveRDS(res_condition, file = '../../../cases/rat_LT/obj_condition.rds')

obj <- create_miRTalk(sc_data = sc_data_raw, sc_celltype = sc_meta$Name_Celltype,
                      species = "Rat", condition = sc_meta$condition, evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = T, use_n_cores = 16)
cci <- obj@cci
saveRDS(obj, file = "../../../cases/rat_LT/obj_filter.rds")

obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
cci <- obj@cci
saveRDS(obj, file = "../../../cases/rat_LT/obj.rds")
