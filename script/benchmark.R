library(miRTalk)
# bladder cancer
d1 <- Read10X_h5("BLCA_GSE130001_expression.h5")
d2 <- read.table('BLCA_GSE130001_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$Celltype..major.lineage.,
                      species = "Human",condition = d2$Sample, evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "cci_bladder.rds")

# chol
load("sc_data_raw.rda")
load("sc_meta.rda")
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
saveRDS(cci, file = "cci_chol.rds")

# ovarian
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
d2$condition <- paste0(d2$treatment_phase,"_", d2$patient_id)
obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$cell_subtype,
                      species = "Human",condition = d2$condition, evbiog = evbiog, risc = risc, if_normalize = T)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "cci_ov.rds")
