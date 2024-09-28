devtools::load_all('~/github/miRTalk/')

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
obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
cci <- obj@cci
saveRDS(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_all.rds")
write.csv(cci, file = "~/workspace/miRTalk/revise/benchmark/cci_ov_all.csv")
