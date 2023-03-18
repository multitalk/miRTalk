sc_data_raw <- rev_gene(sc_data_raw,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sc_meta1 <- sc_meta[sc_meta$Specimen %in% c("Pt14b"),]
sc_data_raw1 <- sc_data_raw[,sc_meta1$Id]

obj <- create_miRTalk(sc_data = sc_data_raw1, sc_celltype = sc_meta1$Name_Celltype, species = "Human")
obj <- find_cci(obj, mir_info = mir_info, mir2tar = mir2tar)

res <- list()
# mouse
sc_data_raw <- rev_gene(sc_data_raw,data_type = 'count',species = 'Mouse',geneinfo = geneinfo)
sampleid <- unique(sc_meta_specimen$Specimen)
for (i in 1:length(sampleid)) {
  print(i)
  d2 <- sc_meta[sc_meta$Specimen == sampleid[i], ]
  d1 <- sc_data_raw[,d2$Id]
  obj <- create_miRTalk(sc_data = d1, sc_celltype = d2$Name_Celltype, species = "Mouse")
  obj <- find_cci(obj, mir_info = mir_info, mir2tar = mir2tar)
  res[[i]] <- obj
}


