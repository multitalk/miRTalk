load("~/github/miRTalk/inst/extdata/example.rda")
devtools::load_all("~/github/miRTalk/")
obj <- create_miRTalk(sc_data = sc_data, sc_celltype = sc_celltype,
                      species = "Human",condition = rep("condition", length(sc_celltype)),
                      evbiog = evbiog,risc = risc)

obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,if_use_human_data = F,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_doParallel = T)

obj_cci <- get_miRTalk_cci(obj)

plot_miRTalk_chord(object = obj,show_type = "EVmiR_score")
plot_miRTalk_chord(object = obj,show_type = "score")
plot_miRTalk_circle(object = obj,show_type = "EVmiR_score")

plot_miRTalk_circle_simple(object = obj, celltype = "Tumor", show_type = "EVmiR_score")

plot_miRTalk_heatmap(object = obj,show_type = "EVmiR_score")
plot_miRTalk_sankey(object = obj,show_type = "EVmiR_score")
plot_miR_bubble(object = obj)
plot_miR_heatmap(object = obj)

plot_target_heatmap(object = obj,celltype = "Tumor")

plot_miR2tar_chord(object = obj, celltype_sender = "Tumor", celltype_receiver = "Bcell")
plot_miR2tar_circle(object = obj, celltype_sender = "Tumor", celltype_receiver = "Bcell")
plot_miR2tar_heatmap(object = obj, celltype_sender = "Tumor", celltype_receiver = "Bcell")

obj_filter <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)

object <- obj
min_cell_num = 10
min_percent = 0.05
pvalue = 0.05
per_num = 1000
if_consider_condition = TRUE
if_doParallel = FALSE
# > load("~/miRTalk/inst/extdata/example.rda")
#
# > obj <- create_miRTalk(sc_data = sc_data,
#                         sc_celltype = sc_celltype,
#                         species = "Human",
#                         condition = c(rep("Control", 250), rep("Disease", 265)),
#                         evbiog = evbiog,risc = risc)
#
# Warning: The following features are not present in the object: AGO2, not searching for symbol synonyms
#
# > obj <- find_miRNA(object = obj,mir_info = mir_info)
#
# > obj <- find_hvtg(object = obj)
#
# > obj <- find_miRTalk(obj, mir2tar = mir2tar, if_doParallel = F)
# [Control]
# [++++++++++++++++++++++++++++++] Finished:100% time:00:03:24
# [Disease]
# [++++++++++++++++++++++++++++++] Finished:100% time:00:02:24

