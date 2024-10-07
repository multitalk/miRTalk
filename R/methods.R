#' @title Create miRTalk object
#'
#' @description create miRTalk object using single-cell transcriptomics data
#' @param sc_data A data.frame or matrix or dgCMatrix containing raw counts of single-cell RNA-seq data. see \code{\link{demo_sc_data}}
#' @param sc_celltype A character containing the cell type of the single-cell RNA-seq data with the same length as the number of cells.
#' @param species A character meaning species of the single-cell transcriptomics data.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @param condition A character with the same length as the number of cells, e.g., control/disease/treatment, phase 1/2/3, men/women.
#' @param if_normalize Normalize sc_data with Seurat LogNormalize. Set it \code{FLASE} when sc_data has been normalized.
#' @param evbiog A data.frame of the system data containing extracellular vesicle biogenesis genes of \code{"Human"}, \code{"Mouse"}, and \code{"Rat"}.
#' @param risc A data.frame of the system data containing RNA-induced silencing complex related genes of \code{"Human"}, \code{"Mouse"}, and \code{"Rat"}.
#' @return miRTalk object
#' @import Matrix methods Seurat
#' @export

create_miRTalk <- function(sc_data, sc_celltype, species, condition, if_normalize = TRUE, evbiog, risc) {
    if (is(sc_data, "data.frame")) {
        sc_data <- .get_dgCMatrix(as.matrix(sc_data))
    }
    if (is(sc_data, "matrix")) {
        sc_data <- .get_dgCMatrix(as.matrix(sc_data))
    }
    if (!is(sc_data, "dgCMatrix")) {
        stop("sc_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (!is.character(sc_celltype)) {
        stop("sc_celltype is not a character!")
    }

    if (length(unique(sc_celltype)) == 1) {
      stop("Only one cell type. Please take it as bulk data and run find_miRTalk_bulk() instead!")
    }
    if (!is.character(condition)) {
        stop("condition is not a character!")
    }
    if (ncol(sc_data) != length(sc_celltype)) {
        stop("ncol(sc_data) is not consistent with length(sc_celltype)!")
    }
    if (ncol(sc_data) != length(condition)) {
        stop("ncol(sc_data) is not consistent with length(sc_celltype)!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    evbiog_genes <- evbiog[evbiog$species == species, ]$gene
    risc_genes <- risc[risc$species == species, ]$gene
    sc_celltype_new <- .rename_chr(sc_celltype)
    warning_info <- .show_warning(sc_celltype, sc_celltype_new)
    if (!is.null(warning_info)) {
        warning(warning_info)
    }
    colnames(sc_data) <- .rename_chr(colnames(sc_data))
    sc_meta <- data.frame(cell = colnames(sc_data), celltype = sc_celltype_new, condition = condition, stringsAsFactors = FALSE)
    sc_data <- sc_data[rowSums(sc_data) > 0, ]
    sc_data <- Seurat::CreateSeuratObject(sc_data)
    # if_normalize
    if (if_normalize) {
        sc_data <- Seurat::NormalizeData(sc_data,verbose = FALSE)
    }
    sc_data <- Seurat::AddModuleScore(sc_data, features = list(evbiog = evbiog_genes, risc = risc_genes))
    sc_meta$evbiog_score <- .minmax_normalize(sc_data@meta.data$Cluster1)
    sc_meta$risc_score <- .minmax_normalize(sc_data@meta.data$Cluster2)
    ver <- packageVersion("Seurat")
    ver <- substr(ver,1,1)
    if (ver >= 5) {
        sc_data <- sc_data@assays$RNA$data
    } else {
        sc_data <- sc_data[["RNA"]]@data
    }
    # generate miRTalk object
    object <- new("miRTalk", data = list(data = sc_data), meta = sc_meta, species = species)
    return(object)
}
#' @title Find highly variable target genes
#'
#' @description Find highly variable target genes with DEGs and HVGs
#' @param object miRTalk object after \code{\link{find_miRNA}}
#' @param pvalue Cutoff of p value. Default is \code{0.05}
#' @param log2fc log2 fold change for identifying the highly expressed genes in each cell type. Default is \code{0.5}
#' @param min_cell_num Min cell number for each cell type. Default is \code{10}
#' @param nfeatures Number of features to select as top variable features. Default is \code{3000}
#' @return miRTalk object containing highly variable target genes without the cell-type-specific potential marker genes
#' @import Matrix Seurat
#' @importFrom utils installed.packages
#' @export

find_hvtg <- function(object, pvalue = 0.05, log2fc = 0.5, min_cell_num = 10, nfeatures = 3000) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    # get var_genes
    sc_data <- object@data$data
    sc_meta <- object@meta
    sc_data <- Seurat::CreateSeuratObject(counts = sc_data)
    ver <- packageVersion("Seurat")
    ver <- substr(ver,1,1)
    if (ver >= 5) {
        sc_data@assays$RNA$data <- sc_data@assays$RNA$counts
    }
    sc_data <- Seurat::FindVariableFeatures(sc_data, verbose = FALSE, nfeatures = nfeatures)
    if (ver >= 5) {
        var_genes <- sc_data@assays[["RNA"]]@meta.data[["var.features"]]
    } else {
        var_genes <- sc_data@assays$RNA@var.features
    }
    Seurat::Idents(sc_data) <- sc_meta$celltype
    # select cell types with min_cell_num
    celltype_meta <- as.data.frame(table(sc_meta$celltype), stringsAsFactors = FALSE)
    colnames(celltype_meta) <- c("celltype", "cell_number")
    celltype_meta <- celltype_meta[celltype_meta$cell_number >= min_cell_num, ]
    if (nrow(celltype_meta) == 0) {
        stop("No cell types found with cells more than min_cell_num!")
    }
    if (nrow(celltype_meta) > 1) {
        marker_genes <- NULL
        for (i in 1:nrow(celltype_meta)) {
            marker_genes_tmp  <- .get_markers(sc_data, celltype_meta$celltype[i], celltype_meta[-i, ]$celltype, pvalue, log2fc)
            marker_genes <- c(marker_genes, marker_genes_tmp)
        }
        if (!is.null(marker_genes)) {
            marker_genes <- unique(marker_genes)
        }
        var_genes <- union(var_genes, marker_genes)
    }
    object@data$var_genes <- var_genes
    return(object)
}

#' @title Find expressed miRNAs
#'
#' @description Find expressed miRNAs among all cells and generate background distribution for permutation test
#' @param object miRTalk object after \code{\link{create_miRTalk}}
#' @param mir_info A data.frame of the system data containing information of EV-derived miRNA of \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_mir_info}}
#' @param mir2tar A data.frame of the system data containing miRNA-target interactions for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_mir2tar}}
#' @param min_percent Min percent of expressed cells for target genes of miRNA. Default is \code{0.05}
#' @param database Which database of miRNA-target interactions to use, "miRTarBase" and/or "TarBase". Default is the "miRTarBase". It can also be "TarBase" or c("miRTarBase", "TarBase")
#' @param resolution Correct to precursor or mature miRNAs. Use 'precursor' or 'mature'. Default is \code{'mature'}
#' @param regulation Inference of negative or positive regulation. Default is "negative". Set it as "positive" and set database as "TarBase" for inferring positive regulation.
#' @param EXOmotif A sequence called EXOmotif to help miRNA secretion in EVs such as "CAUG", "CGGGAG". Please refer to https://doi.org/10.1038/s41586-021-04234-3
#' @param if_use_human_data Whether to use homologous human data in mir_info and mir2tar for mouse or rat scRNA-seq data. For human scRNA-seq data, no need to do it. For mouse or rat data, you can set it TRUE.
#' @param if_combine Whether to use combined homologous mir_info and mir2tar when if_use_human_data is TRUE. Default is TRUE.
#' @param gene2gene A data.frame of the system data containing the gene orthologs among human, mouse, and rat. If if_use_human_data is TRUE, please provide it, like "gene2gene = gene2gene"
#' @param per_num Number of permutation test. Default is 1000
#' @return miRTalk object containing the expressed miRNAs
#' @import Matrix stringr
#' @export

find_miRNA <- function(object, mir_info, mir2tar, min_percent = 0.05, database = "miRTarBase", resolution = "mature", regulation = "negative",
                       EXOmotif = NULL, if_use_human_data = FALSE, if_combine = TRUE, gene2gene = NULL, per_num = 1000) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    if (!all(c("miRNA", "miRNA_mature", "gene", "species") %in% colnames(mir_info))) {
        stop("Please provide a correct mir_info data.frame! See demo_mir_info()!")
    }
    # check miRNA
    sc_data <- object@data$data
    sc_meta <- object@meta
    species <- object@species
    var_genes <- object@data$var_genes
    geneinfo_species <- geneinfo[geneinfo$species == species,]
    if (length(database) == 1) {
        if (database == "miRTarBase") {
            if (regulation != "negative") {
                stop("Please set database as TarBase for positive regulation!")
            }
            mir2tar <- mir2tar[grep(pattern = "miRTarBase", x = mir2tar$source),]
        } else {
            mir2tar <- mir2tar[grep(pattern = "TarBase-v9.0", x = mir2tar$source),]
            if (regulation != "negative") {
                mir2tar <- mir2tar[mir2tar$regulation == "Positive",]
            }
        }
    }
    object@type <- c(resolution, regulation)
    if (!is.null(EXOmotif)) {
        mir_info <- mir_info[grep(pattern = EXOmotif,x = mir_info$seq), ]
    }
    if (nrow(mir_info) == 0) {
        stop("No miRNAs found with EXOmotif!")
    }
    mir_list <- .use_human_data(mir_info, mir2tar, geneinfo_species, gene2gene, species, if_use_human_data, if_combine)
    mir_info <- mir_list[[1]]
    mir2tar <- mir_list[[2]]
    mir_gene <- unique(mir_info$gene)
    mir_gene <- mir_gene[mir_gene %in% rownames(sc_data)]
    if (length(mir_gene) == 0) {
        stop("No miRNAs found in sc_data!")
    }
    sc_data <- sc_data[mir_gene, ]
    if (length(mir_gene) == 1) {
        mir_gene_num <- .number_cell(sc_data)
        names(mir_gene_num) <- mir_gene
        mir_gene <- mir_gene_num
    } else {
        mir_gene <- base::apply(sc_data, 1, .number_cell)
    }
    mir_gene <- mir_gene[mir_gene > 0]
    if (length(mir_gene) == 0) {
        stop("No expressed miRNAs found in sc_data!")
    }
    mir_info <- mir_info[mir_info$gene %in% names(mir_gene), ]
    if (resolution == "mature") {
        mir2tar <- mir2tar[mir2tar$miRNA_mature %in% mir_info$miRNA_mature, ]
    } else {
        mir2tar <- mir2tar[mir2tar$miRNA %in% mir_info$miRNA, ]
    }
    sc_data <- object@data$data
    mir2tar <- mir2tar[mir2tar$target_gene %in% rownames(sc_data), ]
    if (nrow(mir2tar) == 0) {
        stop("No miRNA target genes found in sc_data!")
    }
    mir2tar$percent <- as.numeric(base::apply(X = sc_data[mir2tar$target_gene, ], MARGIN = 1, FUN = .percent_cell))
    mir2tar <- mir2tar[mir2tar$percent >= min_percent, ]
    if (nrow(mir2tar) == 0) {
        stop("No miRNA target genes with min_percent found in sc_data!")
    }
    mir2tar <- mir2tar[mir2tar$target_gene %in% var_genes, ]
    if (nrow(mir2tar) == 0) {
        stop("No miRNA target genes overlapped with var_genes!")
    }
    if (resolution == "mature") {
        mir_info <- mir_info[mir_info$miRNA_mature %in% mir2tar$miRNA_mature, ]
    } else {
        mir_info <- mir_info[mir_info$miRNA %in% mir2tar$miRNA, ]
    }
    mir2tar <- .get_miR2tar(mir_info, mir2tar, resolution)
    object@miR <- mir_info
    object@miR2tar <- mir2tar
    # condition_per
    condition <- unique(sc_meta$condition)
    per_test_list <- list()
    for (c in 1:length(condition)) {
        sc_meta_tmp <- sc_meta[sc_meta$condition == condition[c], ]
        per_test_list_tmp <- .get_per_test_list(sc_meta_tmp, per_num)
        per_test_list[[c]] <- per_test_list_tmp
    }
    names(per_test_list) <- condition
    object@per_test_list <- per_test_list
    return(object)
}

#' @title Infer cell-cell communications mediated by EV-derived miRNAs
#'
#' @description Infer cell-cell communications mediated by exosomal miRNAs from senders to receivers
#' @param object miRTalk object after \code{\link{create_miRTalk}}
#' @param min_cell_num Min cell number for each cell type and expressed miRNA. Default is \code{10}
#' @param min_percent Min percent of expressed cells for target genes of miRNA. Default is \code{0.05}
#' @param pvalue Cutoff of p value. Default is \code{0.05}
#' @param if_filter_miRNA Whether to filter the significantly highly expressed miRNAs. Default is FALSE
#' @param if_doParallel Use doParallel. Default is TRUE
#' @param use_n_cores Number of CPU cores to use. Default is 4
#' @return miRTalk object containing the inferred cell-cell communications mediated by EV-derived miRNAs
#' @import Matrix methods progress doParallel parallel foreach
#' @importFrom stats median wilcox.test sd
#' @importFrom correlation correlation
#' @export

find_miRTalk <- function(object, min_cell_num = 10, min_percent = 0.05, pvalue = 0.05, per_num = 1000, if_filter_miRNA = FALSE, if_doParallel = TRUE, use_n_cores = 4) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    if (!all(c("miRNA", "miRNA_mature", "target_gene", "species") %in% colnames(mir2tar))) {
        stop("Please provide a correct pathways mir2tar! See demo_mir2tar()!")
    }
    # check miRNA
    sc_data <- object@data$data
    sc_meta <- object@meta
    species <- object@species
    regulation <- object@type[2]
    condition <- unique(sc_meta$condition)
    miR_info <- object@miR
    miR2tar <- object@miR2tar
    per_test_list_all <- object@per_test_list
    if (nrow(miR_info) == 0) {
        stop("No expressed miRNA found! Please run find_miRNA(). Ignore this if you have run find_miRNA()")
    }
    cci_all <- data.frame()
    for (c in 1:length(condition)) {
        cat(paste0("[", condition[c], "]"),"\n")
        per_test_list <- per_test_list_all[[condition[c]]]
        sc_meta_tmp <- sc_meta[sc_meta$condition == condition[c], ]
        sc_data_tmp <- sc_data[ ,sc_meta_tmp$cell]
        # select cell types with min_cell_num
        celltype_sender_meta <- as.data.frame(table(sc_meta_tmp$celltype), stringsAsFactors = FALSE)
        colnames(celltype_sender_meta) <- c("celltype_sender", "cell_number")
        rownames(celltype_sender_meta) <- celltype_sender_meta$celltype_sender
        celltype_sender_meta <- celltype_sender_meta[celltype_sender_meta$cell_number >= min_cell_num, ]
        if (nrow(celltype_sender_meta) == 0) {
            warning(paste0("No cell types found with cells more than min_cell_num for the condition of ", condition[c], "!"))
            next
        }
        # generate pair-wise cell types
        cellname <- celltype_sender_meta$celltype_sender
        celltype_pair <- NULL
        for (i in 1:length(cellname)) {
            d1 <- data.frame(celltype_sender = rep(cellname[i], length(cellname)), celltype_receiver = cellname, stringsAsFactors = FALSE)
            celltype_pair <- rbind(celltype_pair, d1)
        }
        celltype_pair$condition <- condition[c]
        # perform cci
        if (if_doParallel) {
            n_cores <- parallel::detectCores()
            if (n_cores <= 4) {
                n_cores <- 1
            } else {
                n_cores <- use_n_cores
            }
            cl <- parallel::makeCluster(n_cores)
            doParallel::registerDoParallel(cl)
            all_res <- foreach::foreach(k = 1:nrow(celltype_pair), .packages = "Matrix", .export = c(".percent_cell", ".get_ndata_receiver",".number_cell",
                ".get_percent_cell", ".get_miR_gene_percent", ".get_cci", ".get_S", ".prepare_evbiog", ".get_evmirna_mean", ".get_evmirna_score",".get_sc_data_evbiog", ".get_sc_data_risc")) %dopar% {
                cci <- data.frame()
                # celltype_sender
                celltype_sender <- celltype_pair$celltype_sender[k]
                cell_sender <- sc_meta_tmp[sc_meta_tmp$celltype == celltype_sender, ]
                # select miR_gene with min_expressed_percent in senders
                miR_gene <- unique(miR2tar$gene)
                miR_gene <- data.frame(miR_gene = miR_gene, stringsAsFactors = FALSE)
                ndata_sender <- sc_data_tmp[miR_gene$miR_gene, cell_sender$cell]
                ndata_sender_evbiog <- .prepare_evbiog(ndata_sender, cell_sender, miR_gene) # added
                if (nrow(miR_gene) == 1) {
                    miR_gene$percent <- .number_cell(ndata_sender)
                    miR_gene$EVmiR_score <- as.numeric(.get_evmirna_mean(ndata_sender_evbiog)) # added
                } else {
                    miR_gene$percent <- as.numeric(apply(ndata_sender, 1, .number_cell))
                    miR_gene$EVmiR_score <- as.numeric(apply(ndata_sender_evbiog, 1, .get_evmirna_mean)) # added
                }
                miR_gene <- miR_gene[miR_gene$percent >= min_cell_num, ]
                if (nrow(miR_gene) > 0) {
                    miR_gene <- .get_evmirna_score(miR_gene, ndata_sender, sc_data_tmp)
                    # celltype_receiver
                    celltype_receiver <- celltype_pair$celltype_receiver[k]
                    cell_receiver <- sc_meta_tmp[sc_meta_tmp$celltype == celltype_receiver, ]
                    # infer miRNA
                    miR_gene <- .get_miR_gene_percent(miR_gene, cell_sender, cell_receiver, sc_data_tmp, pvalue, if_filter_miRNA) # revise
                    if (nrow(miR_gene) > 0) {
                        miR2tar_sub <- miR2tar[miR2tar$gene %in% miR_gene$miR_gene, ]
                        miR_gene <- miR_gene[miR_gene$miR_gene %in% miR2tar_sub$gene,]
                        miR_name <- unique(miR2tar_sub$miRNA)
                        # ndata_receiver
                        ndata_receiver <- sc_data_tmp[, cell_receiver$cell]
                        ndata_other <- sc_data_tmp[unique(miR2tar_sub$target_gene), !colnames(sc_data_tmp) %in% colnames(ndata_receiver)]
                        ndata_receiver_mean <- .get_ndata_receiver(ndata_receiver)
                        ndata_receiver_mean_median <- stats::median(ndata_receiver_mean$activity)
                        ndata_receiver_mean <- ndata_receiver_mean[unique(miR2tar_sub$target_gene), ]
                        ndata_receiver <- ndata_receiver[unique(miR2tar_sub$target_gene), ]
                        ndata_receiver_mean <- .get_percent_cell(ndata_receiver_mean, ndata_receiver, ndata_other, pvalue, regulation)
                        if (min(ndata_receiver_mean$target_gene_percent) < min_percent) {
                            ndata_receiver_mean[ndata_receiver_mean$target_gene_percent < min_percent, ]$sig <- "NO"
                        }
                        cci_temp <- data.frame()
                        # get sc_data_evbiog and sc_data_risc # added
                        sc_data_evbiog <- .get_sc_data_evbiog(sc_data_tmp, sc_meta_tmp, miR2tar_sub) # added
                        sc_data_risc <- .get_sc_data_risc(sc_data_tmp, sc_meta_tmp, miR2tar_sub, regulation) # added
                        for (i in 1:length(miR_name)) {
                            miR_gene2tar <- miR2tar_sub[miR2tar_sub$miRNA == miR_name[i], ]
                            miR_gene2tar <- base::merge(miR_gene2tar, ndata_receiver_mean)
                            miR_gene2tar <- .get_S(miR_gene2tar, sc_data_evbiog, sc_data_risc, cell_sender, cell_receiver, per_test_list) # revise
                            miR_gene2tar <- miR_gene2tar[miR_gene2tar$sig == "YES", ]
                            if (nrow(miR_gene2tar) > 0) {
                                miR_gene2tar <- miR_gene2tar[miR_gene2tar$prob < pvalue, ] # revise
                                if (nrow(miR_gene2tar) > 0) { # revise
                                    cci_temp <- rbind(cci_temp, miR_gene2tar) # revise
                                }
                            }
                        }
                        if (nrow(cci_temp) > 0) {
                            cci_temp <- .get_cci(cci_temp, celltype_sender, celltype_receiver, miR_gene)
                            cci <- rbind(cci, cci_temp)
                        }
                    }
                }
                list(cci = cci)
            }
            doParallel::stopImplicitCluster()
            parallel::stopCluster(cl)
            cci <- data.frame()
            for (i in 1:length(all_res)) {
                cci_temp <- all_res[[i]]$cci
                if (nrow(cci_temp) > 0) {
                    cci <- rbind(cci, cci_temp)
                }
            }
        } else {
            pb <- progress::progress_bar$new(format = "[:bar] Finished::percent time::elapsedfull",
                total = nrow(celltype_pair), clear = FALSE, width = 60, complete = "+", incomplete = "-")
            cci <- data.frame()
            for (k in 1:nrow(celltype_pair)) {
                pb$tick()
                # celltype_sender
                celltype_sender <- celltype_pair$celltype_sender[k]
                cell_sender <- sc_meta_tmp[sc_meta_tmp$celltype == celltype_sender, ]
                # select miR_gene with min_expressed_percent in senders
                miR_gene <- unique(miR2tar$gene)
                miR_gene <- data.frame(miR_gene = miR_gene, stringsAsFactors = FALSE)
                ndata_sender <- sc_data_tmp[miR_gene$miR_gene, cell_sender$cell]
                ndata_sender_evbiog <- .prepare_evbiog(ndata_sender, cell_sender, miR_gene) # added
                if (nrow(miR_gene) == 1) { # revise
                    miR_gene$percent <- .number_cell(ndata_sender)
                    miR_gene$EVmiR_score <- as.numeric(.get_evmirna_mean(ndata_sender_evbiog)) # added
                } else {
                    miR_gene$percent <- as.numeric(apply(ndata_sender, 1, .number_cell))
                    miR_gene$EVmiR_score <- as.numeric(apply(ndata_sender_evbiog, 1, .get_evmirna_mean)) # added
                }
                miR_gene <- miR_gene[miR_gene$percent >= min_cell_num, ]
                if (nrow(miR_gene) == 0) {
                    next
                }
                miR_gene <- .get_evmirna_score(miR_gene, ndata_sender, sc_data_tmp)
                # celltype_receiver
                celltype_receiver <- celltype_pair$celltype_receiver[k]
                cell_receiver <- sc_meta_tmp[sc_meta_tmp$celltype == celltype_receiver, ]
                # infer miRNA
                miR_gene <- .get_miR_gene_percent(miR_gene, cell_sender, cell_receiver, sc_data_tmp, pvalue, if_filter_miRNA) # revise
                if (nrow(miR_gene) == 0) { # added
                    next # added
                } # added
                miR2tar_sub <- miR2tar[miR2tar$gene %in% miR_gene$miR_gene, ]
                miR_gene <- miR_gene[miR_gene$miR_gene %in% miR2tar_sub$gene,]
                miR_name <- unique(miR2tar_sub$miRNA)
                # ndata_receiver
                ndata_receiver <- sc_data_tmp[, cell_receiver$cell]
                ndata_other <- sc_data_tmp[unique(miR2tar_sub$target_gene), !colnames(sc_data_tmp) %in% colnames(ndata_receiver)]
                ndata_receiver_mean <- .get_ndata_receiver(ndata_receiver)
                ndata_receiver_mean <- ndata_receiver_mean[unique(miR2tar_sub$target_gene), ]
                ndata_receiver <- ndata_receiver[unique(miR2tar_sub$target_gene), ]
                ndata_receiver_mean <- .get_percent_cell(ndata_receiver_mean, ndata_receiver, ndata_other, pvalue, regulation)
                if (min(ndata_receiver_mean$target_gene_percent) < min_percent) {
                    ndata_receiver_mean[ndata_receiver_mean$target_gene_percent < min_percent, ]$sig <- "NO"
                }
                cci_temp <- data.frame()
                # get sc_data_evbiog and sc_data_risc # added
                sc_data_evbiog <- .get_sc_data_evbiog(sc_data_tmp, sc_meta_tmp, miR2tar_sub) # added
                sc_data_risc <- .get_sc_data_risc(sc_data_tmp, sc_meta_tmp, miR2tar_sub, regulation) # added
                for (i in 1:length(miR_name)) {
                    miR_gene2tar <- miR2tar_sub[miR2tar_sub$miRNA == miR_name[i], ]
                    miR_gene2tar <- base::merge(miR_gene2tar, ndata_receiver_mean)
                    miR_gene2tar <- .get_S(miR_gene2tar, sc_data_evbiog, sc_data_risc, cell_sender, cell_receiver, per_test_list) # revise
                    miR_gene2tar <- miR_gene2tar[miR_gene2tar$sig == "YES", ]
                    if (nrow(miR_gene2tar) > 0) {
                        miR_gene2tar <- miR_gene2tar[miR_gene2tar$prob < pvalue, ] #revise
                        if (nrow(miR_gene2tar) > 0) { # revise
                            cci_temp <- rbind(cci_temp, miR_gene2tar) # revise
                        }
                    }
                }
                if (nrow(cci_temp) > 0) {
                    cci_temp <- .get_cci(cci_temp, celltype_sender, celltype_receiver, miR_gene)
                    cci <- rbind(cci, cci_temp)
                }
            }
        }
        if (nrow(cci) > 0) {
            cci$condition <- condition[c]
            cci_all <-rbind(cci_all, cci)
        }
    }
    object@cci <- cci_all
    return(object)
}

#' @title Infer EV-derived miR-target interactions for bulk data
#'
#' @description Infer EV-derived miR-target interactions for paired bulk RNA-seq and miRNA-seq data
#' @param rna_data RNA-seq data with rows containg genes and columns containing samples.
#' @param mirna_data miRNA-seq data with rows containg miRNA and columns containing samples.
#' @param type Which types of mirna_data, miRNA genes, precursor or mature miRNAs. Use 'gene', 'precursor' or 'mature'.
#' @param resolution Correct to precursor or mature miRNAs. Use 'precursor' or 'mature'. Default is \code{'mature'}
#' @param species A character meaning species of the single-cell transcriptomics data.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @param mir_info A data.frame of the system data containing information of EV-derived miRNA of \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_mir_info}}
#' @param mir2tar A data.frame of the system data containing relationship of miRNA and its target genes for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_mir2tar}}
#' @param if_normalize Normalize sc_data with Seurat LogNormalize. Set it \code{FLASE} when rna_data and mirna_data have been normalized.
#' @param if_use_evbiog_risc if considering module score of extracellular vesicle biogenesis genes and RNA-induced silencing complex related genes. Default is TRUE. Consider set it FALSE for the comparion of scores between various conditions.
#' @param evbiog A data.frame of the system data containing extracellular vesicle biogenesis genes of \code{"Human"}, \code{"Mouse"}, and \code{"Rat"}.
#' @param risc A data.frame of the system data containing RNA-induced silencing complex related genes of \code{"Human"}, \code{"Mouse"}, and \code{"Rat"}.
#' @param score_scale_method Methods for scale the Seurat scores of evbiog and risc signatures, "1" for min max scale used in scRNA-seq data by default, and "2" for rank scale. For small sample size, set it to "2" to reduce the zero scores.
#' @param target_scale_method Methods for scale the target gene expression, "1" for rank scale used in scRNA-seq data by default, "2" for the scale with values divided by the max value for each sample, "3" for min max scale. For bulk RNA-seq, consider set it to "2" (all values >= 0) or "3" to reduce the significant heterogeneity of samples, especially human samples.
#' @param database Which database of miRNA-target interactions to use, "miRTarBase" and/or "TarBase". Default is the "miRTarBase". It can also be "TarBase" or c("miRTarBase", "TarBase")
#' @param regulation Inference of negative or positive regulation. Default is "negative". Set it as "positive" and set database as "TarBase" for inferring positive regulation.
#' @param if_use_human_data Whether to use homologous human data in mir_info and mir2tar for mouse or rat scRNA-seq data. For human scRNA-seq data, no need to do it. For mouse or rat data, you can set it TRUE.
#' @param if_combine Whether to use combined homologous mir_info and mir2tar when if_use_human_data is TRUE. Default is TRUE.
#' @param gene2gene A data.frame of the system data containing the gene orthologs among human, mouse, and rat. If if_use_human_data is TRUE, please provide it, like "gene2gene = gene2gene"
#' @return A data.frame containing score for each miRNA-target interaction acoss samples
#' @import Matrix methods Seurat progress
#' @export

find_miRTalk_bulk <- function(rna_data, mirna_data, type, resolution = "mature", species,
    mir_info, mir2tar, if_normalize = TRUE, if_use_evbiog_risc = TRUE, evbiog = NULL, risc = NULL, score_scale_method = "1", target_scale_method = "1",
    database = "miRTarBase", regulation = "negative", if_use_human_data = FALSE, if_combine = TRUE, gene2gene = NULL){
    if (is(rna_data, "data.frame")) {
        rna_data <- .get_dgCMatrix(as.matrix(rna_data))
    }
    if (is(rna_data, "matrix")) {
        rna_data <- .get_dgCMatrix(as.matrix(rna_data))
    }
    if (!is(rna_data, "dgCMatrix")) {
        stop("rna_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (is(mirna_data, "data.frame")) {
        mirna_data <- .get_dgCMatrix(as.matrix(mirna_data))
    }
    if (is(mirna_data, "matrix")) {
        mirna_data <- .get_dgCMatrix(as.matrix(mirna_data))
    }
    if (!is(mirna_data, "dgCMatrix")) {
        stop("mirna_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (!ncol(rna_data) == ncol(mirna_data)) {
        stop("ncol(rna_data) is not consistent with ncol(mirna_data)!")
    }
    if (length(type) > 1 | !type %in% c("gene", "precursor", "mature")) {
        stop("Please provide a correct type, i.e., 'gene', 'precursor' or 'mature'!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    if(if_use_evbiog_risc){
        if (is.null(evbiog)) {
            stop("Please provide the data containing extracellular vesicle biogenesis genes!")
        }
        if (is.null(risc)) {
            stop("Please provide the data containing RNA-induced silencing complex related genes!")
        }
        evbiog_genes <- evbiog[evbiog$species == species, ]$gene
        risc_genes <- risc[risc$species == species, ]$gene
    }
    sample_meta <- data.frame(sample = colnames(rna_data), stringsAsFactors = F)
    rna_data <- rna_data[which(rowSums(rna_data) > 0),]
    rna_data <- Seurat::CreateSeuratObject(rna_data)
    mirna_data <- mirna_data[which(rowSums(mirna_data) > 0),]
    mirna_data <- Seurat::CreateSeuratObject(mirna_data)
    # if_normalize
    if (if_normalize) {
        rna_data <- Seurat::NormalizeData(rna_data,verbose = FALSE)
        mirna_data <- Seurat::NormalizeData(mirna_data,verbose = FALSE)
    }
    if(if_use_evbiog_risc){
        rna_data <- Seurat::AddModuleScore(rna_data, features = list(evbiog = evbiog_genes, risc = risc_genes))
        if (scale_method == "1") {
            sample_meta$evbiog_score <- .minmax_normalize(rna_data$Cluster1)
            sample_meta$risc_score <- .minmax_normalize(rna_data$Cluster2)
        }
        if (scale_method == "2") {
            sample_meta$evbiog_score <- .rank_normalize(rna_data$Cluster1)
            sample_meta$risc_score <- .rank_normalize(rna_data$Cluster2)
        }
    }
    ver <- packageVersion("Seurat")
    ver <- substr(ver,1,1)
    if (ver >= 5) {
        genenames <- rownames(rna_data)
        cellnames <- colnames(rna_data)
        rna_data <- rna_data[["RNA"]]@layers$data
        rownames(rna_data) <- genenames
        colnames(rna_data) <- cellnames
        genenames <- rownames(mirna_data)
        cellnames <- colnames(mirna_data)
        mirna_data <- mirna_data[["RNA"]]@layers$data
        rownames(mirna_data) <- genenames
        colnames(mirna_data) <- cellnames
    } else {
        rna_data <- rna_data[["RNA"]]@data
        mirna_data <- mirna_data[["RNA"]]@data
    }
    # map
    mir_list <- .use_human_data(mir_info, mir2tar, geneinfo_species, gene2gene, species, if_use_human_data, if_combine)
    mir_info <- mir_list[[1]]
    if (type == "gene") {
        mir_info <- mir_info[mir_info$gene %in% rownames(mirna_data),]
    }
    if (type == "precursor") {
        mir_info <- mir_info[mir_info$miRNA %in% rownames(mirna_data),]
    }
    if (type == "mature") {
        mir_info <- mir_info[mir_info$miRNA_mature %in% rownames(mirna_data),]
    }
    if (nrow(mir_info) == 0) {
        stop("No miRNA found in mirna_data!")
    }
    mir2tar <- mir_list[[2]]
    mir2tar <- mir2tar[ ,c(1,3,4)]
    if (resolution == "mature") {
        mir2tar <- mir2tar[mir2tar$miRNA_mature %in% mir_info$miRNA_mature, ]
        if (nrow(mir2tar) == 0) {
            stop("No mature miRNA found!")
        }
        mir2tar <- mir2tar[,-1]
        mir2tar <- merge(mir2tar, mir_info)
    } else {
        mir2tar <- mir2tar[mir2tar$miRNA %in% mir_info$miRNA, ]
        if (nrow(mir2tar) == 0) {
            stop("No precursor miRNA found!")
        }
        mir2tar <- mir2tar[,-2]
        mir2tar <- merge(mir2tar, mir_info)
    }
    mir2tar <- mir2tar[mir2tar$target_gene %in% rownames(rna_data), ]
    if (nrow(mir2tar) == 0) {
        stop("No miRNA target genes found in rna_data!")
    }
    if (regulation == "negative") {
        if (target_scale_method == "1") {
            rna_target <- apply(rna_data, 2, function(x){
                x_rank <- rank(-x)
                x_rank <- x_rank/max(x_rank)
                return(x_rank)
            })
        }
        if (target_scale_method == "2") {
            rna_target <- apply(rna_data, 2, function(x){
                x_scale <- x/max(x)
                x_scale <- 1-x_scale
                return(x_scale)
            })
        }
        if (target_scale_method == "3") {
            rna_target <- apply(rna_data, 2, function(x){
                x_minmax <- (x-min(x))/(max(x)-min(x))
                x_minmax <- 1-x_minmax
                return(x_minmax)
            })
        }
    } else {
        if (target_scale_method == "1") {
            rna_target <- apply(rna_data, 2, function(x){
                x_rank <- rank(x)
                x_rank <- x_rank/max(x_rank)
                return(x_rank)
            })
        }
        if (target_scale_method == "2") {
            rna_target <- apply(rna_data, 2, function(x){
                x_scale <- x/max(x)
                return(x_scale)
            })
        }
        if (target_scale_method == "3") {
            rna_target <- apply(rna_data, 2, function(x){
                x_minmax <- (x-min(x))/(max(x)-min(x))
                return(x_minmax)
            })
        }
    }
    res_list <- list()
    if (if_use_evbiog_risc) {
        evbiog_score <- sample_meta$evbiog_score
        risc_score <- sample_meta$risc_score
    }
    pb <- progress::progress_bar$new(format = "[:bar] Finished::percent time::elapsedfull",
        total = nrow(mir2tar), clear = FALSE, width = 60, complete = "+", incomplete = "-")
    for (i in 1:nrow(mir2tar)) {
        pb$tick()
        if (type == "gene") {
            mirname_exp <- mirna_data[mir2tar$gene[i], ]
        }
        if (type == "precursor") {
            mirname_exp <- mirna_data[mir2tar$miRNA[i], ]
        }
        if (type == "mature") {
            mirname_exp <- mirna_data[mir2tar$miRNA_mature[i], ]
        }
        targetgene_exp <- rna_target[mir2tar$target_gene[i], ]
        if (if_use_evbiog_risc) {
            sscore_mir2gene <- evbiog_score*mirname_exp*targetgene_exp*risc_score
        } else {
            sscore_mir2gene <- mirname_exp*targetgene_exp
        }
        res_list[[i]] <- sscore_mir2gene
    }
    res_list <- as.data.frame(t(as.data.frame(res_list)))
    rownames(res_list) <- 1:nrow(res_list)
    colnames(res_list) <- colnames(rna_data)
    rownames(mir2tar) <- 1:nrow(mir2tar)
    mir2tar <- cbind(mir2tar, res_list)
    return(mir2tar)
}

#' @title Get miRNA-target interactions
#'
#' @description Get simple results of miRNA-target interactions and specificity.
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param simple Whether to show the simple results. Default is \code{TRUE}
#' @return A data.frame containing all miRNA-target interactions.
#' @import Matrix
#' @export

get_miRTalk_cci <- function(object, simple = TRUE) {
    # check object
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (simple) {
        cci <- cci[,c(1:3, 7:8, 10, 17, 19)]
        condition <- unique(cci$condition)
        cci_all <- data.frame()
        for (c in 1:length(condition)) {
            cci_tmp <- cci[cci$condition == condition[c], ]
            cci_tmp <- .get_specifity(cci_tmp)
            cci_all <- rbind(cci_all, cci_tmp)
        }
    }
    return(cci_all)
}

#' @title Get circulating score of inferred miRNAs
#'
#' @description Get circulating score of inferred miRNAs.
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @return A data.frame containing all potential circulating miRNAs.
#' @import Matrix
#' @export

get_miRTalk_circulating_score <- function(object) {
    # check object
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    miR_data <- object@miR
    miR_data <- miR_data[which(miR_data$circulating_miRNA != "NA"), ]
    if (nrow(miR_data) > 0) {
        miR_data$avg_rpm <- log2(miR_data$avg_rpm+1)
        miR_data$score <- miR_data$avg_rpm/max(miR_data$avg_rpm)
        type <- object@type[1]
        if (type == "mature") {
            miR_data <- unique(miR_data[,c("miRNA_mature", "tissue_TarBase", "score")])
            colnames(miR_data)[1] <- "miRNA"
        } else {
            miR_data <- unique(miR_data[,c("miRNA", "tissue_TarBase", "score")])
        }
        cci <- object@cci
        cci <- unique(cci[,c("miRNA", "miR_gene", "celltype_receiver","target_gene", "condition")])
        cci <- merge(miR_data, cci)
        return(cci)
    } else {
        print("No avaliable circlulating miRNAs!")
    }
}

#' @title Get overlapped pathways
#'
#' @description Get overlapped pathways between miRNA and target gene related pathways
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param gene2path A data.frame of the system data containing gene-related pathways from KEGG, Reactome, GO_BP, Wikipathways for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}.
#' @param mir2path A data.frame of the system data containing miRNA-related pathways from KEGG, Reactome, GO_BP, Wikipathways for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}.
#' @param miRNA which miRNAs to analyze. Default is all inferred miRNAs in senders.
#' @param targetgenes which targetgenes to analyze. Default is all inferred target genes in receivers.
#' @return A list of pathways for miRNAs and target genes.
#' @import Matrix
#' @export

get_miRTalk_pathway <- function(object, gene2path, mir2path, miRNA = NULL, targetgenes = NULL) {
    # check object
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    speices <- object@species
    cci <- object@cci
    cci_list <- list()
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(miRNA)) {
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (!is.null(targetgenes)) {
        cci <- cci[cci$target_gene %in% targetgenes, ]
    }
    species <- object@species
    type <- object@type[1]
    gene2path <- gene2path[gene2path$species == species,]
    miRNA_mature_other <- mir2path$miRNA
    if (species == "Mouse") {
        miRNA_mature_other <- stringr::str_replace(string = miRNA_mature_other, pattern = "hsa",replacement = "mmu")
    }
    if (species == "Rat"){
        miRNA_mature_other <- stringr::str_replace(string = miRNA_mature_other, pattern = "hsa",replacement = "rno")
    }
    if (type != "mature") {
        miRNA_mature_other <- stringr::str_remove(string = miRNA_mature_other, pattern = "-3p")
        miRNA_mature_other <- stringr::str_remove(string = miRNA_mature_other, pattern = "-5p")
    }
    mir2path$miRNA <- miRNA_mature_other
    cci <- unique(cci[,c("miRNA","target_gene")])
    colnames(mir2path)[2] <- "miRNA"
    cci1 <- merge(cci, mir2path)
    colnames(cci1)[3:5] <- paste0("miRNA_",colnames(cci1)[3:5])
    cci_list[[1]] <- cci1
    colnames(gene2path)[2] <- "target_gene"
    cci2 <- merge(cci, gene2path)
    cci2 <- cci2[,c(2,1,3,4,6)]
    colnames(cci2)[3:5] <- paste0("target_gene_",colnames(cci2)[3:5])
    cci_list[[2]] <- cci2
    cci2 <- cci2[,-1]
    cci1 <- merge(cci1, cci2)
    cci1 <- cci1[,c(2,1,3:ncol(cci1))]
    cci1 <- cci1[cci1$miRNA_term == cci1$target_gene_term, ]
    cci_list[[3]] <- cci1
    names(cci_list) <- c("miRNA-pathways","target-gene-pathways","overlapped")
    return(cci_list)
}
