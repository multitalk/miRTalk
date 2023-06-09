#' @title Create miRTalk object
#'
#' @description create miRTalk object using single-cell transcriptomics data
#' @param sc_data A data.frame or matrix or dgCMatrix containing raw counts of single-cell RNA-seq data. see \code{\link{demo_sc_data}}
#' @param sc_celltype A character containing the cell type of the single-cell RNA-seq data
#' @param species A character meaning species of the single-cell transcriptomics data.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @param if_normalize Normalize sc_data with Seurat LogNormalize. Set it \code{FLASE} when sc_data has been normalized.
#' @return miRTalk object
#' @import Matrix methods Seurat
#' @export

create_miRTalk <- function(sc_data, sc_celltype, species, if_normalize = TRUE) {
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
    if (ncol(sc_data) != length(sc_celltype)) {
        stop("ncol(sc_data) is not consistent with length(sc_celltype)!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    sc_celltype_new <- .rename_chr(sc_celltype)
    warning_info <- .show_warning(sc_celltype, sc_celltype_new)
    if (!is.null(warning_info)) {
        warning(warning_info)
    }
    colnames(sc_data) <- .rename_chr(colnames(sc_data))
    sc_meta <- data.frame(cell = colnames(sc_data), celltype = sc_celltype_new, stringsAsFactors = FALSE)
    sc_data <- sc_data[rowSums(sc_data) > 0, ]
    # if_normalize
    if (if_normalize) {
        sc_data <- Seurat::CreateSeuratObject(sc_data)
        sc_data <- Seurat::NormalizeData(sc_data,verbose = FALSE)
        sc_data <- sc_data[["RNA"]]@data
    }
    # generate miRTalk object
    object <- new("miRTalk", data = list(data = sc_data), meta = sc_meta, species = species)
    return(object)
}

#' @title Find expressed miRNAs
#'
#' @description Find expressed miRNAs among all cells
#' @param object miRTalk object after \code{\link{create_miRTalk}}
#' @param mir_info A data.frame of the system data containing information of EV-derived miRNA of \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_mir_info}}
#' @return miRTalk object containing the expressed miRNAs
#' @import Matrix
#' @export

find_miRNA <- function(object, mir_info) {
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
    mir_info <- mir_info[mir_info$species == species, ]
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
        stop("No miRNAs with min_expressed_num found in sc_data!")
    }
    mir_info <- mir_info[mir_info$gene %in% names(mir_gene), ]
    object@miR <- mir_info
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
    sc_data <- Seurat::FindVariableFeatures(sc_data, verbose = FALSE, nfeatures = nfeatures)
    var_genes <- sc_data@assays$RNA@var.features
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

#' @title Infer cell-cell communications mediated by EV-derived miRNAs
#'
#' @description Infer cell-cell communications mediated by exosomal miRNAs from senders to receivers
#' @param object miRTalk object after \code{\link{create_miRTalk}}
#' @param mir2tar A data.frame of the system data containing relationship of miRNA and its target genes for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_mir2tar}}
#' @param min_cell_num Min cell number for each cell type and expressed miRNA. Default is \code{10}
#' @param pvalue Cutoff of p value. Default is \code{0.05}
#' @param resolution Correct to precursor or mature miRNAs. Use 'precursor' or 'mature'. Default is \code{'mature'}
#' @param min_percent Min percent of expressed cells for target genes of miRNA. Default is \code{0.05}
#' @param if_doParallel Use doParallel. Default is TRUE
#' @param use_n_cores Number of CPU cores to use. Default is 4
#' @return miRTalk object containing the inferred cell-cell communications mediated by EV-derived miRNAs
#' @import Matrix methods progress doParallel parallel foreach
#' @importFrom stats median wilcox.test sd
#' @importFrom correlation correlation
#' @export

find_miRTalk <- function(object, mir2tar, min_cell_num = 10, pvalue = 0.05, resolution = "mature", min_percent = 0.05, if_doParallel = TRUE, use_n_cores = 4) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    if (!all(c("miRNA", "miRNA_mature", "target_gene", "species") %in% colnames(mir2tar))) {
        stop("Please provide a correct pathways mir2tar! See demo_mir2tar()!")
    }
    # check miRNA
    sc_data <- object@data$data
    var_genes <- object@data$var_genes
    sc_meta <- object@meta
    species <- object@species
    mir_info <- object@miR
    if (nrow(mir_info) == 0) {
        stop("No expressed miRNA found! Please run find_miRNA(). Ignore this if you have run find_miRNA()")
    }
    mir2tar <- mir2tar[mir2tar$species == species, ]
    if (resolution == "mature") {
        mir2tar <- mir2tar[mir2tar$miRNA_mature %in% mir_info$miRNA_mature, ]
    } else {
        mir2tar <- mir2tar[mir2tar$miRNA %in% mir_info$miRNA, ]
    }
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
    object@miR <- mir_info
    object@miR2tar <- mir2tar
    # select cell types with min_cell_num
    celltype_sender_meta <- as.data.frame(table(sc_meta$celltype), stringsAsFactors = FALSE)
    colnames(celltype_sender_meta) <- c("celltype_sender", "cell_number")
    rownames(celltype_sender_meta) <- celltype_sender_meta$celltype_sender
    celltype_sender_meta <- celltype_sender_meta[celltype_sender_meta$cell_number >= min_cell_num, ]
    if (nrow(celltype_sender_meta) == 0) {
        stop("No cell types found with cells more than min_cell_num!")
    }
    # generate pair-wise cell types
    cellname <- celltype_sender_meta$celltype_sender
    celltype_pair <- NULL
    for (i in 1:length(cellname)) {
        d1 <- data.frame(celltype_sender = rep(cellname[i], length(cellname)), celltype_receiver = cellname, stringsAsFactors = FALSE)
        celltype_pair <- rbind(celltype_pair, d1)
    }
    miR2tar <- .get_miR2tar(object, resolution)
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
        all_res <- foreach::foreach(k = 1:nrow(celltype_pair), .packages = "Matrix", .export = c(".percent_cell", ".get_ndata_receiver",
            ".number_cell", ".get_percent_cell", ".get_miR_gene_percent", ".get_cci", ".get_P")) %dopar% {
            cci <- data.frame()
            # celltype_sender
            celltype_sender <- celltype_pair$celltype_sender[k]
            cell_sender <- sc_meta[sc_meta$celltype == celltype_sender, ]
            # select miR_gene with min_expressed_percent in senders
            miR_gene <- unique(miR2tar$gene)
            miR_gene <- data.frame(miR_gene = miR_gene, stringsAsFactors = FALSE)
            ndata_sender <- sc_data[miR_gene$miR_gene, cell_sender$cell]
            if (nrow(miR_gene) == 1) {
                miR_gene$percent <- .number_cell(ndata_sender)
            } else {
                miR_gene$percent <- as.numeric(apply(ndata_sender, 1, .number_cell))
            }
            miR_gene <- miR_gene[miR_gene$percent >= min_cell_num, ]
            if (nrow(miR_gene) > 0) {
                miR_gene$percent <- log10(miR_gene$percent)
                miR_gene$percent <- miR_gene$percent/max(miR_gene$percent)
                # celltype_receiver
                celltype_receiver <- celltype_pair$celltype_receiver[k]
                cell_receiver <- sc_meta[sc_meta$celltype == celltype_receiver, ]
                # infer miRNA
                miR_gene <- .get_miR_gene_percent(miR_gene, cell_sender, cell_receiver, sc_data)
                miR2tar_sub <- miR2tar[miR2tar$gene %in% miR_gene$miR_gene, ]
                miR_gene <- miR_gene[miR_gene$miR_gene %in% miR2tar_sub$gene,]
                miR_name <- unique(miR2tar_sub$miRNA)
                # ndata_receiver
                ndata_receiver <- sc_data[, cell_receiver$cell]
                ndata_other <- sc_data[unique(miR2tar_sub$target_gene), !colnames(sc_data) %in% colnames(ndata_receiver)]
                ndata_receiver_mean <- .get_ndata_receiver(ndata_receiver)
                ndata_receiver_mean_median <- stats::median(ndata_receiver_mean$activity)
                ndata_receiver_mean <- ndata_receiver_mean[unique(miR2tar_sub$target_gene), ]
                ndata_receiver <- ndata_receiver[unique(miR2tar_sub$target_gene), ]
                ndata_receiver_mean <- .get_percent_cell(ndata_receiver_mean, ndata_receiver, ndata_other, pvalue)
                if (min(ndata_receiver_mean$target_gene_percent) < min_percent) {
                    ndata_receiver_mean[ndata_receiver_mean$target_gene_percent < min_percent, ]$sig <- "NO"
                }
                cci_temp <- data.frame()
                for (i in 1:length(miR_name)) {
                    miR_gene2tar <- miR2tar_sub[miR2tar_sub$miRNA == miR_name[i], ]
                    miR_gene2tar <- base::merge(miR_gene2tar, ndata_receiver_mean)
                    miR_gene2tar <- .get_P(miR_gene2tar, sc_data, cell_sender, cell_receiver, pvalue)
                    miR_gene2tar <- miR_gene2tar[miR_gene2tar$sig == "YES", ]
                    if (nrow(miR_gene2tar) > 0) {
                        cci_temp <- rbind(cci_temp, miR_gene2tar)
                    }
                }
                if (nrow(cci_temp) > 0) {
                    cci_temp <- .get_cci(cci_temp, celltype_sender, celltype_receiver, miR_gene)
                    cci <- rbind(cci, cci_temp)
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
        cci <- data.frame()
        pb <- progress::progress_bar$new(format = "[:bar] Finished::percent time::elapsedfull",
            total = nrow(celltype_pair), clear = FALSE, width = 60, complete = "+", incomplete = "-")
        for (k in 1:nrow(celltype_pair)) {
            # celltype_sender
            celltype_sender <- celltype_pair$celltype_sender[k]
            cell_sender <- sc_meta[sc_meta$celltype == celltype_sender, ]
            # select miR_gene with min_expressed_percent in senders
            miR_gene <- unique(miR2tar$gene)
            miR_gene <- data.frame(miR_gene = miR_gene, stringsAsFactors = FALSE)
            ndata_sender <- sc_data[miR_gene$miR_gene, cell_sender$cell]
            if (nrow(miR_gene) == 1) {
                miR_gene$percent <- as.numeric(.number_cell(ndata_sender))
            } else {
                miR_gene$percent <- as.numeric(apply(ndata_sender, 1, .number_cell))
            }
            miR_gene <- miR_gene[miR_gene$percent >= min_cell_num, ]
            if (nrow(miR_gene) == 0) {
                next
            }
            miR_gene$percent <- log10(miR_gene$percent)
            miR_gene$percent <- miR_gene$percent/max(miR_gene$percent)
            # celltype_receiver
            celltype_receiver <- celltype_pair$celltype_receiver[k]
            cell_receiver <- sc_meta[sc_meta$celltype == celltype_receiver, ]
            # infer miRNA
            miR_gene <- .get_miR_gene_percent(miR_gene, cell_sender, cell_receiver, sc_data)
            miR2tar_sub <- miR2tar[miR2tar$gene %in% miR_gene$miR_gene, ]
            miR_gene <- miR_gene[miR_gene$miR_gene %in% miR2tar_sub$gene,]
            miR_name <- unique(miR2tar_sub$miRNA)
            # ndata_receiver
            ndata_receiver <- sc_data[, cell_receiver$cell]
            ndata_other <- sc_data[unique(miR2tar_sub$target_gene), !colnames(sc_data) %in% colnames(ndata_receiver)]
            ndata_receiver_mean <- .get_ndata_receiver(ndata_receiver)
            ndata_receiver_mean <- ndata_receiver_mean[unique(miR2tar_sub$target_gene), ]
            ndata_receiver <- ndata_receiver[unique(miR2tar_sub$target_gene), ]
            ndata_receiver_mean <- .get_percent_cell(ndata_receiver_mean, ndata_receiver, ndata_other, pvalue)
            if (min(ndata_receiver_mean$target_gene_percent) < min_percent) {
                ndata_receiver_mean[ndata_receiver_mean$target_gene_percent < min_percent, ]$sig <- "NO"
            }
            cci_temp <- data.frame()
            for (i in 1:length(miR_name)) {
                miR_gene2tar <- miR2tar_sub[miR2tar_sub$miRNA == miR_name[i], ]
                miR_gene2tar <- base::merge(miR_gene2tar, ndata_receiver_mean)
                miR_gene2tar <- .get_P(miR_gene2tar, sc_data, cell_sender, cell_receiver, pvalue)
                miR_gene2tar <- miR_gene2tar[miR_gene2tar$sig == "YES", ]
                if (nrow(miR_gene2tar) > 0) {
                    cci_temp <- rbind(cci_temp, miR_gene2tar)
                }
            }
            if (nrow(cci_temp) > 0) {
                cci_temp <- .get_cci(cci_temp, celltype_sender, celltype_receiver, miR_gene)
                cci <- rbind(cci, cci_temp)
            }
            pb$tick()
        }
    }
    object@type <- resolution
    object@cci <- cci
    return(object)
}

#' @title Get miRNA-target interactions
#'
#' @description Get simple results of miRNA-target interactions
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
        cci <- cci[,c(1:4,7:9,16:17)]
        cci$pair <- paste0(cci$celltype_sender, cci$celltype_receiver, cci$miRNA, cci$target_gene)
        cci_pair <- as.data.frame(table(cci$pair), stringsAsFactors = F)
        cci_pair <- cci_pair[cci_pair$Freq > 1,]
        if (nrow(cci_pair) > 0) {
            pair_name <- unique(cci_pair$Var1)
            cci1 <- cci[!cci$pair %in% pair_name,]
            cci2 <- cci[cci$pair %in% pair_name,]
            cci3 <- NULL
            for (i in 1:length(pair_name)) {
                cci_tmp <- cci2[cci2$pair == pair_name[i], ]
                cci_tmp$miRNA_activity <- sum(cci_tmp$miRNA_activity)
                cci_tmp$prob <- max(cci_tmp$prob)
                cci_tmp$score <- max(cci_tmp$score)
                mir_genes <- cci_tmp$miR_gene
                mir_gene <- mir_genes[1]
                for (j in 2:length(mir_genes)) {
                    mir_gene <- paste0(mir_gene, ", ", mir_genes[j])
                }
                cci_tmp$miR_gene <- mir_gene
                cci3 <- rbind(cci3, cci_tmp[1,])
            }
            cci <- rbind(cci1, cci3)
        }
        cci <- cci[ ,-10]
        cci <- unique(cci)
    }
    return(cci)
}

#' @title Get pathways
#'
#' @description Get pathways for target genes
#' @param target_genes Character of one or more target genes
#' @param pathways A data.frame of the system data containing gene-gene interactions and pathways from KEGG and Reactome for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_pathways}}
#' @param species A character meaning species of the target genes.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @return Pathways for one or more target genes
#' @import Matrix
#' @export

get_pathways <- function(target_genes, pathways, species) {
    if (!all(c("src", "dest", "pathway", "species") %in% colnames(pathways))) {
        stop("Please provide a correct pathways data.frame! See demo_pathways()!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    pathways <- pathways[pathways$species == species, ]
    pathways <- pathways[pathways$src %in% target_genes | pathways$dest %in% target_genes, ]
    return(pathways)
}

#' @title Get GO terms
#'
#' @description Get GO terms for target genes
#' @param target_genes Character of one or more target genes
#' @param gene2go A data.frame of the system data containing GO terms for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_gene2go}}
#' @param species A character meaning species of the target genes.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @return GO terms for one or more target genes
#' @import Matrix
#' @export

get_gene2go <- function(target_genes, gene2go, species) {
    if (!all(c("symbol", "GO_term", "species") %in% colnames(gene2go))) {
        stop("Please provide a correct gene2go data.frame! See demo_gene2go()!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    gene2go <- gene2go[gene2go$species == species, ]
    gene2go <- gene2go[gene2go$symbol %in% target_genes, ]
    return(gene2go)
}
