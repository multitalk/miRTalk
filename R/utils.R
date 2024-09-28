.rename_chr <- function(x){
    x <- strsplit(x, split = " ")
    x_new <- NULL
    for (i in 1:length(x)) {
        x1 <- x[[i]]
        if (length(x1) > 1) {
            x2 <- x1[1]
            for (j in 2:length(x1)) {
                x2 <- paste(x2, x1[j], sep = "_")
            }
            x1 <- x2
        }
        x_new <- c(x_new, x1)
    }
    x <- strsplit(x_new, split = "-")
    x_new <- NULL
    for (i in 1:length(x)) {
        x1 <- x[[i]]
        if (length(x1) > 1) {
            x2 <- x1[1]
            for (j in 2:length(x1)) {
                x2 <- paste(x2, x1[j], sep = "_")
            }
            x1 <- x2
        }
        x_new <- c(x_new, x1)
    }
    return(x_new)
}

.show_warning <- function(celltype, celltype_new){
    sc_meta <- data.frame(celltype = celltype, celltype_new = celltype_new, stringsAsFactors = FALSE)
    sc_meta <- unique(sc_meta)
    sc_meta <- sc_meta[sc_meta$celltype != sc_meta$celltype_new, ]
    warning_info <- NULL
    if (nrow(sc_meta) > 0) {
        warning_info <- "celltype of "
        if (nrow(sc_meta) == 1) {
            warning_info <- paste0(warning_info, sc_meta$celltype[1], " has been replaced by ", sc_meta$celltype_new[1])
        } else {
            for (i in 1:nrow(sc_meta)) {
                if (i == nrow(sc_meta)) {
                    warning_info <- paste0(warning_info, "and ",sc_meta$celltype[i])
                } else {
                    warning_info <- paste0(warning_info, sc_meta$celltype[i], ", ")
                }
            }
            warning_info <- paste0(warning_info, " have been replaced by ")
            for (i in 1:nrow(sc_meta)) {
                if (i == nrow(sc_meta)) {
                    warning_info <- paste0(warning_info, "and ",sc_meta$celltype_new[i])
                } else {
                    warning_info <- paste0(warning_info, sc_meta$celltype_new[i], ", ")
                }
            }
        }
    }
    return(warning_info)
}

.get_dgCMatrix <- function(sc_data){
    all_pac <- installed.packages()
    all_pac <- as.data.frame(all_pac, stringsAsFactors = FALSE)
    all_pac <- all_pac[all_pac$Package == "Matrix", ]
    mat_version <- all_pac$Version
    mat_version <- strsplit(mat_version, split = '-')
    mat_version <- mat_version[[1]][1]
    mat_version <- strsplit(mat_version, split = '\\.')
    mat_version <- mat_version[[1]][2]
    mat_version <- substr(mat_version,1,1)
    mat_version <- as.numeric(mat_version)
    if (mat_version >= 5) {
        sc_data <- as(object = sc_data, Class = "CsparseMatrix")
    } else {
        sc_data <- as(object = sc_data, Class = "dgCMatrix")
    }
    return(sc_data)
}

.minmax_normalize <- function(x){
    x_minmax <- (x-min(x))/(max(x)-min(x))
    return(x_minmax)
}

.rank_normalize <- function(x){
    x_rank <- rank(x)
    x_rank <- x_rank/max(x_rank)
    return(x_rank)
}

.get_percent <- function(x) {
    x_name <- unique(names(x))
    x_num <- rep(0, length(x_name))
    names(x_num) <- x_name
    x_gene <- as.data.frame(table(names(x[x > 0])), stringsAsFactors = FALSE)
    x_num[x_gene$Var1] <- x_gene$Freq
    return(list(x_num))
}

.percent_cell <- function(x) {
    return(length(x[x > 0])/length(x))
}

.number_cell <- function(x) {
    return(length(x[x > 0]))
}

.prepare_evbiog <- function(ndata_sender, cell_sender, miR_gene){
    ndata_sender1 <- ndata_sender
    if (nrow(miR_gene) == 1) {
        ndata_sender1 <- as.numeric(cell_sender$evbiog_score)
        names(ndata_sender1) <- paste0("ev", names(ndata_sender1))
        ndata_sender_evbiog <- c(ndata_sender, ndata_sender1)
    } else {
        ndata_sender1[,] <- rep(cell_sender$evbiog_score, each = nrow(ndata_sender))
        colnames(ndata_sender1) <- paste0("ev", colnames(ndata_sender1))
        ndata_sender_evbiog <- cbind(ndata_sender, ndata_sender1)
    }
    return(ndata_sender_evbiog)
}

.get_evmirna_mean <- function(x){
    x_mir <- x[1:(length(x)/2)]
    return(mean(x_mir))
}

.get_evmirna_score <- function(miR_gene, ndata_sender, sc_data_tmp){
    miR_gene$ratio <- 0
    ndata_sender1 <- ndata_sender[miR_gene$miR_gene,]
    sc_data_tmp1 <- sc_data_tmp[miR_gene$miR_gene,]
    if (nrow(miR_gene) == 1) {
        miR_gene$ratio <- sum(ndata_sender1)/sum(sc_data_tmp1)
    }
    else {
        miR_gene$ratio <- rowSums(ndata_sender1)/rowSums(sc_data_tmp1)
    }
    miR_gene$EVmiR_score <- miR_gene$EVmiR_score * miR_gene$ratio
    miR_gene <- miR_gene[,-4]
    return(miR_gene)
}

.get_sc_data_evbiog <- function(sc_data_tmp, sc_meta_tmp, miR2tar_sub){
    sc_data_tmp1 <- sc_data_tmp[miR2tar_sub$gene, ]
    sc_data_tmp2 <- sc_data_tmp1
    if (nrow(miR2tar_sub) == 1) {
        sc_data_tmp2 <- as.numeric(sc_meta_tmp$evbiog_score)
        names(sc_data_tmp2) <- paste0("ev", names(sc_data_tmp2))
        sc_data_weighted <- sc_data_tmp1 * sc_data_tmp2
        sc_data_weighted <- t(as.matrix(sc_data_weighted))
        rownames(sc_data_weighted) <- miR2tar_sub$gene
    } else {
        sc_data_tmp2[,] <- rep(sc_meta_tmp$evbiog_score, each = nrow(sc_data_tmp2))
        colnames(sc_data_tmp2) <- paste0("ev", colnames(sc_data_tmp2))
        sc_data_tmp12 <- cbind(sc_data_tmp1, sc_data_tmp2)
        sc_data_weighted <- apply(sc_data_tmp12, 1, function(x){
            x1 <- x[1:(length(x)/2)]
            x2 <-x[(length(x)/2+1):length(x)]
            x12 <- x1 * x2
            return(x12)
        })
        sc_data_weighted <- t(sc_data_weighted)
    }
    return(sc_data_weighted)
}

.get_sc_data_risc <- function(sc_data_tmp, sc_meta_tmp, miR2tar_sub, regulation){
    sc_data_tmp1 <- sc_data_tmp[miR2tar_sub$target_gene, ]
    sc_data_tmp2 <- sc_data_tmp1
    if (nrow(miR2tar_sub) == 1) {
        sc_data_tmp2 <- as.numeric(sc_meta_tmp$risc_score)
        names(sc_data_tmp2) <- paste0("ev", names(sc_data_tmp2))
        if (regulation == "negative") {
            x1_rank <- rank(-sc_data_tmp1)
        } else {
            x1_rank <- rank(sc_data_tmp1)
        }
        x1_rank <- x1_rank/max(x1_rank)
        sc_data_weighted <- x1_rank * sc_data_tmp2
        sc_data_weighted <- t(as.matrix(sc_data_weighted))
        rownames(sc_data_weighted) <- miR2tar_sub$target_gene
    } else {
        sc_data_tmp2[,] <- rep(sc_meta_tmp$risc_score, each = nrow(sc_data_tmp2))
        colnames(sc_data_tmp2) <- paste0("ev", colnames(sc_data_tmp2))
        sc_data_tmp12 <- cbind(sc_data_tmp1, sc_data_tmp2)
        sc_data_weighted <- apply(sc_data_tmp12, 1, function(x){
            x1 <- x[1:(length(x)/2)]
            if (regulation == "negative") {
                x1_rank <- rank(-x1)
            } else {
                x1_rank <- rank(x1)
            }
            x1_rank <- x1_rank/max(x1_rank)
            x2 <-x[(length(x)/2+1):length(x)]
            x12 <- x1_rank * x2
            return(x12)
        })
        sc_data_weighted <- t(sc_data_weighted)
    }
    return(sc_data_weighted)
}

.get_markers <- function(sc_data, celltype, celltype_other, pvalue, log2fc) {
    marker_genes <- NULL
    for (i in 1:length(celltype_other)) {
        marker_genes_tmp <- Seurat::FindMarkers(sc_data, ident.1 = celltype, ident.2 = celltype_other[i],
            logfc.threshold = log2fc, verbose = FALSE)
        marker_genes_tmp <- marker_genes_tmp[marker_genes_tmp$p_val < pvalue,]
        if (nrow(marker_genes_tmp) > 0) {
            marker_genes <- c(marker_genes, rownames(marker_genes_tmp))
        }
    }
    if (!is.null(marker_genes)) {
        marker_genes <- as.data.frame(table(marker_genes), stringsAsFactors = FALSE)
        marker_genes <- marker_genes[marker_genes$Freq == length(celltype_other),]
        if (nrow(marker_genes) == 0) {
            marker_genes <- NULL
        } else {
            marker_genes <- marker_genes$marker_genes
        }
    }
    return(marker_genes)
}

.get_miR_gene_percent <- function(miR_gene, cell_sender, cell_receiver, sc_data_tmp, pvalue, if_filter_miRNA) {
    .number_cell <- function(x) {
        return(length(x[x > 0]))
    }
    .wilcox_pvalue <- function(x){
        x_sender <- x[names(x) == "sender"]
        x_other <- x[names(x) == "other"]
        x_pvalue <- wilcox.test(x_sender, x_other, alternative = "greater")
        return(as.numeric(x_pvalue$p.value))
    }
    miR_gene$pvalue <- 1
    ndata_temp <- sc_data_tmp[miR_gene$miR_gene, ]
    if (nrow(miR_gene) == 1) {
        miR_gene$percent_sender <- .number_cell(ndata_temp[cell_sender$cell])/.number_cell(ndata_temp)
        miR_gene$percent_receiver <- .number_cell(ndata_temp[cell_receiver$cell])/.number_cell(ndata_temp)
        if (if_filter_miRNA) {
            ndata_sender <- ndata_temp[cell_sender$cell]
            ndata_other <- ndata_temp[!names(ndata_temp) %in% cell_sender$cell]
            miR_gene$pvalue <- as.numeric(wilcox.test(ndata_sender, ndata_other, alternative = "greater")$p.value)
            miR_gene <- miR_gene[miR_gene$pvalue < pvalue, ]
        }
    } else {
        miR_gene$percent_sender <- as.numeric(apply(ndata_temp[, cell_sender$cell], 1, .number_cell))/as.numeric(apply(ndata_temp, 1, .number_cell))
        miR_gene$percent_receiver <- as.numeric(apply(ndata_temp[, cell_receiver$cell], 1, .number_cell))/as.numeric(apply(ndata_temp, 1, .number_cell))
        if (if_filter_miRNA) {
            ndata_sender <- ndata_temp[ ,cell_sender$cell]
            ndata_other <- ndata_temp[ ,!colnames(ndata_temp) %in% cell_sender$cell]
            colnames(ndata_sender) <- rep("sender", ncol(ndata_sender))
            colnames(ndata_other) <- rep("other", ncol(ndata_other))
            miR_gene$pvalue <- base::apply(X = cbind(ndata_sender, ndata_other), MARGIN = 1, FUN = .wilcox_pvalue)
            miR_gene <- miR_gene[miR_gene$pvalue < pvalue, ]
        }
    }
    miR_gene <- miR_gene[ ,-4]
    return(miR_gene)
}

.get_percent_cell <- function(ndata_receiver_mean, ndata_receiver, ndata_other, pvalue, regulation) {
    .percent_cell <- function(x) {
        return(length(x[x > 0])/length(x))
    }
    .wilcox_pvalue_negative <- function(x){
        x_receiver <- x[names(x) == "receiver"]
        x_other <- x[names(x) == "other"]
        x_pvalue <- wilcox.test(x_receiver, x_other, alternative = "less")
        return(as.numeric(x_pvalue$p.value))
    }
    .wilcox_pvalue_positive <- function(x){
        x_receiver <- x[names(x) == "receiver"]
        x_other <- x[names(x) == "other"]
        x_pvalue <- wilcox.test(x_receiver, x_other, alternative = "greater")
        return(as.numeric(x_pvalue$p.value))
    }
    if (nrow(ndata_receiver_mean) == 1) {
        ndata_receiver_mean$target_gene_mean_exp <- as.numeric(mean(ndata_receiver))
        ndata_receiver_mean$target_gene_mean_exp_other <- as.numeric(mean(ndata_other))
        ndata_receiver_mean$target_gene_percent <- .percent_cell(ndata_receiver)
        ndata_receiver_mean$target_gene_percent_other <- .percent_cell(ndata_other)
        if (regulation == "negative") {
            ndata_receiver_mean$pvalue <- as.numeric(wilcox.test(ndata_receiver, ndata_other, alternative = "less")$p.value)
        } else {
            ndata_receiver_mean$pvalue <- as.numeric(wilcox.test(ndata_receiver, ndata_other, alternative = "greater")$p.value)
        }
    } else {
        ndata_receiver_mean$target_gene_mean_exp <- as.numeric(rowMeans(ndata_receiver))
        ndata_receiver_mean$target_gene_mean_exp_other <- as.numeric(rowMeans(ndata_other))
        ndata_receiver_mean$target_gene_percent <- base::apply(X = ndata_receiver, MARGIN = 1, FUN = .percent_cell)
        ndata_receiver_mean$target_gene_percent_other <- base::apply(X = ndata_other, MARGIN = 1, FUN = .percent_cell)
        colnames(ndata_receiver) <- rep("receiver", ncol(ndata_receiver))
        colnames(ndata_other) <- rep("other", ncol(ndata_other))
        if (regulation == "negative") {
            ndata_receiver_mean$pvalue <- base::apply(X = cbind(ndata_receiver, ndata_other), MARGIN = 1, FUN = .wilcox_pvalue_negative)
        } else {
            ndata_receiver_mean$pvalue <- base::apply(X = cbind(ndata_receiver, ndata_other), MARGIN = 1, FUN = .wilcox_pvalue_positive)
        }
    }
    ndata_receiver_mean$sig <- "YES"
    if (nrow(ndata_receiver_mean[is.na(ndata_receiver_mean$pvalue), ]) > 0) {
        ndata_receiver_mean[is.na(ndata_receiver_mean$pvalue), ]$pvalue <- 1
    }
    if (max(ndata_receiver_mean$pvalue) >= pvalue) {
        ndata_receiver_mean[ndata_receiver_mean$pvalue >= pvalue, ]$sig <- "NO"
    }
    exp_fc <- ndata_receiver_mean$target_gene_mean_exp - ndata_receiver_mean$target_gene_mean_exp_other
    if (regulation == "negative") {
        if (max(exp_fc) > 0) {
            ndata_receiver_mean[which(exp_fc > 0), ]$sig <- "NO"
        }
    } else {
        if (min(exp_fc) < 0) {
            ndata_receiver_mean[which(exp_fc < 0), ]$sig <- "NO"
        }
    }
    return(ndata_receiver_mean)
}

.use_human_data <- function(mir_info, mir2tar, geneinfo_species, gene2gene, species, if_use_human_data, if_combine){
    if (if_use_human_data) {
        if (species != "Human") {
            # mir_info
            mir_info_human <- mir_info[mir_info$species == "Human",]
            mir_info_human$gene_other <- stringr::str_to_title(mir_info_human$gene)
            mir_info_human[!mir_info_human$gene_other %in% geneinfo_species$symbol,]$gene_other <- "NA"
            miRNA_mature_other <- mir_info_human[mir_info_human$gene_other != "NA",]$miRNA_mature
            if (species == "Mouse") {
                miRNA_mature_other <- stringr::str_replace(string = miRNA_mature_other, pattern = "hsa",replacement = "mmu")
            } else {
                miRNA_mature_other <- stringr::str_replace(string = miRNA_mature_other, pattern = "hsa",replacement = "rno")
            }
            mir_info_human$miRNA_mature_other <- "NA"
            mir_info_human$miRNA_other <- "NA"
            mir_info_human[mir_info_human$gene_other != "NA",]$miRNA_mature_other <- miRNA_mature_other
            miRNA_mature_other <- stringr::str_remove(string = miRNA_mature_other, pattern = "-3p")
            miRNA_mature_other <- stringr::str_remove(string = miRNA_mature_other, pattern = "-5p")
            mir_info_human[mir_info_human$gene_other != "NA",]$miRNA_other <- miRNA_mature_other
            mir_info_human <- mir_info_human[which(mir_info_human$gene_other != "NA"),]
            # mir2tar
            mir2tar_human <- mir2tar[mir2tar$species == "Human",]
            mir2tar_human <- mir2tar_human[mir2tar_human$miRNA_mature %in% mir_info_human$miRNA_mature,]
            miRNA_mature_other <- mir2tar_human$miRNA_mature
            if (species == "Mouse") {
                miRNA_mature_other <- stringr::str_replace(string = miRNA_mature_other, pattern = "hsa",replacement = "mmu")
            } else {
                miRNA_mature_other <- stringr::str_replace(string = miRNA_mature_other, pattern = "hsa",replacement = "rno")
            }
            mir2tar_human$miRNA_mature <- miRNA_mature_other
            miRNA_mature_other <- stringr::str_remove(string = miRNA_mature_other, pattern = "-3p")
            miRNA_mature_other <- stringr::str_remove(string = miRNA_mature_other, pattern = "-5p")
            mir2tar_human$miRNA <- miRNA_mature_other
            # map target gene
            gene2gene <- gene2gene[gene2gene$species == species, ]
            mir2tar_human <- mir2tar_human[mir2tar_human$target_gene %in% gene2gene$gene,]
            gene2gene <- gene2gene[gene2gene$gene %in% mir2tar_human$target_gene, -3]
            colnames(gene2gene)[1] <- "target_gene"
            mir2tar_human <- merge(mir2tar_human, gene2gene)
            mir2tar_human$target_gene <- mir2tar_human$gene_other
            mir2tar_human <- mir2tar_human[ ,c(2:4,1,5:9)]
            mir_info_human <- mir_info_human[, c("miRNA_other","miRNA_mature_other","gene_other","species","EV_evidence","tissue_TarBase","celltype_TarBase")]
            colnames(mir_info_human)[1:3] <- c("miRNA", "miRNA_mature", "gene")
            mir_info_human$species <- paste0("Human-->", species)
            if(if_combine){
                mir_info_raw <- mir_info[mir_info$species == species, ]
                mir2tar_raw <- mir2tar[mir2tar$species == species, ]
                mir_info <- rbind(mir_info_raw, mir_info_human)
                mir2tar <- rbind(mir2tar_raw, mir2tar_human)
            } else {
                mir_info <- mir_info_human
                mir2tar <- mir2tar_human
            }
        } else {
            mir_info <- mir_info[mir_info$species == species, ]
            mir2tar <- mir2tar[mir2tar$species == species, ]
        }
    } else {
        mir_info <- mir_info[mir_info$species == species, ]
        mir2tar <- mir2tar[mir2tar$species == species, ]
    }
    return(list(mir_info, mir2tar))
}

.get_ndata_receiver <- function(ndata_receiver) {
    ndata_receiver_mean <- data.frame(gene = rownames(ndata_receiver), mean_exp = as.numeric(rowMeans(ndata_receiver)), stringsAsFactors = FALSE)
    ndata_receiver_mean$rank <- rank(-ndata_receiver_mean$mean_exp)
    ndata_receiver_mean$activity <- ndata_receiver_mean$mean_exp/max(ndata_receiver_mean$mean_exp)
    rownames(ndata_receiver_mean) <- ndata_receiver_mean$gene
    colnames(ndata_receiver_mean)[1] <- "target_gene"
    ndata_receiver_mean <- ndata_receiver_mean[, -2]
    return(ndata_receiver_mean)
}

.get_miR2tar <- function(mir_info, miR2tar, resolution) {
    if (nrow(mir_info) == 0) {
        stop("No shared miRNAs found in sc_data!")
    }
    if (nrow(miR2tar) == 0) {
        stop("No shared miRNAs and targets found in sc_data!")
    }
    if (resolution == "mature") {
        mir_info <- unique(mir_info[, c("miRNA_mature", "gene")])
        miR2tar <- unique(miR2tar[, c("miRNA_mature", "target_gene")])
        miR2tar <- base::merge(mir_info, miR2tar)
        colnames(miR2tar)[1] <- "miRNA"
    } else {
        mir_info <- unique(mir_info[, c("miRNA", "gene")])
        miR2tar <- unique(miR2tar[, c("miRNA", "target_gene")])
        miR2tar <- base::merge(mir_info, miR2tar)
    }
    return(miR2tar)
}

.get_per_test_list <- function(sc_meta_tmp, per_num){
    cellnames <- sc_meta_tmp$cell
    per_test_list <- list()
    for(i in 1:per_num){
        set.seed(i)
        cell_index <- sample(x = 1:length(cellnames), size = length(cellnames), replace = F)
        per_test_list[[i]] <- cellnames[cell_index]
    }
    return(per_test_list)
}

.get_S <- function(miR_gene2tar, sc_data_evbiog, sc_data_risc, cell_sender, cell_receiver, per_test_list){
    get_sscore <- function(x) {
        x_miR_sender <- x[names(x) == "miR_sender"]
        x_tar_receiver <- x[names(x) == "tar_receiver"]
        x_miR_other_mean <- x[names(x) == "miR_other_mean"]
        x_tar_other_mean <- x[names(x) == "tar_other_mean"]
        x_miR <- x_miR_sender[x_miR_sender > x_miR_other_mean]
        score1 <- sum(x_miR)/length(x_miR_sender)
        x_tar <- x_tar_receiver[x_tar_receiver > x_tar_other_mean]
        score2 <- sum(x_tar)/length(x_tar_receiver)
        sscore <- score1 * score2
        return(sscore)
    }
    if (nrow(miR_gene2tar) == 1) {
        # sender
        miR_other_mean <- mean(sc_data_evbiog[miR_gene2tar$gene, colnames(sc_data_evbiog[,!colnames(sc_data_evbiog) %in% cell_sender$cell])])
        miR_sender1 <- sc_data_evbiog[miR_gene2tar$gene, cell_sender$cell]
        names(miR_sender1) <- rep("miR_sender", length(miR_sender1))
        miR_sender2 <- miR_sender1[c(1,2)]
        miR_sender2[1] <- as.numeric(miR_other_mean)
        names(miR_sender2) <- c("miR_other_mean", "tmp_sender")
        # receiver
        tar_other_mean <- mean(sc_data_risc[miR_gene2tar$target_gene, colnames(sc_data_risc[,!colnames(sc_data_risc) %in% cell_receiver$cell])])
        tar_receiver1 <- sc_data_risc[miR_gene2tar$target_gene, cell_receiver$cell]
        names(tar_receiver1) <- rep("tar_receiver", length(tar_receiver1))
        tar_receiver2 <- tar_receiver1[c(1,2)]
        tar_receiver2[1] <- as.numeric(tar_other_mean)
        names(tar_receiver2) <- c("tar_other_mean", "tmp_receiver")
        miRtar_sender_receiver <- c(miR_sender1, miR_sender2, tar_receiver1, tar_receiver2)
        miR_gene2tar$score <- get_sscore(miRtar_sender_receiver)
        per_test_res <- list()
        for (i in 1:length(per_test_list)) {
            colnames(sc_data_evbiog) <- per_test_list[[i]]
            colnames(sc_data_risc) <- per_test_list[[i]]
            # sender
            miR_other_mean <- mean(sc_data_evbiog[miR_gene2tar$gene, colnames(sc_data_evbiog[,!colnames(sc_data_evbiog) %in% cell_sender$cell])])
            miR_sender1 <- sc_data_evbiog[miR_gene2tar$gene, cell_sender$cell]
            names(miR_sender1) <- rep("miR_sender", length(miR_sender1))
            miR_sender2 <- miR_sender1[c(1,2)]
            miR_sender2[1] <- as.numeric(miR_other_mean)
            names(miR_sender2) <- c("miR_other_mean", "tmp_sender")
            # receiver
            tar_other_mean <- mean(sc_data_risc[miR_gene2tar$target_gene, colnames(sc_data_risc[,!colnames(sc_data_risc) %in% cell_receiver$cell])])
            tar_receiver1 <- sc_data_risc[miR_gene2tar$target_gene, cell_receiver$cell]
            names(tar_receiver1) <- rep("tar_receiver", length(tar_receiver1))
            tar_receiver2 <- tar_receiver1[c(1,2)]
            tar_receiver2[1] <- as.numeric(tar_other_mean)
            names(tar_receiver2) <- c("tar_other_mean", "tmp_receiver")
            miRtar_sender_receiver <- c(miR_sender1, miR_sender2, tar_receiver1, tar_receiver2)
            per_test_res[[i]] <- get_sscore(miRtar_sender_receiver)
        }
    } else {
        # sender
        miR_other_mean <- rowMeans(sc_data_evbiog[miR_gene2tar$gene, colnames(sc_data_evbiog[,!colnames(sc_data_evbiog) %in% cell_sender$cell])])
        miR_sender1 <- sc_data_evbiog[miR_gene2tar$gene, cell_sender$cell]
        colnames(miR_sender1) <- rep("miR_sender", ncol(miR_sender1))
        miR_sender2 <- miR_sender1[ ,c(1,2)]
        miR_sender2[ ,1] <- as.numeric(miR_other_mean)
        colnames(miR_sender2) <- c("miR_other_mean", "tmp_sender")
        # receiver
        tar_other_mean <- rowMeans(sc_data_risc[miR_gene2tar$target_gene, colnames(sc_data_risc[,!colnames(sc_data_risc) %in% cell_receiver$cell])])
        tar_receiver1 <- sc_data_risc[miR_gene2tar$target_gene, cell_receiver$cell]
        colnames(tar_receiver1) <- rep("tar_receiver", ncol(tar_receiver1))
        tar_receiver2 <- tar_receiver1[ ,c(1,2)]
        tar_receiver2[ ,1] <- as.numeric(tar_other_mean)
        colnames(tar_receiver2) <- c("tar_other_mean", "tmp_receiver")
        miRtar_sender_receiver <- cbind(miR_sender1, miR_sender2, tar_receiver1, tar_receiver2)
        miR_gene2tar$score <- as.numeric(base::apply(miRtar_sender_receiver, 1, FUN = get_sscore))
        per_test_res <- list()
        for (i in 1:length(per_test_list)) {
            colnames(sc_data_evbiog) <- per_test_list[[i]]
            colnames(sc_data_risc) <- per_test_list[[i]]
            # sender
            miR_other_mean <- rowMeans(sc_data_evbiog[miR_gene2tar$gene, colnames(sc_data_evbiog[,!colnames(sc_data_evbiog) %in% cell_sender$cell])])
            miR_sender1 <- sc_data_evbiog[miR_gene2tar$gene, cell_sender$cell]
            colnames(miR_sender1) <- rep("miR_sender", ncol(miR_sender1))
            miR_sender2 <- miR_sender1[ ,c(1,2)]
            miR_sender2[ ,1] <- as.numeric(miR_other_mean)
            colnames(miR_sender2) <- c("miR_other_mean", "tmp_sender")
            # receiver
            tar_other_mean <- rowMeans(sc_data_risc[miR_gene2tar$target_gene, colnames(sc_data_risc[,!colnames(sc_data_risc) %in% cell_receiver$cell])])
            tar_receiver1 <- sc_data_risc[miR_gene2tar$target_gene, cell_receiver$cell]
            colnames(tar_receiver1) <- rep("tar_receiver", ncol(tar_receiver1))
            tar_receiver2 <- tar_receiver1[ ,c(1,2)]
            tar_receiver2[ ,1] <- as.numeric(tar_other_mean)
            colnames(tar_receiver2) <- c("tar_other_mean", "tmp_receiver")
            miRtar_sender_receiver <- cbind(miR_sender1, miR_sender2, tar_receiver1, tar_receiver2)
            per_test_res[[i]] <- as.numeric(base::apply(miRtar_sender_receiver, 1, FUN = get_sscore))
        }
    }
    per_test_res <- as.data.frame(per_test_res)
    colnames(per_test_res) <- paste0("V",1:ncol(per_test_res))
    per_test_res$real <- miR_gene2tar$score
    miR_gene2tar$prob <- apply(per_test_res, 1, function(x){
        x_per <- x[names(x) != "real"]
        x_real <- x[names(x) == "real"]
        return(length(x_per[x_per > x_real])/length(x_per))
    })
    return(miR_gene2tar)
}

.get_cci <- function(cci_temp, celltype_sender, celltype_receiver, miR_gene) {
    cci_temp$celltype_sender <- celltype_sender
    cci_temp$celltype_receiver <- celltype_receiver
    miR_gene <- miR_gene[, -2]
    colnames(cci_temp)[3] <- "miR_gene"
    cci_temp <- base::merge(cci_temp, miR_gene)
    cci_temp <- cci_temp[, c(14:15, 3, 1, 17:18, 16, 2, 4:13)]
    colnames(cci_temp)[9] <- "target_gene_rank"
    colnames(cci_temp)[10] <- "target_gene_activity"
    cci_temp <- cci_temp[order(-cci_temp$prob), ]
    return(cci_temp)
}

.get_specifity <- function(cci_tmp){
    cci_tmp$miR2tar <- paste0(cci_tmp$miRNA, ":", cci_tmp$target_gene)
    celltype_receiver <- unique(cci_tmp$celltype_receiver)
    res_specifity <- data.frame()
    cci_tmp$pair <- paste0(cci_tmp$celltype_sender, ":", cci_tmp$miRNA)
    for (i in 1:length(celltype_receiver)) {
        cci_tmp1 <- cci_tmp[cci_tmp$celltype_receiver == celltype_receiver[i],]
        target_genes <- unique(cci_tmp1$target_gene)
        for (j in 1:length(target_genes)) {
            cci_tmp1_target <- cci_tmp1[cci_tmp1$target_gene == target_genes[j],]
            raw_nrow <- nrow(cci_tmp1_target)
            cci_pair <- unique(cci_tmp1_target$pair)
            for (k in 1:length(cci_pair)) {
                cci_tmp1_target1 <- cci_tmp1_target[cci_tmp1_target$pair == cci_pair[k],]
                if (nrow(cci_tmp1_target1) > 1) {
                    cci_tmp1_target1 <- cci_tmp1_target1[which.max(cci_tmp1_target1$EVmiR_score), ]
                }
                cci_tmp1_target <- rbind(cci_tmp1_target, cci_tmp1_target1)
            }
            cci_tmp1_target <- cci_tmp1_target[-c(1:raw_nrow), ]
            cci_tmp1_target$specifity <- cci_tmp1_target$EVmiR_score/sum(cci_tmp1_target$EVmiR_score)
            res_specifity <- rbind(res_specifity, cci_tmp1_target)
        }
    }
    res_specifity <- res_specifity[,-10]
    return(res_specifity)
}

.get_coord <- function(cci, cellpair, y_len, miR_name, x_len) {
    y <- 1:y_len
    x <- 1:x_len
    cellpair <- rep(cellpair, x_len)
    y <- rep(y, x_len)
    miR_name <- rep(miR_name, each = y_len)
    x <- rep(x, each = y_len)
    cci_miR_temp <- data.frame(cellpair = cellpair, miRNA = miR_name, x = x, y = y, stringsAsFactors = FALSE)
    cci_miR_temp$EVmiR_score <- 0
    cci_miR_temp$score <- 0
    for (i in 1:nrow(cci)) {
        cci_miR_temp[cci_miR_temp$cellpair == cci$cellpair[i] & cci_miR_temp$miRNA == cci$miRNA[i], ]$EVmiR_score <- cci$EVmiR_score[i]
        cci_miR_temp[cci_miR_temp$cellpair == cci$cellpair[i] & cci_miR_temp$miRNA == cci$miRNA[i], ]$score <- cci$score[i]
    }
    return(cci_miR_temp)
}

.get_bubble <- function(cci) {
    cci$pair <- paste0(cci$cellpair, cci$miRNA)
    cci_tmp <- unique(cci$pair)
    cci_res <- data.frame()
    for (i in 1:length(cci_tmp)) {
        cci1 <- cci[cci$pair == cci_tmp[i], ]
        cci1 <- cci1[cci1$score == max(cci1$score), ]
        cci_res <- rbind(cci_res, cci1[1, -5])
    }
    return(cci_res)
}

.get_miR2tar_circle <- function(cci) {
    cci$pair <- paste0(cci$miRNA, cci$target_gene)
    cci_tmp <- unique(cci$pair)
    cci_res <- data.frame()
    for (i in 1:length(cci_tmp)) {
        cci1 <- cci[cci$pair == cci_tmp[i], ]
        cci1 <- cci1[cci1$score == max(cci1$score), ]
        cci_res <- rbind(cci_res, cci1[1,])
    }
    cci_res <- cci_res[,-6]
    return(cci_res)
}

#' @title Show miRTalk object
#'
#' @param object miRTalk object
#' @return miRTalk object
#' @import Matrix
#' @importFrom methods show
#'
#' @export

setMethod(
    f = 'show',
    signature = 'miRTalk',
    definition = function(object) {
        cat("An object of class miRTalk", "\n")
        cci <- object@cci
        cat(paste0(nrow(cci), " EV-derived miRNA-target interactions"), "\n")
        return(invisible(x = NULL))
    }
)
