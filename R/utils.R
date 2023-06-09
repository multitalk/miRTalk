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

.get_miR_gene_percent <- function(miR_gene, cell_sender, cell_receiver, sc_data) {
    .number_cell <- function(x) {
        return(length(x[x > 0]))
    }
    ndata_temp <- sc_data[miR_gene$miR_gene, ]
    if (nrow(miR_gene) == 1) {
        miR_gene$percent_sender <- .number_cell(ndata_temp[cell_sender$cell])/.number_cell(ndata_temp)
        miR_gene$percent_receiver <- .number_cell(ndata_temp[cell_receiver$cell])/.number_cell(ndata_temp)
    } else {
        miR_gene$percent_sender <- as.numeric(apply(ndata_temp[, cell_sender$cell], 1, .number_cell))/as.numeric(apply(ndata_temp, 1, .number_cell))
        miR_gene$percent_receiver <- as.numeric(apply(ndata_temp[, cell_receiver$cell], 1, .number_cell))/as.numeric(apply(ndata_temp, 1, .number_cell))
    }
    return(miR_gene)
}

.get_percent_cell <- function(ndata_receiver_mean, ndata_receiver, ndata_other, pvalue) {
    .percent_cell <- function(x) {
        return(length(x[x > 0])/length(x))
    }
    .wilcox_pvalue <- function(x){
        x_receiver <- x[names(x) == "receiver"]
        x_other <- x[names(x) == "other"]
        x_pvalue <- wilcox.test(x_receiver, x_other, alternative = "less")
        return(as.numeric(x_pvalue$p.value))
    }
    if (nrow(ndata_receiver_mean) == 1) {
        ndata_receiver_mean$target_gene_mean_exp <- as.numeric(mean(ndata_receiver))
        ndata_receiver_mean$target_gene_mean_exp_other <- as.numeric(mean(ndata_other))
        ndata_receiver_mean$target_gene_percent <- .percent_cell(ndata_receiver)
        ndata_receiver_mean$target_gene_percent_other <- .percent_cell(ndata_other)
        ndata_receiver_mean$pvalue <- as.numeric(wilcox.test(ndata_receiver, ndata_other, alternative = "less")$p.value)
    } else {
        ndata_receiver_mean$target_gene_mean_exp <- as.numeric(rowMeans(ndata_receiver))
        ndata_receiver_mean$target_gene_mean_exp_other <- as.numeric(rowMeans(ndata_other))
        ndata_receiver_mean$target_gene_percent <- base::apply(X = ndata_receiver, MARGIN = 1, FUN = .percent_cell)
        ndata_receiver_mean$target_gene_percent_other <- base::apply(X = ndata_other, MARGIN = 1, FUN = .percent_cell)
        colnames(ndata_receiver) <- rep("receiver", ncol(ndata_receiver))
        colnames(ndata_other) <- rep("other", ncol(ndata_other))
        ndata_receiver_mean$pvalue <- base::apply(X = cbind(ndata_receiver, ndata_other), MARGIN = 1, FUN = .wilcox_pvalue)
    }
    ndata_receiver_mean$sig <- "YES"
    if (nrow(ndata_receiver_mean[is.na(ndata_receiver_mean$pvalue), ]) > 0) {
        ndata_receiver_mean[is.na(ndata_receiver_mean$pvalue), ]$pvalue <- 1
    }
    if (max(ndata_receiver_mean$pvalue) >= pvalue) {
        ndata_receiver_mean[ndata_receiver_mean$pvalue >= pvalue, ]$sig <- "NO"
    }
    exp_fc <- ndata_receiver_mean$target_gene_mean_exp - ndata_receiver_mean$target_gene_mean_exp_other
    if (max(exp_fc) > 0) {
        ndata_receiver_mean[which(exp_fc > 0), ]$sig <- "NO"
    }
    return(ndata_receiver_mean)
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

.get_miR2tar <- function(object, resolution) {
    miR_info <- object@miR
    miR2tar <- object@miR2tar
    if (nrow(miR_info) == 0) {
        stop("No shared miRNAs found in sc_data!")
    }
    if (nrow(miR2tar) == 0) {
        stop("No shared miRNAs and targets found in sc_data!")
    }
    if (resolution == "mature") {
        miR_info <- unique(miR_info[, c("miRNA_mature", "gene")])
        miR2tar <- unique(miR2tar[, c("miRNA_mature", "target_gene")])
        miR2tar <- base::merge(miR_info, miR2tar)
        colnames(miR2tar)[1] <- "miRNA"
    } else {
        miR_info <- unique(miR_info[, c("miRNA", "gene")])
        miR2tar <- unique(miR2tar[, c("miRNA", "target_gene")])
        miR2tar <- base::merge(miR_info, miR2tar)
    }
    return(miR2tar)
}

.get_P <- function(miR_gene2tar, sc_data, cell_sender, cell_receiver, pvalue){
    .get_PI <- function(x) {
        x_sender_miR <- x[names(x) == "miR_sender"]
        x_receiver_miR <- x[names(x) == "miR_reciver"]
        x_receiver <- x[names(x) == "receiver"]
        x_sender_other <- x[names(x) == "other_sender"]
        x_receiver_other <- x[names(x) == "other_receiver"]
        p_s <- length(x_sender_miR[x_sender_miR > x_sender_other])/length(x_sender_miR)
        x_receiver_miR <- x_receiver_miR[which(x_receiver < x_receiver_other)]
        if (length(x_receiver_miR) > 0) {
            x_p <- 0
            for (i in 1:length(x_sender_miR)) {
                x_p1 <- length(x_receiver_miR[x_receiver_miR < x_sender_miR[i]])
                x_p <- x_p + x_p1
            }
            p_i <- x_p / (length(x_sender_miR) * length(x_receiver))
        } else {
            p_i <- 0
        }
        p_i <- sqrt(p_s*p_i)
        return(p_i)
    }
    if (nrow(miR_gene2tar) == 1) {
        # sender
        miR_gene <- mean(sc_data[miR_gene2tar$gene, colnames(sc_data[,!colnames(sc_data) %in% cell_sender$cell])])
        ndata_sender1 <- sc_data[miR_gene2tar$gene, cell_sender$cell]
        names(ndata_sender1) <- rep("miR_sender", length(ndata_sender1))
        ndata_sender2 <- sc_data[miR_gene2tar$gene, cell_receiver$cell]
        names(ndata_sender2) <- rep("miR_receiver", length(ndata_sender2))     
        ndata_sender3 <- ndata_sender2[1:2]
        ndata_sender3[1] <- as.numeric(miR_gene)
        ndata_sender3[2] <- pvalue
        names(ndata_sender3) <- c("other_sender", "pvalue_sender")
        # receiver
        miR_targene <- mean(sc_data[miR_gene2tar$target_gene, colnames(sc_data[,!colnames(sc_data) %in% cell_receiver$cell])])
        ndata_receiver1 <- sc_data[miR_gene2tar$target_gene, cell_receiver$cell]
        names(ndata_receiver1) <- rep("receiver", length(ndata_receiver1))
        ndata_receiver2 <- ndata_receiver1[1:2]
        ndata_receiver2[1] <- as.numeric(miR_targene)
        ndata_receiver2[2] <- pvalue
        names(ndata_receiver2) <- c("other_receiver", "pvalue_receiver")
        ndata_sender_receiver <- c(ndata_sender1, ndata_receiver1, ndata_receiver2)
        miR_gene2tar$prob <- .get_PI(ndata_sender_receiver)
    } else {
        # sender
        miR_gene <- rowMeans(sc_data[miR_gene2tar$gene, colnames(sc_data[,!colnames(sc_data) %in% cell_sender$cell])])
        ndata_sender1 <- sc_data[miR_gene2tar$gene, cell_sender$cell]
        colnames(ndata_sender1) <- rep("miR_sender", ncol(ndata_sender1))
        ndata_sender2 <- sc_data[miR_gene2tar$gene, cell_receiver$cell]
        colnames(ndata_sender2) <- rep("miR_reciver", ncol(ndata_sender2))
        ndata_sender3 <- ndata_sender1[,1:2]
        ndata_sender3[,1] <- as.numeric(miR_gene)
        ndata_sender3[,2] <- pvalue
        colnames(ndata_sender3) <- c("other_sender", "pvalue_sender")
        # receiver
        miR_targene <- rowMeans(sc_data[miR_gene2tar$target_gene, colnames(sc_data[,!colnames(sc_data) %in% cell_receiver$cell])])
        ndata_receiver1 <- sc_data[miR_gene2tar$target_gene, cell_receiver$cell]
        colnames(ndata_receiver1) <- rep("receiver", ncol(ndata_receiver1))
        ndata_receiver2 <- ndata_receiver1[,1:2]
        ndata_receiver2[,1] <- as.numeric(miR_targene)
        ndata_receiver2[,2] <- pvalue
        colnames(ndata_receiver2) <- c("other_receiver", "pvalue_receiver")
        ndata_sender_receiver <- cbind(ndata_sender1, ndata_sender2, ndata_sender3, ndata_receiver1, ndata_receiver2)
        miR_gene2tar$prob <- as.numeric(base::apply(ndata_sender_receiver, 1, FUN = .get_PI))
    }
    return(miR_gene2tar)
}

.get_cci <- function(cci_temp, celltype_sender, celltype_receiver, miR_gene) {
    cci_temp$celltype_sender <- celltype_sender
    cci_temp$celltype_receiver <- celltype_receiver
    miR_gene$miRNA_activity <- sqrt(miR_gene$percent_sender * miR_gene$percent)
    miR_gene <- miR_gene[, -2]
    colnames(cci_temp)[3] <- "miR_gene"
    cci_temp <- base::merge(cci_temp, miR_gene)
    cci_temp <- cci_temp[, c(13:14, 3, 1, 15:17, 2, 5:12)]
    colnames(cci_temp)[9] <- "target_gene_activity"
    cci_temp <- cci_temp[order(-cci_temp$prob), ]
    cci_temp$score <- cci_temp$miRNA_activity * (1 - cci_temp$target_gene_activity)
    return(cci_temp)
}

.get_coord <- function(cci, cellpair, y_len, miR_name, x_len) {
    y <- 1:y_len
    x <- 1:x_len
    cellpair <- rep(cellpair, x_len)
    y <- rep(y, x_len)
    miR_name <- rep(miR_name, each = y_len)
    x <- rep(x, each = y_len)
    cci_miR_temp <- data.frame(cellpair = cellpair, miRNA = miR_name, x = x, y = y, stringsAsFactors = FALSE)
    cci_miR_temp$miRNA_activity <- 0
    cci_miR_temp$prob <- 0
    for (i in 1:nrow(cci)) {
        cci_miR_temp[cci_miR_temp$cellpair == cci$cellpair[i] & cci_miR_temp$miRNA == cci$miRNA[i], ]$miRNA_activity <- cci$miRNA_activity[i]
        cci_miR_temp[cci_miR_temp$cellpair == cci$cellpair[i] & cci_miR_temp$miRNA == cci$miRNA[i], ]$prob <- cci$prob[i]
    }
    return(cci_miR_temp)
}

.get_bubble <- function(cci) {
    cci$pair <- paste0(cci$cellpair, cci$miRNA, cci$miRNA_activity)
    cci_tmp <- unique(cci$pair)
    cci_res <- data.frame()
    for (i in 1:length(cci_tmp)) {
        cci1 <- cci[cci$pair == cci_tmp[i], ]
        cci1 <- cci1[cci1$prob == max(cci1$prob), ]
        cci_res <- rbind(cci_res, cci1[1, -5])
    }
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