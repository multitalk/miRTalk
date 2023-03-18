#' @title Pre-processing step: revising gene symbols
#'
#' @description Revise genes according to NCBI Gene symbols updated in June 19, 2022 for count matrix, user-custom mir_info data.frame, and mir2tar data.frame
#' @param data A data.frame or matrix or dgCMatrix
#' @param data_type A character to define the type of \code{data}, select \code{'count'} for the data matrix, \code{'mir_info'} for the mir_info data.frame, \code{'mir2tar'} for the mir2tar data.frame, \code{'pathways'} for the pathways data.frame, \code{'GO_BP'} for the GO_BP data.frame
#' @param species Species of the data.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @param geneinfo A data.frame of the system data containing gene symbols of \code{'Human'}, \code{'Mouse'} and \code{'Rat'} updated on June 19, 2022 for revising gene symbols
#' @return A new data.frame, matrix, or dgCMatrix.
#' @import Matrix
#' @importFrom crayon cyan
#' @export

rev_gene <- function(data = NULL, data_type = NULL, species = NULL, geneinfo = NULL) {
    if (is.null(data)) {
        stop("Please provide the data for revsing gene symbols!")
    }
    if (is.null(data_type) | !is.character(data_type)) {
        stop("Please provide a correct data_type, i.e., 'count', 'mir_info', 'mir2tar', 'pathways', or 'GO_BP'!")
    }
    if (is.null(geneinfo)) {
        stop("Please provide geneinfo for revsing gene symbols, or use the system data like 'geneinfo = geneinfo'")
    }
    if (length(data_type) > 1 | !data_type %in% c("count", "mir_info", "mir2tar", "pathways", "GO_BP")) {
        stop("Please provide a correct data_type, i.e., 'count', 'mir_info', 'mir2tar', 'pathways', or 'GO_BP'!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    # define the species
    if (species == "Human") {
        geneinfo <- geneinfo[geneinfo$species == "Human", ]
    }
    if (species == "Mouse") {
        geneinfo <- geneinfo[geneinfo$species == "Mouse", ]
    }
    if (species == "Rat") {
        geneinfo <- geneinfo[geneinfo$species == "Rat", ]
    }
    # revise matrix
    if (data_type == "count") {
        if (is(data, "data.frame")) {
            data <- methods::as(as.matrix(data), "dgCMatrix")
        }
        if (is(data, "matrix")) {
            data <- methods::as(data, "dgCMatrix")
        }
        if (!is(data, "dgCMatrix")) {
            stop("st_data must be a data.frame or matrix or dgCMatrix!")
        }
        Sys.sleep(1)
        # revise gene symbols
        genename <- rownames(data)
        genename1 <- genename[genename %in% geneinfo$symbol]
        if (length(genename1) == 0) {
            stop("Please ensure the rownames of data are gene symbols! See demo_st_data()!")
        }
        genename2 <- genename[!genename %in% geneinfo$symbol]
        if (length(genename2) > 0) {
            genename3 <- genename2[genename2 %in% geneinfo$synonyms]
            if (length(genename3) > 0) {
                genename4 <- rep("NA", length(genename3))
                for (i in 1:length(genename3)) {
                    d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                    if (length(d1) == 1) {
                        genename4[i] <- d1
                    }
                }
                genename3 <- c(genename1, genename3)
                genename4 <- c(genename1, genename4)
                genedata <- data.frame(raw_name = genename3, new_name = genename4, stringsAsFactors = F)
                genedata <- genedata[!genedata$new_name == "NA", ]
                genedata1 <- as.data.frame(table(genedata$new_name), stringsAsFactors = F)
                genedata1 <- genedata1[genedata1$Freq == 1, ]
                genedata <- genedata[genedata$new_name %in% genedata1$Var1, ]
                data <- data[genedata$raw_name, ]
                rownames(data) <- genedata$new_name
            }
        } else {
            data <- data[rownames(data) %in% geneinfo$symbol, ]
        }
    }
    # revise mir_info
    if (data_type == "mir_info") {
        if (!is.data.frame(data)) {
            stop("data must be a data.frame when data_type is 'mir_info'!")
        }
        cat(crayon::cyan("Revising gene symbols for mir_info data.frame", "\n"))
        Sys.sleep(1)
        if (all(c("miRNA", "miRNA_mature", "gene", "species") %in% colnames(data))) {
            # gene
            genename <- unique(data$gene)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the gene of data are gene symbols! See demo_mir_info()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                    for (i in 1:length(genename3)) {
                        d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                        if (length(d1) == 1) {
                            data[data$gene == genename3[i], ]$gene <- d1
                        }
                    }
                }
            }
            data <- data[data$gene %in% geneinfo$symbol, ]
        } else {
            stop("Please provide a correct mir_info data.frame! See demo_mir_info()!")
        }
    }
    # revise mir2tar
    if (data_type == "mir2tar") {
        if (!is.data.frame(data)) {
            stop("data must be a data.frame when data_type is 'mir2tar'!")
        }
        cat(crayon::cyan("Revising gene symbols for mir2tar data.frame", "\n"))
        Sys.sleep(1)
        if (all(c("miRNA", "miRNA_mature", "target_gene", "species") %in% colnames(data))) {
            # gene
            genename <- unique(data$target_gene)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
              stop("Please ensure the target_gene of data are gene symbols! See demo_mir_info()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                    for (i in 1:length(genename3)) {
                        d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                        if (length(d1) == 1) {
                            data[data$target_gene == genename3[i], ]$target_gene <- d1
                        }
                    }
                }
            }
            data <- data[data$target_gene %in% geneinfo$symbol, ]
        } else {
            stop("Please provide a correct mir2tar data.frame! See demo_mir2tar()!")
        }
    }
    # revise pathways
    if (data_type == "pathways") {
        if (!is.data.frame(data)) {
            stop("data must be a data.frame when data_type is 'pathways'!")
        }
        cat(crayon::cyan("Revising gene symbols for pathways data.frame", "\n"))
        Sys.sleep(1)
        if (all(c("src", "dest", "pathway", "species") %in%
            colnames(data))) {
            # src
            genename <- unique(data$src)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the src of data are gene symbols! See demo_lrpairs()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                    for (i in 1:length(genename3)) {
                        d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                        if (length(d1) == 1) {
                            data[data$src == genename3[i], ]$src <- d1
                        }
                    }
                }
            }
            data <- data[data$src %in% geneinfo$symbol, ]
            # dest
            genename <- unique(data$dest)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the dest of data are gene symbols! See demo_pathways()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                    for (i in 1:length(genename3)) {
                        d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                        if (length(d1) == 1) {
                            data[data$dest == genename3[i], ]$dest <- d1
                        }
                    }
                }
            }
            data <- data[data$dest %in% geneinfo$symbol, ]
        } else {
            stop("Please provide a correct pathways data.frame! See demo_pathways()!")
        }
    }
    return(data)
}
