#' @title Demo data of sc_data
#'
#' @description Demo data of sc_data.
#' @details \code{sc_data} can be a data.frame, matrix, or dgCMatrix object, each column representing a cell, each row representing a gene.
#' @return A dgCMatrix object.
#' @import Matrix
#' @export
#' @examples sc_data_demo <- demo_sc_data()

demo_sc_data <- function() {
    cellname <- paste0("cell", 1:6)
    genename <- c("A1BG", "A2M", "A2MP", "NAT1", "NAT2", "NATP")
    countdata <- sample(x = c(rep(0, 50), 1:50), size = 36, replace = TRUE)
    sc_data <- matrix(countdata, nrow = 6, ncol = 6)
    rownames(sc_data) <- genename
    colnames(sc_data) <- cellname
    sc_data <- as(object = sc_data, Class = "dgCMatrix")
    return(sc_data)
}

#' @title Demo data of geneinfo
#'
#' @description Demo data of geneinfo
#' @details \code{geneinfo} must be a \code{data.frame} object with three columns, namely \code{'symbol'}, \code{'synonyms'}, \code{'species'}.
#' @export
#' @examples geneinfo_demo <- demo_geneinfo()

demo_geneinfo <- function() {
    return(head(geneinfo))
}

#' @title Demo data of mir_info
#'
#' @description Demo data of mir_info
#' @details \code{mir_info} must be a \code{data.frame} object with four columns, namely \code{'miRNA'}, \code{'miRNA_mature'}, \code{'gene'}, \code{'species'}
#' @export
#' @examples mir_info_demo <- demo_mir_info()

demo_mir_info <- function() {
    return(head(mir_info))
}

#' @title Demo data of mir2tar
#'
#' @description Demo data of mir2tar
#' @details \code{mir2tar} must be a \code{data.frame} object with four columns, namely \code{'miRNA'}, \code{'miRNA_mature'}, \code{'target_gene'}, \code{'species'}
#' @export
#' @examples mir2tar_demo <- demo_mir2tar()

demo_mir2tar <- function() {
  return(head(mir2tar))
}
