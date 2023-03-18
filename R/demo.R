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
    gene1 <- c("A1BG", "A1BG", "A2MP1", "Aco1", "Alb1")
    gene2 <- c("A1B", "ABG", "A2MP", "Aco", "Alb")
    species <- c("Human", "Human", "Human", "Mouse", "Rat")
    geneinfo_demo <- data.frame(symbol = gene1, synonyms = gene2, species = species, stringsAsFactors = FALSE)
    return(geneinfo_demo)
}

#' @title Demo data of gene2go
#'
#' @description Demo data of gene2go
#' @details \code{gene2go} must be a \code{data.frame} object with three columns, namely \code{'symbol'}, \code{'GO_term'}, \code{'species'}.
#' @export
#' @examples gene2go_demo <- demo_gene2go()

demo_gene2go <- function() {
    gene <- c("A1BG", "A1BG", "A1BG", "Zzz3", "Zyx")
    GO_term <- c("molecular_function", "extracellular region", "extracellular space", "DNA binding", "metal ion binding")
    species <- c("Human", "Human", "Human", "Mouse", "Rat")
    gene2go_demo <- data.frame(symbol = gene, GO_term = GO_term, species = species, stringsAsFactors = FALSE)
    return(gene2go_demo)
}

#' @title Demo data of mir_info
#'
#' @description Demo data of mir_info
#' @details \code{mir_info} must be a \code{data.frame} object with four columns, namely \code{'miRNA'}, \code{'miRNA_mature'}, \code{'gene'}, \code{'species'}
#' @export
#' @examples mir_info_demo <- demo_mir_info()

demo_mir_info <- function() {
    miRNA <- c("hsa-miR-1", "hsa-miR-1", "hsa-miR-1", "hsa-miR-1", "mmu-miR-105", "rno-miR-106b")
    miRNA_mature <- c("hsa-miR-1-5p", "hsa-miR-1-5p", "hsa-miR-1-3p", "hsa-miR-1-3p", "mmu-miR-105", "rno-miR-106b-5p")
    gene <- c("MIR1-1", "MIR1-2","MIR1-1", "MIR1-2", "Mir105", "Mir106b")
    species <- c("Human", "Human", "Human", "Human", "Mouse", "Rat")
    mir_info_demo <- data.frame(miRNA = miRNA, miRNA_mature = miRNA_mature, gene = gene, species = species, stringsAsFactors = FALSE)
    return(mir_info_demo)
}

#' @title Demo data of mir2tar
#'
#' @description Demo data of mir2tar
#' @details \code{mir2tar} must be a \code{data.frame} object with four columns, namely \code{'miRNA'}, \code{'miRNA_mature'}, \code{'target_gene'}, \code{'species'}
#' @export
#' @examples mir2tar_demo <- demo_mir2tar()

demo_mir2tar <- function() {
    miRNA <- c("hsa-miR-1", "hsa-miR-1", "mmu-miR-105", "rno-miR-106b")
    miRNA_mature <- c("hsa-miR-1-5p", "hsa-miR-1-3p", "mmu-miR-105", "rno-miR-106b-5p")
    target_gene <- c("BDNF", "RBM28", "Abl2", "Mcl1")
    species <- c("Human", "Human", "Mouse", "Rat")
    mir2tar_demo <- data.frame(miRNA = miRNA, miRNA_mature = miRNA_mature, target_gene = target_gene, species = species, stringsAsFactors = FALSE)
    return(mir2tar_demo)
}

#' @title Demo data of pathways
#'
#' @description Demo data of pathways
#' @details \code{pathways} must be a \code{data.frame} object with four columns, namely \code{'src'}, \code{'dest'}, \code{'pathway'}, \code{'species'}
#' @export
#' @examples pathways_demo <- demo_pathways()

demo_pathways <- function() {
    src <- c("CDKN1A", "CDKN1A", "CDK2", "Akt1", "Tcirg1")
    dest <- c("CDK2", "CDK4", "TP53", "Atf2", "Ppa1")
    pathway <- c("p53 signaling pathway", "p53 signaling pathway", "p53 signaling pathway",
                 "PI3K-Akt signaling pathway", "Oxidative phosphorylation")
    species <- c("Human", "Human", "Human", "Mouse", "Rat")
    pathways_demo <- data.frame(src = src, dest = dest, pathway = pathway, species = species, stringsAsFactors = FALSE)
    return(pathways_demo)
}
