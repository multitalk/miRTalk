#' @title Definition of 'miRTalk' class
#'
#' @description An S4 class containing the data, meta, and results of inferred cell-cell communications mediated by EV-derived miRNAs.
#' @slot data A list containing the data and variable genes.
#' @slot meta A data.frame containing the meta data.
#' @slot species A character containing the species.
#' @slot miR A data.frame containing expressed miRNA genes.
#' @slot miR2tar A data.frame containing expressed miRNAs and their target genes.
#' @slot type A character containing the type of miRNA.
#' @slot cci A data.frame containing the significantly enriched EV-derived miRNAs and their target genes.
#' @import methods
#' @name miRTalk
#' @rdname miRTalk
#' @aliases miRTalk-class
#' @exportClass miRTalk

setClass("miRTalk", representation(data = "list", meta = "data.frame", species = "character",
    miR = "data.frame", miR2tar = "data.frame", type = "character", cci = "data.frame"), prototype(data = list(), meta = data.frame(), species = character(),
    miR = data.frame(), miR2tar = data.frame(), type = character(), cci = data.frame()))
