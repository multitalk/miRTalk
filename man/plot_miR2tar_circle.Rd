% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_miR2tar_circle}
\alias{plot_miR2tar_circle}
\title{Circle plot of EV-derived miRNAs and target genes}
\usage{
plot_miR2tar_circle(
  object,
  condition = NULL,
  celltype_sender,
  celltype_receiver,
  celltype_color = NULL,
  miRNA = NULL,
  node_size = 3,
  edge_color = NULL,
  edge_width = 0.5,
  text_size = 3,
  if_show_legend = F
)
}
\arguments{
\item{object}{miRTalk object after \code{\link{find_miRTalk}}}

\item{condition}{which conditions to plot. Default is plot all conditions.}

\item{celltype_sender}{Name of celltype_sender. One cell type}

\item{celltype_receiver}{Name of celltype_receiver. One cell type}

\item{celltype_color}{Colors for the celltype_sender nodes and celltype_receiver nodes, or use \code{"NO"} to make it simple}

\item{miRNA}{which miRNAs to use. Default is to plot all inferred miRNAs in senders.}

\item{node_size}{Size of node. Default is \code{3}}

\item{edge_color}{Colors for the edges from the sender cell type}

\item{edge_width}{Width of edge. Default is \code{0.5}}

\item{text_size}{Size of text labels. Default is \code{3}}

\item{if_show_legend}{Whether to show legends. Default is FALSE}
}
\value{
ggplot2 object for Circle plot of EV-derived miRNAs and target genes
}
\description{
Chord plot of EV-derived miRNAs and target genes from senders to receivers.
}
