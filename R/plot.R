#' @title Chord plot of cell-cell communications
#'
#' @description Chord plot of cell-cell communications from senders to receivers with the sum of inferred miRNAs number, EVmiR_score, or score
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types to plot by order. Default is to plot all cell types.
#' @param celltype_color Colors for the cell types, whose length must be equal to \code{celltype}
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param edge_color Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}
#' @param edge_type Types for the edges from the sender cell type. Default is \code{"big.arrow"}. \code{"ellipse"} for ellipse, "triangle" for triangle, "curved" for curved. Details see \code{\link[circlize]{chordDiagram}}
#' @param show_type which type of miRNAs to show, \code{"number"}, \code{"EVmiR_score"}, or \code{"score"} for sum of inferred miRNAs number, EVmiR_score, and MiTI_score, respectively. Default is \code{"number"}
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param text_size Size of text labels. Default is \code{1.5}
#' @param y_scale y_scale to adjust the text. Default is \code{0.1}
#' @param ... parameters pass to \code{\link[circlize]{chordDiagram}}, e.g., link.arr.width, link.arr.length, link.arr.col
#' @import ggplot2 Matrix circlize
#' @importFrom scales hue_pal
#' @importFrom graphics text
#' @return Chord plot of cell-cell communications mediated by EV-derived miRNA
#' @export

plot_miRTalk_chord <- function(object, condition = NULL, celltype = NULL, celltype_color = NULL, miRNA = NULL, edge_color = NULL, edge_type = "big.arrow", show_type = "number", if_show_autocrine = FALSE, text_size = 1.5, y_scale = 0.1, ...) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            print(celltype_raw)
            stop("Please input the right celltype name as shown above!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    celltype <- celltype[order(celltype)]
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
    } else {
        if (length(celltype_color) != length(celltype)) {
            stop("The length of celltype_color must be equal to celltype!")
        }
        clu_col <- celltype_color
    }
    names(clu_col) <- celltype
    if (is.null(edge_color[1])) {
        link_color <- clu_col
    } else {
        if (length(edge_color) != length(celltype)) {
            stop("The length of edge_color must be equal to celltype!")
        }
        link_color <- edge_color
        names(link_color) <- celltype
    }
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these miRNA!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "miRNA")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "EVmiR_score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "EVmiR_score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$EVmiR_score)
        }
        show_type_new <- "EVmiR_score"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("from", "to", "value")
    chordDiagram(x = cci_pair, grid.col = clu_col[celltype], col = link_color[cci_pair$from], preAllocateTracks = 1, transparency = 0.25, directional = 1,
        direction.type = c("arrows", "diffHeight"), diffHeight = -0.04, annotationTrack = "grid", link.arr.type = edge_type, ...)
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + y_scale, sector.name, facing = "inside", niceFacing = FALSE, adj = c(0.5, 0), cex = text_size)
    }, bg.border = NA)
    graphics::text(0, 1, paste0("Show_type: ", show_type_new), cex = text_size)
}

#' @title Circle plot of cell-cell communications
#'
#' @description Circle plot of cell-cell communications from senders to receivers with the sum of inferred miRNAs number, EVmiR_score, or score
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types to plot. Default is to plot all cell types.
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param celltype_color Colors for the cell types, whose length must be equal to \code{celltype}
#' @param edge_color Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}
#' @param edge_type Types for the edges. \code{"fan"} by default, \code{"link"}, \code{"hive"}
#' @param show_type which type of miRNAs to show, \code{"number"}, \code{"EVmiR_score"}, or \code{"score"} for sum of inferred miRNAs number, EVmiR_score, and MiTI_score, respectively. Default is \code{"number"}
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param edge_alpha Transparency of edge. Default is \code{0.5}
#' @param node_size Size of node. Default is \code{10}
#' @param text_size Size of text. Default is \code{5}
#' @import ggplot2 Matrix ggraph
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @return ggplot2 object for Circle plot of cell-cell communications mediated by EV-derived miRNA
#' @export

plot_miRTalk_circle <- function(object, condition = NULL, celltype = NULL, miRNA = NULL, celltype_color = NULL, edge_color = NULL, edge_type = "fan",
    show_type = "number", if_show_autocrine = FALSE, edge_alpha = 0.5, node_size = 10, text_size = 5) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            print(celltype_raw)
            stop("Please input the right celltype name as shown above!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    celltype <- celltype[order(celltype)]
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
    } else {
        if (length(celltype_color) != length(celltype)) {
            stop("The length of celltype_color must be equal to celltype")
        }
        clu_col <- celltype_color
    }
    names(clu_col) <- celltype
    if (is.null(edge_color[1])) {
        link_color <- clu_col
    } else {
        if (length(edge_color) != length(celltype)) {
            stop("The length of edge_color must be equal to celltype!")
        }
        link_color <- edge_color
        names(link_color) <- celltype
    }
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these miRNA!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "miRNA")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "EVmiR_score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "EVmiR_score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$EVmiR_score)
        }
        show_type_new <- "EVmiR_score"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("from", "to", "value")
    cci_pair$sender <- as.character(cci_pair$from)
    celltype_node <- data.frame(name = celltype, celltype = celltype, id = 1:length(celltype), stringsAsFactors = FALSE)
    # angle
    angle <- 360 * (celltype_node$id - 0.5)/nrow(celltype_node)
    celltype_node$angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)
    celltype_node$hjust <- ifelse(angle > 180, 1, 0)
    celltype_ordered <- celltype[order(celltype)]
    clu_col <- clu_col[celltype_ordered]
    celltype_ordered <- celltype_ordered[celltype_ordered %in% cci_pair$sender]
    link_color <- link_color[celltype_ordered]
    cci_pair_new <- data.frame()
    for (i in 1:length(celltype_ordered)) {
        cci_pair1 <- cci_pair[cci_pair$sender == celltype_ordered[i],]
        cci_pair_new <- rbind(cci_pair_new, cci_pair1)
    }
    mygraph <- graph_from_data_frame(cci_pair_new, vertices = celltype_node, directed = FALSE)
    p <- ggraph(mygraph, layout = "linear", circular = TRUE)
    if (edge_type == "fan") {
        p <- p + geom_edge_fan(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (edge_type == "link") {
        p <- p + geom_edge_link(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (edge_type == "hive") {
        p <- p + geom_edge_hive(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (if_show_autocrine) {
        p <- p + geom_edge_loop(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    p + geom_node_point(aes(color = celltype), size = node_size) + scale_color_manual(values = clu_col) +
        geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = name, angle = angle, hjust = hjust), size = text_size) +
        scale_edge_color_manual(values = link_color) + expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6)) +
        coord_fixed() + theme_minimal() + theme(legend.position = "none", panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"))
}

#' @title Circle plot of cell-cell communications by retaining all cell type nodes
#'
#' @description Circle plot of cell-cell communications from senders to receivers with the sum of inferred miRNAs number, EVmiR_score, or score by retaining all cell type nodes
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types to plot. one or more cell types
#' @param celltype_dir which direction to plot, \code{"sender"} or \code{"receiver"}. Default is as \code{"sender"}.
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param celltype_color Colors for the cell types, whose length must be equal to \code{celltype}
#' @param edge_color Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}
#' @param edge_type Types for the edges. \code{"fan"} by default, \code{"link"}, \code{"hive"}
#' @param show_type which type of miRNAs to show, \code{"number"}, \code{"EVmiR_score"}, or \code{"score"} for sum of inferred miRNAs number, EVmiR_score, and MiTI_score, respectively. Default is \code{"number"}
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param edge_alpha Transparency of edge. Default is \code{0.5}
#' @param node_size Size of node. Default is \code{10}
#' @param text_size Size of text. Default is \code{5}
#' @import ggplot2 Matrix ggraph
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @return ggplot2 object for Circle plot of cell-cell communications mediated by EV-derived miRNA
#' @export

plot_miRTalk_circle_simple <- function(object, condition = NULL, celltype, celltype_dir = "sender", miRNA = NULL, celltype_color = NULL, edge_color = NULL, edge_type = "fan",
    show_type = "number", if_show_autocrine = FALSE, edge_alpha = 0.5, node_size = 10, text_size = 5) {
    # check input
    if (!is(object, "miRTalk")) {
      stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
      stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (is.null(celltype[1])) {
        stop("Please input the celltype!")
    } else {
        if (!all(celltype %in% celltype_raw)) {
            print(celltype_raw)
            stop("Please input the right celltype name as shown above!")
        }
        if (celltype_dir == "sender") {
            cci <- cci[cci$celltype_sender %in% celltype, ]
        } else {
            cci <- cci[cci$celltype_receiver %in% celltype, ]
        }
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    celltype_raw <- celltype_raw[order(celltype_raw)]
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype_raw))
    } else {
        if (length(celltype_color) != length(celltype_raw)) {
            stop("The length of celltype_color must be equal to all celltype in object@cci")
        }
        clu_col <- celltype_color
    }
    names(clu_col) <- celltype_raw
    if (is.null(edge_color[1])) {
        link_color <- clu_col
    } else {
        if (length(edge_color) != length(celltype_raw)) {
            stop("The length of edge_color must be equal to celltype!")
        }
        link_color <- edge_color
        names(link_color) <- celltype_raw
    }
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "miRNA")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "EVmiR_score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "EVmiR_score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$EVmiR_score)
        }
        show_type_new <- "EVmiR_score"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("from", "to", "value")
    cci_pair$sender <- as.character(cci_pair$from)
    celltype_node <- data.frame(name = celltype_raw, celltype = celltype_raw, id = 1:length(celltype_raw), stringsAsFactors = FALSE)
    # angle
    angle <- 360 * (celltype_node$id - 0.5)/nrow(celltype_node)
    celltype_node$angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)
    celltype_node$hjust <- ifelse(angle > 180, 1, 0)
    celltype_ordered <- celltype_raw[order(celltype_raw)]
    clu_col <- clu_col[celltype_ordered]
    celltype_ordered <- celltype_ordered[celltype_ordered %in% cci_pair$sender]
    link_color <- link_color[celltype_ordered]
    cci_pair_new <- data.frame()
    for (i in 1:length(celltype_ordered)) {
        cci_pair1 <- cci_pair[cci_pair$sender == celltype_ordered[i],]
        cci_pair_new <- rbind(cci_pair_new, cci_pair1)
    }
    mygraph <- graph_from_data_frame(cci_pair_new, vertices = celltype_node, directed = FALSE)
    p <- ggraph(mygraph, layout = "linear", circular = TRUE)
    if (edge_type == "fan") {
        p <- p + geom_edge_fan(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (edge_type == "link") {
        p <- p + geom_edge_link(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (edge_type == "hive") {
        p <- p + geom_edge_hive(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (if_show_autocrine) {
        p <- p + geom_edge_loop(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    p + geom_node_point(aes(color = celltype_raw), size = node_size) + scale_color_manual(values = clu_col) +
        geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = name, angle = angle, hjust = hjust), size = text_size) +
        scale_edge_color_manual(values = link_color) + expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6)) +
        coord_fixed() + theme_minimal() + theme(legend.position = "none", panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"))
}

#' @title Sankey plot of cell-cell communications
#'
#' @description Sankey plot of cell-cell communications from senders to receivers with the sum of inferred miRNAs number, EVmiR_score, or score
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types to plot. Default is to plot all cell types
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param celltype_color Colors for the cell types, whose length must be equal to \code{celltype}
#' @param edge_color Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}, Or use \code{"NO"} to cancel it
#' @param show_type which type of miRNAs to show, \code{"number"}, \code{"EVmiR_score"}, or \code{"score"} for sum of inferred miRNAs number, EVmiR_score, and MiTI_score, respectively. Default is \code{"number"}
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param edge_alpha Transparency of edge. Default is \code{0.5}
#' @param node_size Size of node. Default is \code{40}
#' @param text_size Size of text. Default is \code{15}
#' @param node_pad Size of node padding. Numeric essentially influences the width height. Default is \code{20}
#' @param ... parameters pass to \code{\link[networkD3]{sankeyNetwork}}
#' @import networkD3 Matrix
#' @importFrom scales hue_pal
#' @return Sankey plot of cell-cell communications mediated by EV-derived miRNA
#' @export

plot_miRTalk_sankey <- function(object, condition = NULL, celltype = NULL, miRNA = NULL, celltype_color = NULL, edge_color = NULL, show_type = "number", if_show_autocrine = FALSE,
    edge_alpha = 0.5, node_size = 40, text_size = 15, node_pad = 20, ...) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            print(celltype_raw)
            stop("Please input the right celltype name as shown above!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    celltype <- celltype[order(celltype)]
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
    } else {
        if (length(celltype_color) != length(celltype)) {
            stop("The length of celltype_color must be equal to celltype")
        }
        clu_col <- celltype_color
    }
    names(clu_col) <- celltype
    if (is.null(edge_color[1])) {
        link_color <- clu_col
        names(link_color) <- celltype
        edge_color_new <- TRUE
    } else {
        if (edge_color[1] == "NO") {
            edge_color_new <- FALSE
            link_color <- clu_col
            names(link_color) <- celltype
        } else {
            if (length(edge_color) != length(celltype)) {
                stop("The length of edge_color must be equal to celltype!")
            }
            edge_color_new <- TRUE
            link_color <- edge_color
            names(link_color) <- celltype
        }
    }
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these miRNA!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "miRNA")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "EVmiR_score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "EVmiR_score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$EVmiR_score)
        }
        show_type_new <- "EVmiR_score"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("from", "to", "value")
    cci_pair$signal <- cci_pair$from
    cci_pair$fromID <- 0
    cci_pair$toID <- 0
    celltype_sender <- data.frame(name = unique(cci_pair$from), stringsAsFactors = FALSE)
    celltype_sender$id <- (1:nrow(celltype_sender)) - 1
    celltype_receiver <- data.frame(name = unique(cci_pair$to), stringsAsFactors = FALSE)
    celltype_receiver$id <- (1:nrow(celltype_receiver)) - 1 + nrow(celltype_sender)
    celltype_node <- rbind(celltype_sender, celltype_receiver)
    # sender
    for (i in 1:nrow(celltype_sender)) {
        cci_pair[cci_pair$from == celltype_sender$name[i], ]$fromID <- celltype_sender$id[i]
    }
    # receiver
    for (i in 1:nrow(celltype_receiver)) {
        cci_pair[cci_pair$to == celltype_receiver$name[i], ]$toID <- celltype_receiver$id[i]
    }
    # link color
    celltype_sender <- celltype[celltype %in% cci_pair$from]
    cci_pair$link_type <- paste0("type", 1:nrow(cci_pair))
    cci_pair$link_color <- link_color[cci_pair$from]
    # node color
    celltype_node$node_color <- clu_col[celltype_node$name]
    ColourScal <- 'd3.scaleOrdinal() .domain(['
    cci_pair_link <- unique(cci_pair[,c("link_type","link_color")])
    cci_pair_link_type <- cci_pair_link$link_type
    cci_pair_link_color <- cci_pair_link$link_color
    celltype_node_new <- unique(celltype_node[,c("name","node_color")])
    celltype_node_type <- celltype_node_new$name
    celltype_node_color <- celltype_node_new$node_color
    cci_pair_link_type <- c(cci_pair_link_type, celltype_node_type)
    for (i in 1:length(cci_pair_link_type)) {
        if (i != length(cci_pair_link_type)) {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_type[i], '",')
        } else {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_type[i], '"])')
        }
    }
    ColourScal <- paste0(ColourScal, " .range([")
    cci_pair_link_color <- c(cci_pair_link_color, celltype_node_color)
    for (i in 1:length(cci_pair_link_color)) {
        if (i != length(cci_pair_link_color)) {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_color[i], '",')
        } else {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_color[i], '"])')
        }
    }
    if (edge_color_new) {
        sankeyNetwork(Links = cci_pair, Nodes = celltype_node, Source = "fromID", Target = "toID", Value = "value", NodeID = "name",colourScale = ColourScal,
            nodeWidth = node_size, fontSize = text_size, nodePadding = node_pad, fontFamily = "Arial", sinksRight = FALSE, LinkGroup = "link_type", NodeGroup = "name", ...)
    } else {
        sankeyNetwork(Links = cci_pair, Nodes = celltype_node, Source = "fromID", Target = "toID", Value = "value", NodeID = "name",colourScale = ColourScal,
            nodeWidth = node_size, fontSize = text_size, nodePadding = node_pad, fontFamily = "Arial", sinksRight = FALSE, NodeGroup = "name", ...)
    }
}

#' @title Heatmap plot of cell-cell communications
#'
#' @description Heatmap plot of cell-cell communications from senders to receivers with the sum of inferred miRNAs number, EVmiR_score, or score
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types to plot by order. Default is to plot all cell types
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param show_type which type of miRNAs to show, \code{"number"}, \code{"EVmiR_score"}, or \code{"score"} for sum of inferred miRNAs number, EVmiR_score, and MiTI_score, respectively. Default is \code{"number"}
#' @param text_size Size of text labels. Default is \code{10}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @param ... parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color, grid_width
#' @import heatmaply viridis reshape2 Matrix
#' @return Heatmap plot of cell-cell communications mediated by EV-derived miRNA
#' @export

plot_miRTalk_heatmap <- function(object, condition = NULL, celltype = NULL, miRNA = NULL, show_type = "number", text_size = 10, viridis_option = "D", ...) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            print(celltype_raw)
            stop("Please input the right celltype name as shown above!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    celltype <- celltype[order(celltype)]
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these miRNA!")
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "miRNA")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "EVmiR_score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "EVmiR_score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$EVmiR_score)
        }
        show_type_new <- "EVmiR_score"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("sender", "receiver", "value")
    cci_pair <- reshape2::dcast(data = cci_pair, formula = sender ~ receiver, value.var = "value", fill = 0)
    rownames(cci_pair) <- cci_pair$sender
    cci_pair <- cci_pair[, -1]
    heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = viridis_option)
    heatmaply::heatmaply(x = as.matrix(cci_pair), colors = heat_col, dendrogram = "none", xlab = "Receiver", ylab = "Sender", main = paste0("Show type: ", show_type_new), margins = c(60,100,40,20),
        branches_lwd = 0.1, fontsize_row = text_size, fontsize_col = text_size, labCol = colnames(cci_pair), labRow = rownames(cci_pair),
        heatmap_layers = theme(axis.line=element_blank()), label_names = c("Sender","Receiver",show_type_new), ...)
}

#' @title Heatmap plot of inferred miRNA
#'
#' @description heatmap plot of inferred miRNA for each sender. Rows for cell types, and columns for miRNAs by default
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types to plot. Default is to plot all cell types
#' @param miRNA which miRNAs to plot. Default is to plot all inferred miRNAs in senders.
#' @param text_size Size of text labels. Default is \code{10}
#' @param if_horizontal Whether to plot with the horizontal direction. Default is \code{TRUE}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @param ... parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color, grid_width
#' @import heatmaply Matrix
#' @importFrom reshape2 dcast
#' @return Heatmap plot of inferred miRNA
#' @export

plot_miR_heatmap <- function(object, condition = NULL, celltype = NULL, miRNA = NULL, text_size = 10, if_horizontal = TRUE, viridis_option = "D", ...) {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- get_miRTalk_cci(object)
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_sender <- unique(cci$celltype_sender)
    miR_name <- unique(cci$miRNA)
    if (!is.null(celltype[1])) {
        if (!all(celltype %in% celltype_sender)) {
            print(celltype_sender)
            stop("Please input the right celltype name as shown above!")
        }
        cci <- cci[cci$celltype_sender %in% celltype, ]
    }
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    cci <- unique(cci[, c("celltype_sender", "miRNA", "EVmiR_score")])
    if (length(unique(cci$celltype_sender)) < 2 | length(unique(cci$miRNA)) < 2) {
        stop("Limited celltype_sender and miRNA to plot!")
    }
    cci <- reshape2::dcast(data = cci, formula = celltype_sender ~ miRNA, fun.aggregate = mean, value.var = "EVmiR_score", fill = 0)
    rownames(cci) <- cci$celltype_sender
    cci <- cci[, -1]
    if (!if_horizontal) {
        cci <- as.data.frame(t(cci))
    }
    heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = viridis_option)
    if(if_horizontal) {
        heatmaply::heatmaply(x = as.matrix(cci), colors = heat_col, dendrogram = "none", xlab = "miRNA", ylab = "Senders", main = "EVmiR_score in senders", margins = c(60,100,40,20),
            titleX = FALSE, branches_lwd = 0.1, fontsize_row = text_size,
            fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci), heatmap_layers = theme(axis.line=element_blank()), label_names = c("Sender","miRNA","EVmiR_score"), ...)
    } else {
        heatmaply::heatmaply(x = as.matrix(cci), colors = heat_col, dendrogram = "none", xlab = "Senders", ylab = "miRNA", main = "EVmiR_score in senders", margins = c(60,100,40,20),
            titleX = FALSE, branches_lwd = 0.1, fontsize_row = text_size,
            fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci), heatmap_layers = theme(axis.line=element_blank()), label_names = c("miRNA","Sender","EVmiR_score"), ...)
    }
}

#' @title Heatmap plot of inferred targets in receivers
#'
#' @description heatmap plot of inferred targets in receivers. Rows for cell-type-specific miRNAs, and columns for targets in receivers by default
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types as receivers to plot, one or more cell types.
#' @param targetgenes which targetgenes to plot. Default is to plot all inferred target genes in receivers.
#' @param limits A parameter \code{\link[heatmaply]{heatmaply}}, a two dimensional numeric vector specifying the data range for the scale. Default is 0-1
#' @param text_size Size of text labels. Default is \code{10}
#' @param if_horizontal Whether to plot with the horizontal direction. Default is \code{TRUE}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @param ... parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color, grid_width
#' @import heatmaply Matrix
#' @importFrom reshape2 dcast
#' @return Heatmap plot of inferred targets
#' @export

plot_target_heatmap <- function(object, condition = NULL, celltype, targetgenes = NULL, limits = c(0,1), text_size = 10, if_horizontal = TRUE, viridis_option = "D", ...) {
    # get simple output
    cci_simple <- get_miRTalk_cci(object)
    if (nrow(cci_simple) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_sender <- unique(cci_simple$celltype_sender)
    celltype_receiver <- unique(cci_simple$celltype_receiver)
    if (!is.null(celltype[1])) {
        if (!all(celltype %in% celltype_receiver)) {
            print(celltype_receiver)
            stop("Please input the right celltype name as shown above!")
        }
        cci_simple <- cci_simple[cci_simple$celltype_receiver %in% celltype,]
    }
    if (!is.null(targetgenes)) {
        cci_simple <- cci_simple[cci_simple$target_gene %in% targetgenes,]
    }
    cci_simple$celltype_sender <- paste0(cci_simple$celltype_sender,":", cci_simple$miRNA)
    cci_simple$target_gene <- paste0(cci_simple$celltype_receiver,":", cci_simple$target_gene)
    cci <- unique(cci_simple[, c("celltype_sender", "target_gene", "specifity")])
    if (length(unique(cci$celltype_sender)) < 2 | length(unique(cci$target_gene)) < 2) {
        stop("Limited celltype_sender and target_gene to plot!")
    }
    cci <- reshape2::dcast(data = cci, formula = celltype_sender ~ target_gene, fun.aggregate = mean, value.var = "specifity", fill = 0)
    rownames(cci) <- cci$celltype_sender
    cci <- cci[, -1]
    if (!if_horizontal) {
      cci <- as.data.frame(t(cci))
    }
    heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = viridis_option)
    if(if_horizontal) {
        heatmaply::heatmaply(x = as.matrix(cci), colors = heat_col, dendrogram = "none", xlab = paste0("Target genes in ",celltype), ylab = "Senders", main = "Specifity of miRNAs in senders on target genes", margins = c(60,100,40,20),
            titleX = FALSE, branches_lwd = 0.1, fontsize_row = text_size,
            fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci), heatmap_layers = theme(axis.line=element_blank()), label_names = c("Sender", "Receiver", "Specifity"), limits = limits, ...)
    } else {
        heatmaply::heatmaply(x = as.matrix(cci), colors = heat_col, dendrogram = "none", xlab = "Senders", ylab = paste0("Target genes in ",celltype), main = "Specifity of miRNAs in senders on target genes", margins = c(60,100,40,20),
            titleX = FALSE, branches_lwd = 0.1, fontsize_row = text_size,
            fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci), heatmap_layers = theme(axis.line=element_blank()), label_names = c("Receiver", "Sender", "Specifity"), limits = limits, ...)
    }
}

#' @title Bubble plot of inferred miRNA
#'
#' @description Bubble plot of inferred miRNA from senders top receivers. Rows for cell pairs, and columns for miRNAs by default.
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype which cell types to plot. Default is to plot all cell types
#' @param miRNA which miRNAs to plot. Default is to plot all inferred miRNAs in senders.
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param if_horizontal Whether to plot with the horizontal direction. Default is \code{TRUE}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @import ggplot2 Matrix
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @return ggplot2 object for Bubble plot of inferred miRNA
#' @export

plot_miR_bubble <- function(object, condition = NULL, celltype = NULL, miRNA = NULL, if_show_autocrine = FALSE, if_horizontal = TRUE, viridis_option = "D") {
    # check input
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- get_miRTalk_cci(object)
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    celltype_sender <- unique(cci$celltype_sender)
    if (!is.null(celltype[1])) {
        if (!all(celltype %in% celltype_sender)) {
            print(celltype_sender)
            stop("Please input the right celltype name!")
        }
        cci <- cci[cci$celltype_sender %in% celltype, ]
    }
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    cci$cellpair <- paste0(cci$celltype_sender, " | ", cci$celltype_receiver)
    cci <- unique(cci[, c("cellpair", "miRNA", "EVmiR_score", "score")])
    if (length(unique(cci$cellpair)) < 2 | length(unique(cci$miRNA)) < 2) {
        stop("Limited cellpair and miRNA to plot!")
    }
    cci <- .get_bubble(cci)
    cellpair <- unique(cci$cellpair)
    y_len <- length(cellpair)
    miR_name <- unique(cci$miRNA)
    x_len <- length(miR_name)
    cci <- .get_coord(cci, cellpair, y_len, miR_name, x_len)
    if (if_horizontal) {
        ggplot2::ggplot(data = cci) + ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = score,
            size = EVmiR_score)) + viridis::scale_colour_viridis(option = viridis_option) + ggplot2::scale_y_continuous(name = "Senders | Receivers",
            breaks = 1:y_len, labels = cellpair, limits = c(1, y_len)) + ggplot2::scale_x_continuous(name = "EV-derived miRNAs",
            breaks = 1:x_len, labels = miR_name, limits = c(1, x_len)) + ggplot2::labs(color = "max_score",
            size = "EVmiR_score") + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    } else {
        ggplot2::ggplot(data = cci) + ggplot2::geom_point(ggplot2::aes(x = y, y = x, color = score,
            size = EVmiR_score)) + viridis::scale_colour_viridis(option = viridis_option) + ggplot2::scale_x_continuous(name = "Senders | Receivers",
            breaks = 1:y_len, labels = cellpair, limits = c(1, y_len)) + ggplot2::scale_y_continuous(name = "EV-derived miRNAs",
            breaks = 1:x_len, labels = miR_name, limits = c(1, x_len)) + ggplot2::labs(color = "max_score",
            size = "EVmiR_score") + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
}

#' @title Chord plot of EV-derived miRNAs and target genes
#'
#' @description Chord plot of EV-derived miRNAs and target genes from senders to receivers with communication score
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype_sender Name of celltype_sender. One or more cell types
#' @param celltype_receiver Name of celltype_receiver. One or more cell types
#' @param celltype_color Colors for the celltype_sender nodes and celltype_receiver nodes, or use \code{"NO"} to make it simple
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param edge_color Colors for the edges from the sender cell type
#' @param edge_type Types for the edges from the sender cell type. Default is \code{"circle"}. \code{"big.arrow"} for big arrow, "triangle" for triangle, "ellipse" for ellipse, "curved" for curved. Details see \code{\link[circlize]{chordDiagram}}
#' @param text_size Size of text labels. Default is \code{0.5}
#' @param y_scale y_scale to adjust the text. Default is \code{1}
#' @param ... parameters pass to \code{\link[circlize]{chordDiagram}}, e.g., link.arr.width, link.arr.length, link.arr.col
#' @import ggplot2 Matrix circlize
#' @importFrom scales hue_pal
#' @importFrom graphics text
#' @return Chord plot of EV-derived miRNAs and target genes
#' @export

plot_miR2tar_chord <- function(object, condition = NULL, celltype_sender, celltype_receiver, celltype_color = NULL, miRNA = NULL, edge_color = NULL, edge_type = "circle", text_size = 0.5, y_scale = 1, ...) {
    # check object
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- get_miRTalk_cci(object)
    celltype <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    # check celltype_sender
    if (!all(celltype_sender %in% cci$celltype_sender)) {
        print(celltype)
        stop("Please provide the correct name of celltype_sender as shown above!")
    }
    # check celltype_receiver
    if (!all(celltype_receiver %in% cci$celltype_receiver)) {
        print(celltype)
        stop("Please provide the correct name of celltype_receiver as shown above!")
    }
    cci <- cci[cci$celltype_sender %in% celltype_sender & cci$celltype_receiver %in% celltype_receiver, ]
    if (nrow(cci) == 0) {
        stop("No miRNA found from senders to receivers!")
    }
    # check miRNA
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    if (nrow(cci) == 0) {
        stop("No miRNA found from senders to receivers!")
    }
    celltype_sender_new <- unique(cci$celltype_sender)
    celltype_receiver_new <- unique(cci$celltype_receiver)
    if (length(celltype_sender) != length(celltype_sender_new)) {
        celltype_sender_no <- celltype_sender[!celltype_sender %in% celltype_sender_new]
        if (length(celltype_sender_no) == 1) {
            warning(paste0("celltype_sender of ", celltype_sender_no, " has been removed!"))
        }
        else {
            celltype_sender_no1 <- celltype_sender_no[1]
            for (i in 2:length(celltype_sender_no)) {
                celltype_sender_no1 <- paste0(celltype_sender_no1, ", ", celltype_sender_no[i])
            }
            warning(paste0("celltype_sender of ", celltype_sender_no1, " have been removed!"))
        }
        celltype_sender <- celltype_sender_new
    }
    if (length(celltype_receiver) != length(celltype_receiver_new)) {
        celltype_receiver_no <- celltype_receiver[!celltype_receiver %in% celltype_receiver_new]
        if (length(celltype_sender_no) == 1) {
            warning(paste0("celltype_receiver of ", celltype_receiver_no, " has been removed!"))
        }
        else {
            celltype_receiver_no1 <- celltype_receiver_no[1]
            for (i in 2:length(celltype_receiver_no)) {
                celltype_receiver_no1 <- paste0(celltype_receiver_no1, ", ", celltype_receiver_no[i])
            }
            warning(paste0("celltype_receiver of ", celltype_receiver_no1, " have been removed!"))
        }
        celltype_receiver <- celltype_receiver_new
    }
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
        names(clu_col) <- celltype
        clu_col <- clu_col[unique(c(celltype_sender, celltype_receiver))]
        celltype_color_new <- FALSE
    } else {
        if (celltype_color[1] == "NO") {
            clu_col <- scales::hue_pal()(length(celltype))
            names(clu_col) <- celltype
            clu_col <- clu_col[unique(c(celltype_sender, celltype_receiver))]
            celltype_color_new <- TRUE
        } else {
            if (length(celltype_color) != length(unique(c(celltype_sender, celltype_receiver)))) {
                stop("The length of celltype_color must be equal to sum of celltype_sender and celltype_receiver!")
            }
            celltype_color_new <- TRUE
            clu_col <- celltype_color
            names(clu_col) <- unique(c(celltype_sender, celltype_receiver))
        }
    }
    if (is.null(edge_color[1])) {
        link_color <- clu_col
        link_color <- link_color[celltype_sender]
    } else {
        if (length(edge_color) != length(celltype_sender)) {
            stop("The length of edge_color must be equal to celltype_sender!")
        }
        link_color <- edge_color
        names(link_color) <- celltype_sender
    }
    link_color <- link_color[cci$celltype_sender]
    cci$miRNA <- paste0(cci$miRNA, " (", cci$celltype_sender, ")")
    cci$target_gene <- paste0(cci$target_gene, " (", cci$celltype_receiver, ")")
    celltype <- c(cci$celltype_sender, cci$celltype_receiver)
    grid_col <- clu_col[celltype]
    names(grid_col) <- c(cci$miRNA, cci$target_gene)
    cci <- unique(cci[ ,c("miRNA", "target_gene", "score")])
    if (celltype_color_new) {
        chordDiagram(x = cci, grid.col = grid_col, col = link_color, preAllocateTracks = 1, transparency = 0.25, directional = 1, direction.type = c("arrows", "diffHeight"),
            diffHeight = -0.04, annotationTrack = "grid", link.arr.type = edge_type, ...)
    } else {
        chordDiagram(x = cci, grid.col = grid_col, col = link_color, preAllocateTracks = 1, transparency = 0.25, directional = 1, direction.type = c("arrows", "diffHeight"),
            diffHeight = -0.04, annotationTrack = "grid", link.arr.type = edge_type, ...)
    }
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1]+0.5, sector.name, facing = "clockwise", niceFacing = TRUE, cex = text_size)
    }, bg.border = NA)
}

#' @title Circle plot of EV-derived miRNAs and target genes
#'
#' @description Chord plot of EV-derived miRNAs and target genes from senders to receivers.
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype_sender Name of celltype_sender. One cell type
#' @param celltype_receiver Name of celltype_receiver. One cell type
#' @param celltype_color Colors for the celltype_sender nodes and celltype_receiver nodes, or use \code{"NO"} to make it simple
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param node_size Size of node. Default is \code{3}
#' @param edge_color Colors for the edges from the sender cell type
#' @param text_size Size of text labels. Default is \code{3}
#' @param edge_width Width of edge. Default is \code{0.5}
#' @param if_show_legend Whether to show legends. Default is FALSE
#' @import Matrix ggraph
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @return ggplot2 object for Circle plot of EV-derived miRNAs and target genes
#' @export

plot_miR2tar_circle <- function(object, condition = NULL, celltype_sender, celltype_receiver, celltype_color = NULL, miRNA = NULL, node_size = 3, edge_color = NULL,
    edge_width = 0.5, text_size = 3, if_show_legend = F) {
    # check object
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- get_miRTalk_cci(object)
    celltype <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    # check celltype_sender
    if (!celltype_sender %in% cci$celltype_sender) {
        print(celltype)
        stop("Please provide the correct name of celltype_sender as shown above!")
    }
    # check celltype_receiver
    if (!celltype_receiver %in% cci$celltype_receiver) {
        print(celltype)
        stop("Please provide the correct name of celltype_receiver as shown above!")
    }
    cci <- cci[cci$celltype_sender == celltype_sender & cci$celltype_receiver == celltype_receiver, ]
    if (nrow(cci) == 0) {
        stop("No miRNA found from senders to receivers!")
    }
    # check miRNA
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
        names(clu_col) <- celltype
        clu_col <- clu_col[c(celltype_sender, celltype_receiver)]
        celltype_color_new <- FALSE
    } else {
        if (celltype_color[1] == "NO") {
            clu_col <- scales::hue_pal()(length(celltype))
            names(clu_col) <- celltype
            clu_col <- clu_col[c(celltype_sender, celltype_receiver)]
            celltype_color_new <- TRUE
        } else {
            if (length(celltype_color) != 2) {
                stop("The length of celltype_color must be equal to 2!")
            }
            celltype_color_new <- TRUE
            clu_col <- celltype_color
        }
    }
    if (is.null(edge_color[1])) {
        link_color <- clu_col[1]
    } else {
        if (length(edge_color) != 1) {
            stop("The length of edge_color must be equal to 1!")
        }
        link_color <- edge_color
    }
    cci <- unique(cci[ ,c("miRNA", "target_gene", "EVmiR_score", "target_gene_activity", "score")])
    cci <- .get_miR2tar_circle(cci)
    target_genes<- as.data.frame(table(cci$target_gene), stringsAsFactors = F)
    if (max(target_genes$Freq) > 1) {
        target_genes <- target_genes[target_genes$Freq > 1, ]
        for (i in 1:nrow(target_genes)) {
            target_genes_tmp <- rep(target_genes$Var1[i], target_genes$Freq[i])
            target_genes_tmp <- paste0(target_genes_tmp, "_", 1:length(target_genes_tmp))
            cci[cci$target_gene == target_genes$Var1[i], ]$target_gene <- target_genes_tmp
        }
    }
    cci_sender <- unique(data.frame(name = cci$miRNA, activity = 1, group = "sender", stringsAsFactors = FALSE))
    cci_new <- data.frame()
    for (i in 1:nrow(cci_sender)) {
        cci1 <- cci[cci$miRNA == cci_sender$name[i], ]
        cci_new <- rbind(cci_new, cci1)
    }
    cci <- cci_new
    cci_receiver <- unique(data.frame(name = cci$target_gene, activity = 1, group = "receiver", stringsAsFactors = FALSE))
    cci_sender$id <- 1:nrow(cci_sender)
    cci_receiver$id <- 1:nrow(cci_receiver)
    cci_node <- rbind(cci_sender, cci_receiver)
    cci_node1 <- data.frame(name = celltype_sender, activity = 1, group = "sender", id = 0, stringsAsFactors = FALSE)
    cci_node <- rbind(cci_node1, cci_node)
    cci_node$miR2tar <- "sender"
    for (i in 1:nrow(cci_sender)) {
        cci1 <- cci[cci$miRNA == cci_sender$name[i], ]
        cci_node[cci_node$name %in% cci1$target_gene, ]$miR2tar <- cci_sender$name[i]
    }
    rownames(cci_node) <- cci_node$name
    cci <- cci[ ,c("miRNA","target_gene")]
    cci1 <- data.frame(miRNA = celltype_sender, target_gene = cci$miRNA, stringsAsFactors = FALSE)
    cci<- rbind(cci1, cci)
    colnames(cci)[1:2] <- c("from", "to")
    mygraph <- igraph::graph_from_data_frame(cci, vertices = cci_node)
    if (celltype_color_new) {
        p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
            geom_edge_diagonal(color = link_color, width = edge_width, alpha = 0.8) +
            geom_node_point(aes(x = x*1.07, y = y*1.07, color = group), size = node_size) +
            scale_color_manual(values = as.character(clu_col)[c(2,1)]) +
            geom_node_text(aes(x = x*1.25, y = y*1.25,label = name), size = text_size) +
            expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) + theme_void()
        if (if_show_legend) {
            p
        } else {
            p+NoLegend()
        }
    } else {
        p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
            geom_edge_diagonal(color = link_color, width = edge_width) +
            geom_node_point(aes(x = x*1.07, y = y*1.07, color = miR2tar), size = node_size) +
            geom_node_text(aes(x = x*1.25, y = y*1.25, label = name), size = text_size) +
            expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) + theme_void()
        if (if_show_legend) {
            p
        } else {
            p+NoLegend()
        }
    }
}

#' @title Heatmap plot of EV-derived miRNAs and target genes
#'
#' @description Heatmap plot of EV-derived miRNAs and target genes from senders to receivers with communication score displayed
#' @param object miRTalk object after \code{\link{find_miRTalk}}
#' @param condition which conditions to plot. Default is plot all conditions.
#' @param celltype_sender Name of celltype_sender
#' @param celltype_receiver Name of celltype_receiver
#' @param miRNA which miRNAs to use. Default is to plot all inferred miRNAs in senders.
#' @param text_size Size of text labels. Default is \code{3}
#' @param if_horizontal Whether to plot with the horizontal direction. Default is \code{TRUE}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @param ... parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color
#' @import pheatmap Matrix
#' @importFrom reshape2 dcast
#' @return Heatmap plot of EV-derived miRNAs and target genes
#' @export

plot_miR2tar_heatmap <- function(object, condition = NULL, celltype_sender, celltype_receiver, miRNA = NULL, text_size = 5, if_horizontal = TRUE, viridis_option = "D", ...) {
    # check object
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    cci <- get_miRTalk_cci(object)
    celltype <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    if (!is.null(condition)) {
        cci_condition <- unique(cci$condition)
        if (!all(condition %in% cci_condition)) {
            print(condition)
            stop("Please input the right condition as shown above!")
        }
        cci <- cci[cci$condition %in% condition, ]
    }
    # check celltype_sender
    if (!celltype_sender %in% cci$celltype_sender) {
        print(celltype)
        stop("Please provide the correct name of celltype_sender as shown above!")
    }
    # check celltype_receiver
    if (!celltype_receiver %in% cci$celltype_receiver) {
        print(celltype)
        stop("Please provide the correct name of celltype_receiver as shown above!")
    }
    cci <- cci[cci$celltype_sender == celltype_sender & cci$celltype_receiver == celltype_receiver, ]
    if (nrow(cci) == 0) {
        stop("No miRNA found from senders to receivers!")
    }
    # check miRNA
    miR_name <- unique(cci$miRNA)
    if (!is.null(miRNA[1])) {
        if (!all(miRNA %in% miR_name)) {
            stop("Please input the right miRNA name!")
        }
        cci <- cci[cci$miRNA %in% miRNA, ]
    }
    cci <- unique(cci[ ,c("miRNA", "target_gene", "score")])
    cci <- reshape2::dcast(data = cci, formula = miRNA ~ target_gene, fun.aggregate = mean, value.var = "score", fill = 0)
    main_title <- "score"
    rownames(cci) <- cci$miRNA
    cci <- cci[, -1]
    cci <- as.matrix(cci)
    cci[which(cci == 0)] <- NA
    if (!if_horizontal) {
        cci <- as.data.frame(t(cci))
    }
    if (if_horizontal) {
        heatmaply::heatmaply(x = as.matrix(cci), dendrogram = "none", xlab = paste0(celltype_receiver, ": target"), ylab = paste0(celltype_sender, ": miRNA"), main = main_title,
            fontsize_row = text_size, na.value = "white", fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci),
            heatmap_layers = theme(axis.line=element_blank()), label_names = c("miRNA","target",main_title), ...)
    } else{
        heatmaply::heatmaply(x = as.matrix(cci), dendrogram = "none", xlab = paste0(celltype_sender, ": miRNA"), ylab = paste0(celltype_receiver, ": target"), main = main_title,
            fontsize_row = text_size, na.value = "white", fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci),
            heatmap_layers = theme(axis.line=element_blank()), label_names = c("target","miRNA",main_title), ...)
    }
}
