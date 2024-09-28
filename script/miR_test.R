# data: sc data or spot data matrix
# species, "Human", "Mouse", "Rat"
# mir_info

miR_test <- function(data, species, mir_info){
    mir_info1 <- mir_info[mir_info$species == species,]
    mir_gene <- unique(mir_info1$gene)
    mir_gene <- mir_gene[mir_gene %in% rownames(data)]
    if (length(mir_gene) >= 5) {
        data1 <- data[mir_gene,]
        expcells <- apply(data1, 1, function(x){
            length(x[x>0])/length(x)
        })
        expcells <- expcells[expcells > 0.01]
        if (length(mir_gene) >= 5) {
            return("OKK!")
        } else{
            return("No miRNA!")
        }
    } else {
        return("No miRNA!")
    }
}





