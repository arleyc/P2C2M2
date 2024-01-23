read.gtree.samples2 <- function(file) {

TREE2 <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
   X2 <- TREE2[grep("tree STATE", TREE2)]
   X2 <- gsub("tree STATE_.*= ", "", X2)

## Get translation table
    trans_beg2 <- grep("Translate", TREE2)
    trans_end2 <- grep(";", TREE2)[which(grep(";", TREE2) > trans_beg2)[1]]
	trans_tab2 <- read.table(file=file, skip=trans_beg2, nrows=(trans_end2-trans_beg2-1))
	trans_tab2$Translate <- gsub(",$", "", trans_tab2$Translate)
	
###### a function to produce branch rate vectors and order them appropriately

    extract.order.stats2 <- function(treestring) {

        ntax2 <- length(unlist(strsplit(treestring, "\\["))) /2
        edges2 <- (ntax2 + 2):((2 * ntax2) - 1)

        Y2 <- treestring
        
        #### adds internal node labels
        for (i in 1:length(edges2)) {
            repl2 <- paste(")", edges2[i], "[", sep = "")
            Y2 <- sub(")\\[", repl2, Y2)
        }
        
        meta2 <- unlist(strsplit(Y2, "\\[|\\]"))[grep("rate", unlist(strsplit(Y2, "\\[|\\]")))]
        
        metacols2 <- length(unlist(strsplit(meta2[1], ",")))
        meta2 <- gsub("&|rate=|\\{|\\}", "", meta2)
        meta2<-c(rep("1.0", length(meta2)))
        meta2<-meta2[-length(meta2)]
        
        Z2 <- gsub(";", "", Y2)
        Ysub2 <- gsub("\\[[^]]*\\]", "\\[\\]", Z2)
        Ysub2 <- unlist(strsplit(Ysub2, ",|)"))
        Ysub2 <- gsub("\\(|\\)|;|\\[|\\]", "", Ysub2)
        Ysub2 <- Ysub2[-length(Ysub2)]
        
        brlen2 <- do.call("rbind", strsplit(Ysub2, ":"))
        brrate2 <- do.call("rbind", strsplit(meta2, ","))
        branchdata2 <- cbind(brlen2, brrate2)
        
        if (metacols2 == 1) {
            colnames(branchdata2) <- c("br", "length", "rate")
        }
        rownames(branchdata2) <- branchdata2[, 1]
         
        string2 <- gsub("\\[[^]]*\\]", "", Y2)
        
        stree2 <- ape::read.tree(text = string2)
        
        translate2 <- cbind(c(stree2$node.label[-1], stree2$tip.label),
                            c((ntax2 + 2):edges2[length(edges2)], 1:ntax2))

        rownames(translate2) <- translate2[, 2]
        translate4 <- translate2[as.character(stree2$edge[, 2]), ]
        branchdata4 <- branchdata2[translate4[, 1], ]
        rownames(branchdata4) <- NULL

        stree2$tip.label <- trans_tab2[stree2$tip.label, "Translate"]

        list(tree=stree2, branchdata=branchdata4)
        
        }

        tree_info2 <- lapply(X2, extract.order.stats2)
        tree2 <- lapply(tree_info2, function(x) x$tree)
        branchdata2 <- lapply(tree_info2, function(x) x$branchdata)
      
        ####### append branch rate stats to trees.

        if (length(tree2) > 1) {
           tree2 <- mapply(function(tr2, mat2) {
            mat2 <- mat2[, c(-1, -2), drop=FALSE]
            mat2 <- apply(mat2, 2, function(x) {
                as.numeric(gsub("(.+)\\.([0-9]+\\.[0-9]+E?-?[0-9]?)$",
                                "\\2", x))
            })
            tr2[["rate"]] <- mat2
            tr2
        }, tr=tree2, mat=branchdata2, SIMPLIFY=FALSE)
    } else {
        mode(branchdata2[[1]]) <- "numeric"
        tree2[["rate"]] <- as.numeric(branchdata2[[1]][, c(-1, -2)])
    }
    
    return(tree2)
}
