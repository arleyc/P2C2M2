read.sptree.samples2 <- function(file) {
    
    ##### parts of this function were adapted from a function in package 'phyloch',
    ##### written by Christoph Heibl read in trees without stats
   
    tree2 <- read.nexus(file)
    many2 <- class(tree2) == "multiPhylo"
    
    ##### get tree strings
    X2 <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    X2 <- X2[grep("tree STATE", X2)]
    X2 <- gsub("tree STATE_[[:digit:]]+[[:space:]]= ", "", X2)

    ###### a function to produce branch width vectors and order them appropriately
    
    extract.order.stats2 <- function(treestring) {
        
        ntax <- length(unlist(strsplit(treestring, "\\[")))/2
        nodes <- (ntax + 1):((2 * ntax) - 1)
        Y <- treestring
        
        #### adds internal node labels
        for (i in 1:length(nodes)) {
            repl <- paste(")", nodes[i], "[", sep = "")
            Y <- sub(")\\[", repl, Y)
        }
        
        meta <- unlist(strsplit(Y, "\\[|\\]"))[grep("dmv", unlist(strsplit(Y, "\\[|\\]")))]
        metacols <- length(unlist(strsplit(meta[1], ",")))
        meta <- gsub("&|dm.=|\\{|\\}", "", meta)
        
        Ysub <- gsub("\\[[^]]*\\]", "\\[\\]", Y)
        Ysub <- unlist(strsplit(Ysub, ",|)"))
        Ysub <- gsub("\\(|\\)|;|\\[|\\]", "", Ysub)
        
        # to read sptrees from Beast 2.6
        # Ysub[length(Ysub)] <- paste(Ysub[length(Ysub)], NA, sep = ":") 
        
        branchdata <- array(dim = c(length(meta), 2 + metacols))
        
        for (i in 1:length(meta)) {
            branchdata[i, ] <- c(unlist(strsplit(Ysub[i], ":")), unlist(strsplit(meta[i], 
                ",")))
        }
        
        if (metacols == 3) {
            colnames(branchdata) <- c("br", "length", "dmt", "dmv_start", "dmv_end")
        }
        if (metacols == 2) {
            colnames(branchdata) <- c("br", "length", "dmv_start", "dmv_end")
        }

        if (metacols == 1) {
            colnames(branchdata) <- c("br", "length", "dmv")
        }
        rownames(branchdata) <- branchdata[, 1]
        
        string <- gsub("\\[[^]]*\\]", "", Y)
        stree <- read.tree(text = string)
        translate <- cbind(stree$node.label, (ntax + 1):nodes[length(nodes)])
        translate <- rbind(translate, cbind(stree$tip.label, 1:ntax))
        rownames(translate) <- translate[, 2]
        translate2 <- translate[as.character(stree$edge[, 2]), ]
        branchdata2 <- branchdata[translate2[, 1], ]
        
        branchdata2 <- rbind(branchdata2, branchdata[length(branchdata[, 2]), ])
        
        rownames(branchdata2) <- NULL
        return(branchdata2)
            
    }

   branchdata2 <- lapply(X2, extract.order.stats2)
   
    ####### append branch width stats to trees.  vector 'check' is included as a
    ####### verification and should exactly match vector 'edge.length' this is the LONGEST
    ####### part of the function.  can I do this faster???!?

   
   if (length(branchdata2[[1]][1, ]) == 5) {
         if (class(tree2) == "multiPhylo") {
             for (i in 1:length(tree2)) {
                 # tree2[[i]][['check']]<-as.numeric(branchdata2[[i]][,2])
                 tree2[[i]][["dmt"]] <- as.numeric(branchdata2[[i]][, 3])
                 tree2[[i]][["dmv"]] <- as.numeric(branchdata2[[i]][, 4])
                 tree2[[i]][["dmv_start"]] <- as.numeric(branchdata2[[i]][, 4])
                 tree2[[i]][["dmv_end"]] <- as.numeric(branchdata2[[i]][, 5])
             }
         }
         if (class(tree2) == "phylo") {
             # tree[['check']]<-as.numeric(branchdata[[1]][,2])
             tree2[["dmt"]] <- as.numeric(branchdata2[[1]][, 3])
             tree2[["dmv"]] <- as.numeric(branchdata2[[1]][, 4])
             tree2[["dmv_start"]] <- as.numeric(branchdata2[[1]][, 4])
             tree2[["dmv_end"]] <- as.numeric(branchdata2[[1]][, 5])
         }
     }
     
    if (length(branchdata2[[1]][1, ]) == 4) {
         if (class(tree2) == "multiPhylo") {
             for (i in 1:length(tree2)) {
                 # tree2[[i]][['check']]<-as.numeric(branchdata2[[i]][,2])
                 tree2[[i]][["dmv"]] <- as.numeric(branchdata2[[i]][, 3])
                 tree2[[i]][["dmv_start"]] <- as.numeric(branchdata2[[i]][, 3])
                 tree2[[i]][["dmv_end"]] <- as.numeric(branchdata2[[i]][, 4])
             }
         }
         if (class(tree2) == "phylo") {
             # tree[['check']]<-as.numeric(branchdata[[1]][,2])
             tree2[["dmv"]] <- as.numeric(branchdata2[[1]][, 3])
             tree2[["dmv_start"]] <- as.numeric(branchdata2[[1]][, 3])
             tree2[["dmv_end"]] <- as.numeric(branchdata2[[1]][, 4])
         }
     }

     if (length(branchdata2[[1]][1, ]) == 3) {
         if (class(tree2) == "multiPhylo") {
             for (i in 1:length(tree2)) {
                 # tree[[i]][['check']]<-as.numeric(branchdata[[i]][,2])
                 tree2[[i]][["dmv"]] <- as.numeric(branchdata2[[i]][, 3])
             }
         }
         if (class(tree2) == "phylo") {
             # tree[['check']]<-as.numeric(branchdata[[1]][,2])
             tree2[["dmv"]] <- as.numeric(branchdata2[[1]][, 3])
         }
     }
    return(tree2)
} 