p2c2m2 <- function (xmlfile="beast.xml", combinedfiledirectory, msdir, nreps=1, skip=0, subsample=0, use.mpi=FALSE) {

xmlfile2 <- paste(unlist(c(combinedfiledirectory, "/", xmlfile)), collapse="")
paste(msdir, "ms", sep="/")

read.starbeast2_6 <- function (xmlfile2) {
    
    ##### edit xml before running starbeast2 in order to write branch rate info in gene tree files
    ##### add the following: substitutions="true" branchratemodel="@StrictClock.c:26" (specify branchRateModel id in the dots)
    ##### read in a beast xml file in order to extract the taxon-allele bindings and
    ##### locus names.


######## READ GENE TREES ##########

# source("read.gtree.samples2.R")

read.gtree.samples2 <- function(file) {

TREE2 <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
   X2 <- TREE2[grep("tree STATE", TREE2)]
   X2 <- gsub("tree STATE_.*= ", "", X2)
   
# remove the burn-in trees if skip > 0
if (skip > 0) {
	X2<-X2[-1:-skip]
	}
	
# subsample every subsample-th tree
if (subsample > 0) {
   X2<-X2[1:subsample==subsample]
	}
	
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
        
        library(ape)
        stree2 <- read.tree(text = string2)
        
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

######## READ SPECIES TREE ###############

# source("read.sptree.samples2.R")

read.sptree.samples2 <- function(file) {
    ##### parts of this function were adapted from a function in package 'phyloch',
    ##### written by Christoph Heibl read in trees without stats
   
    tree2 <- read.nexus(file)
    many2 <- class(tree2) == "multiPhylo"
 
# remove the burn-in trees if skip > 0
if (skip > 0) {
	tree2<-tree2[-1:-skip]
	}
   
# subsample every subsample-th tree
if (subsample > 0) {
    tree2 <- tree2[1:subsample==subsample]
    }
    
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
        Ysub[length(Ysub)] <- paste(Ysub[length(Ysub)], NA, sep = ":")
        
        branchdata <- array(dim = c(length(meta), 2 + metacols))
        
        for (i in 1:length(meta)) {
            branchdata[i, ] <- c(unlist(strsplit(Ysub[i], ":")), unlist(strsplit(meta[i], 
                ",")))
        }
        
        if (metacols == 3) {
            colnames(branchdata) <- c("br", "length", "dmt", "dmv_start", "dmv_end")
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

######### END OF READ SPECIES TREE ##############


xmlfile2 <- paste(unlist(c(combinedfiledirectory, "/", xmlfile)), collapse="")

D2 <- scan(file = xmlfile2, what = "", sep = "\n", quiet = TRUE)

D2 <- gsub("\t", "", D2)

Y2 <- D2[grep("<data", D2)[1]:grep("</data>", D2)[length(grep("</data>", D2))]]

Q2 <- D2[grep("<distribution id=\"treeLikelihood", D2)[1]:grep("<branchRateModel.*", D2)[length(grep("<branchRateModel.*", D2))]]

X2 <- D2[grep("<taxon id", D2)[1]:grep("</taxon>", D2)[length(grep("</taxon>", D2))]]

species2 <- X2[grep("<taxon id=.*spec=\"TaxonSet\">", X2)]
species2 <- gsub(".*<taxon id=\"", "", species2)
species2 <- gsub("\" spec=\"TaxonSet\">", "", species2)

alleles2 <- X2[grep("<taxon id=.*spec=\"Taxon\"", X2)]
alleles2 <- gsub(".*<taxon id=\"", "", alleles2)
alleles2 <- gsub("\" spec=\"Taxon\"/>", "", alleles2)

num2 <- (grep("</taxon>", X2) - 1) - grep("<taxon id=.*spec=\"TaxonSet\"", X2)

spec2 <- c()
    for (i in 1:length(species2)) {
        spec2 <- c(spec2, rep(species2[i], times = num2[i]))
    }
    
bindings2 <- cbind(spec2, alleles2)

loci2 <- D2[grep("<distribution id=\"treePrior", D2)]

### ploidy....

ploidy2<-c(rep("NA,length(loci2)"))

for (i in 1:length(loci2)) {
if(length(grep("ploidy", loci2[i]))==1) {
 		locus <- gsub(".*<distribution.*ploidy=\"", "", loci2[i])
 		ploidy <-gsub("\".*>", "", locus)
 	} else {
 		ploidy <-c("2")
   }
ploidy2[i] <- ploidy
}

orderploidy  <- gsub(".*<distribution.*treePrior.t:", "", loci2)
loci2 <- gsub("\" spec.*>", "", orderploidy)
names(ploidy2)  <- loci2
names(loci2) <- loci2

### reading sequence alignments...

p.begin2 <- grep("<distribution id=", Q2)
p.end2 <- grep("</distribution>", Q2)
genes.to.partitions2 <- matrix(ncol = 2, nrow = length(p.begin2))
colnames(genes.to.partitions2) <- c("genes", "partitions")

    for (i in 1:length(p.begin2)) {
        element2 <- Q2[p.begin2[i]]
        gene2 <- gsub(".*<distribution id=\"treeLikelihood.", "", element2)
        gene2 <- gsub("\" spec=.*>", "", gene2)
        part2 <- gsub(".*<distribution.*data=\"@", "", element2)
        part2 <- gsub("\".*>", "", part2)
        genes.to.partitions2[i, ] <- c(gene2, part2)
    }

neworder<-gsub("\"","", gsub(".*id=\"","",Y2[grep("^id=", Y2)]))
rownames(genes.to.partitions2)<-genes.to.partitions2[,1]
genes.to.partitions2 <- genes.to.partitions2[neworder,]
ploidy2 <- ploidy2[neworder]
loci22 <- loci2[neworder]

     alignments2 <- list()

     a.begin2 <- grep("name=\"alignment\"", Y2)
     a.end2 <- grep("</data>", Y2)

for (i in 1:length(genes.to.partitions2[, 1])) {
        Z2 <- Y2[a.begin2[i]:a.end2[i]]
        ID2 <- Z2[grep("taxon=", Z2)]
        ID2 <- gsub(".*<sequence id.*taxon=\"", "", ID2)
        ID2 <- gsub("\" totalcount=.*/>", "", ID2)
        seq2 <- Z2[grep("<sequence id=", Z2)]
        # seq <- as.list(seq)
        names(seq2) <- ID2
        seq2 <- gsub(".*<sequence id=.*value=\"", "", seq2)
        seq2 <- gsub("\"/>", "", seq2)
        seq2 <- sapply(seq2, strsplit, split = "")
        alignments2[[genes.to.partitions2[i, 2]]] <- seq2
    }

cat("bindings, loci, and alignments read\n")
    
    ### read in gene trees from combinedfiledirectory. assumes custom perl script used
    ### to thin and burn in beast results.

library(ape)
    
    gene.trees2 <- list()
    for (i in 1:length(loci2)) {
         treefile2 <- c(combinedfiledirectory, "/", loci2[i], ".trees")
         gene.trees2[[loci2[i]]] <- read.gtree.samples2(paste(unlist(treefile2), collapse = ""))
         cat("locus ", loci2[i], " is done.\n")
    }
    cat("gene trees done\n")   
 
    ### read in species trees from combinedfiledirectory. assumes custom perl script
    ### used to thin and burn in beast results.
    
    species.trees2 <- read.sptree.samples2(paste(unlist(c(combinedfiledirectory, "/", 
        "species.trees")), collapse = ""))
        cat("species trees done\n")
 

### make the various taxon associations.
    empirical <- list()
    for (i in 1:length(loci2)) {
        locus.binding <- c()
        tips <- gene.trees2[[loci2[[i]]]][[1]]$tip.label
        for (j in 1:length(tips)) {
            locus.binding <- rbind(locus.binding, bindings2[which(bindings2[, 2] == 
                tips[j]), ])
        }
        empirical[[loci2[[i]]]] <- as.matrix(as.data.frame(locus.binding))
    }
    
    simulate <- list()
    for (i in 1:length(loci2)) {
        simulate[[loci2[[i]]]] <- as.matrix(as.data.frame(table(empirical[[loci2[[i]]]][, 
            1])))
        # simulate[[loci[[i]]]][,1]<-levels(simulate[[loci[[i]]]][,1])
    }
    
    calculate <- list()
    for (i in 1:length(loci2)) {
        spec <- c()
        for (j in 1:length(simulate[[i]][, 1])) {
            spec <- c(spec, rep(simulate[[i]][j, 1], times = simulate[[i]][j, 2]))
        }
        calculate[[loci2[[i]]]] <- cbind(spec, 1:length(spec))
        calculate[[loci2[[i]]]] <- as.matrix(as.data.frame(calculate[[loci2[[i]]]]))
    }
    
    associations <- list()
    associations[["empirical"]] <- empirical
    associations[["simulate"]] <- simulate
    associations[["calculate"]] <- calculate
    
    #### read in log file
 #   logfile <- c(combinedfiledirectory, "/", logfile)
 #   logdata <- read.table((paste(unlist(logfile), collapse = "")), header = TRUE)
 #   cat("logfile done\n")
    
    #### organize output components
    output <- list()
    output[["associations"]] <- associations
    output[["alignments"]] <- alignments2
    output[["gene.trees"]] <- gene.trees2
    output[["species.trees"]] <- species.trees2
  #  output[["log"]] <- logdata
    output[["genes"]] <- loci2
    output[["genes.to.partitions"]] <- genes.to.partitions2
    output[["ploidy"]] <- as.numeric(ploidy2)
    class(output) <- "starbeast.data"
    return(output)

}
    
output <- read.starbeast2_6()

#################### OUTPUT TO P2C2M ###########################

sbpps_to_p2c2m <- function (output) {

inLwieg = list()

inLwieg$asoc$nme = "ordered allele-species associations"
inLwieg$loci$nme = "loci"
inLwieg$gtre$nme = "gene trees"
inLwieg$pldy$nme = "ploidy level"
inLwieg$ptre$nme = "phybase species trees"
inLwieg$stre$nme = "regular species trees"

inLwieg$asoc$dta[[1]]<-output$associations$empirical$`Full12S+Liolaemus`
colnames(inLwieg$asoc$dta[[1]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[2]]<-output$associations$empirical$`FullA12D+Liolaemus`
colnames(inLwieg$asoc$dta[[2]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[3]]<-output$associations$empirical$`FullA1D+Liolaemus`
colnames(inLwieg$asoc$dta[[3]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[4]]<-output$associations$empirical$`FullA4B+Liolaemus`
colnames(inLwieg$asoc$dta[[4]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[5]]<-output$associations$empirical$`FullA9C+Liolaemus`
colnames(inLwieg$asoc$dta[[5]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[6]]<-output$associations$empirical$`FullCMOS+Liolaemus`
colnames(inLwieg$asoc$dta[[6]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[7]]<-output$associations$empirical$`FullDMXL+Liolaemus`
colnames(inLwieg$asoc$dta[[7]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[8]]<-output$associations$empirical$`FullDNAH3+Liolaemus`
colnames(inLwieg$asoc$dta[[8]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[9]]<-output$associations$empirical$`FullEXPH5+Liolaemus`
colnames(inLwieg$asoc$dta[[9]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[10]]<-output$associations$empirical$`FullKif24+Liolaemus`
colnames(inLwieg$asoc$dta[[10]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[11]]<-output$associations$empirical$`FullMXRA5+Liolaemus`
colnames(inLwieg$asoc$dta[[11]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[12]]<-output$associations$empirical$`FullND4+Liolaemus`
colnames(inLwieg$asoc$dta[[12]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[13]]<-output$associations$empirical$`FullPNN+Liolaemus`
colnames(inLwieg$asoc$dta[[13]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[14]]<-output$associations$empirical$`FullPRLR+Liolaemus`
colnames(inLwieg$asoc$dta[[14]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[15]]<-output$associations$empirical$`FullSNCAIP+Liolaemus`
colnames(inLwieg$asoc$dta[[15]])<-c("speciesName","alleleName")
inLwieg$asoc$dta[[16]]<-output$associations$empirical$`Fullcytb+Liolaemus`
colnames(inLwieg$asoc$dta[[16]])<-c("speciesName","alleleName")


inLwieg$loci$dta<-list("Full12S+Liolaemus","FullA12D+Liolaemus","FullA1D+Liolaemus",
"FullA4B+Liolaemus","FullA9C+Liolaemus","FullCMOS+Liolaemus","FullDMXL+Liolaemus",
"FullDNAH3+Liolaemus","FullEXPH5+Liolaemus","FullKif24+Liolaemus","FullMXRA5+Liolaemus",
"FullND4+Liolaemus","FullPNN+Liolaemus","FullPRLR+Liolaemus","FullSNCAIP+Liolaemus","Fullcytb+Liolaemus")

inLwieg$gtre$dta$`Full12S+Liolaemus` <- output$gene.trees$`Full12S+Liolaemus`
inLwieg$gtre$dta$`FullA12D+Liolaemus` <- output$gene.trees$`FullA12D+Liolaemus`
inLwieg$gtre$dta$`FullA1D+Liolaemus` <- output$gene.trees$`FullA1D+Liolaemus`
inLwieg$gtre$dta$`FullA4B+Liolaemus` <- output$gene.trees$`FullA4B+Liolaemus`
inLwieg$gtre$dta$`FullA9C+Liolaemus` <- output$gene.trees$`FullA9C+Liolaemus`
inLwieg$gtre$dta$`FullCMOS+Liolaemus` <- output$gene.trees$`FullCMOS+Liolaemus`
inLwieg$gtre$dta$`FullDMXL+Liolaemus` <- output$gene.trees$`FullDMXL+Liolaemus`
inLwieg$gtre$dta$`FullDNAH3+Liolaemus` <- output$gene.trees$`FullDNAH3+Liolaemus`
inLwieg$gtre$dta$`FullEXPH5+Liolaemus` <- output$gene.trees$`FullEXPH5+Liolaemus`
inLwieg$gtre$dta$`FullKif24+Liolaemus` <- output$gene.trees$`FullKif24+Liolaemus`
inLwieg$gtre$dta$`FullMXRA5+Liolaemus` <- output$gene.trees$`FullMXRA5+Liolaemus`
inLwieg$gtre$dta$`FullND4+Liolaemus` <- output$gene.trees$`FullND4+Liolaemus`
inLwieg$gtre$dta$`FullPNN+Liolaemus` <- output$gene.trees$`FullPNN+Liolaemus`
inLwieg$gtre$dta$`FullPRLR+Liolaemus` <- output$gene.trees$`FullPRLR+Liolaemus`
inLwieg$gtre$dta$`FullSNCAIP+Liolaemus` <- output$gene.trees$`FullSNCAIP+Liolaemus`
inLwieg$gtre$dta$`Fullcytb+Liolaemus` <- output$gene.trees$`Fullcytb+Liolaemus`

inLwieg$pldy$dta = c(0.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.5, 2.0, 2.0, 2.0, 0.5)

#inLwieg$ptre$dta = pTrees	# si se usa COAL, hay que cargar el árbol con phybase usando función readtree.phybase

for (i in 1:length(output$species.trees)) {
inLwieg$stre$dta[[i]] <- output$species.trees[[i]]
}

return (inLwieg)

}

inLwieg <- sbpps_to_p2c2m (output)

################### P2C2M analyze ############################


.onAttach <-
function(libname = find.package("P2C2M"), pkgname = "P2C2M") {

  # Check if Richard Hudson's `ms` is present in compiled form
  # msdir = system.file("msdir", "", package="P2C2M")
  if (!file.exists(paste(msdir, "ms", sep="/"))) {
    #packageStartupMessage("-- Richard Hudson's 'ms' (Hudson 2002) is being compiled --")
    system(paste("cd", msdir, ";", "gcc -o ms ms.c streec.c rand1.c -lm"))
    #packageStartupMessage(" DONE -- \n\n")
  }

}
.onLoad <-
function(libname = find.package("P2C2M"), pkgname = "P2C2M") {
    #packageStartupMessage("-- Setting up a package-specific environment --")
    assign("P2C2M_globalVars", new.env(), envir=parent.env(environment()))
}

.onLoad()
.onAttach()

xml.file="beast.xml" 
descr.stats=c("NDC","LCWT")
beast.vers="1.8" 
single.allele=c("O")
num.reps=nreps 
use.sorted=FALSE
#use.mpi=FALSE
save.metadata=FALSE
verbose=TRUE
dbg=FALSE
                           
assign("P2C2M_flg_xmlFile", xml.file, envir = P2C2M_globalVars)
  # Note: not descr.stats, but descrStats
  assign("P2C2M_flg_dscrStat", descr.stats, envir = P2C2M_globalVars)
  assign("P2C2M_flg_beastV", beast.vers, envir = P2C2M_globalVars) 
  assign("P2C2M_flg_snglAl", single.allele, envir = P2C2M_globalVars)
  assign("P2C2M_flg_nReps", num.reps, envir = P2C2M_globalVars)
  assign("P2C2M_flg_srtBool", use.sorted, envir = P2C2M_globalVars)
  assign("P2C2M_flg_mpiBool", use.mpi, envir = P2C2M_globalVars)
  assign("P2C2M_flg_vrbBool", verbose, envir = P2C2M_globalVars)
  assign("P2C2M_flg_dbgBool", dbg, envir = P2C2M_globalVars)
  
  loghelpers.prntmngr <-
function(inText, uprFlg=F, nwlFlg=T) {
  ##############################
  # Function "loghelpers.prnt" #
  ##############################
  # Descr:  print text to screen
  # Deps:   -
  # I/p:    inText

  #metadataBool = get("P2C2M_flg_metadataBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  # Logging to screen
  if (verboseBool) {
      if (uprFlg & nwlFlg) {cat("\n", toupper(inText), "\n", sep="")}
      if (!uprFlg & nwlFlg) {cat("\n", inText, "\n", sep="")}
      if (!uprFlg & !nwlFlg) {cat(inText)}
  }

  ## Logging to parameter file
  #if (metadataBool) {
  #  if (tofileFlg) {
  #    writeLines(paste("\n## ", toupper(inText), "\n", sep=""),
  #               prmFHndl)
  #  }
  #}

}

corehelpers.metrics <-
function(gTree, pTree, pTreeNames, sTree, 
                               assoc, ploidy, descrStats, singleAllele) {
  # Descr:    corehelpers.metrics
  # Deps:     calc.lcwt
  #           calc.ndc
  #           calc.coal
  # I/p:      gTree
  #           pTree
  #           pTreeNames
  #           sTree
  #           assoc
  #           ploidy
  #           descrStats
  #           singleAllele

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> corehelpers.metrics", fg="red"), sep="")
  }

  outL = list()

  if ("GSI" %in% descrStats) {
    gsi = calc.gsi(gTree, assoc, singleAllele)
    outL$GSI = frmtMntse(gsi, 4)                                        # Setting a specific number of significands
  }

  if ("LCWT" %in% descrStats) {
    gtp = calc.lcwt(gTree, sTree, assoc, ploidy)
    outL$LCWT = frmtMntse(gtp, 4)
  }

  if ("NDC" %in% descrStats) {
    ndc = calc.ndc(gTree, sTree, assoc)
    outL$NDC = frmtMntse(ndc, 4)
  }

  if ("COAL" %in% descrStats) {
      
    # DEBUGLINES:
    #cat("\ngTree\n"); print(gTree)
    #cat("\npTree\n"); print(pTree)
    #cat("\npTreeNames\n"); print(pTreeNames)
    #cat("\nassoc\n"); print(assoc)
    
    ray = calc.coal(gTree, pTree, pTreeNames, assoc)
    outL$COAL = frmtMntse(ray, 4)
  }

  outD = unlist(outL)
  names(outD) = NULL
  outD = toString(outD)

  return(outD)
}

rtnlc <-
function (m, x) {
  # function for returning listcolumn x in a matrix of lists 
  # ("ReTurN List Column")
  # m = inMatrix
  aFun = function(r, c, m, x) {as.numeric(unlist(strsplit(m[r, c],",")))[x]}
  vect_aFun = Vectorize(aFun, vectorize.args = c('r','c'))
  return(outer(1:nrow(m), 1:ncol(m) , FUN = vect_aFun, m, x))
}

readhelpers.makeCFtable <-
function(ind) {
  # Descr:  generating a constituent-frequency table
  # Deps: -
  # I/p:  ind
  #
  # The command will perform the following conversion:
  # FROM:
  #         Var1            Freq
  #    [1,] "alpinus"        "4"
  # TO:
  #        spec              V2 
  #   Var1 "alpinus"         "1"
  #   Var1 "alpinus"         "2"
  #   Var1 "alpinus"         "3"
  #   Var1 "alpinus"         "4"

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readhelpers.makeCFtable", fg="red"), 
        sep="")
  }

  ind = as.matrix(as.data.frame(table(ind[,1])))
  aList = c()
  for (i in 1:length(ind[,1])) {
    aList = c(aList, rep(ind[i,1], times=ind[i,2]))
  }
  outd = as.matrix(as.data.frame(cbind(aList, 1:length(aList))))

  return(outd)
}

ms.exec <-
function(assoc, stree, nTips, ploidy, prmFHndl) {
  ######################
  # Function "ms.exec" #
  ######################
  # Descr:  conduction simulations via Hudson's ms   
  # Deps:   calchelpers.nodetips
  # I/p:    stree = the species tree, inferred via the multispecies 
  #                coalescent model
  #         nTips = the total number of tips to be simulated
  #         assoc = a dataframe with two columns and n rows
  #         prmFHndle
  #
  # Note:   "assoc": The first column contains the tip labels of the 
  #                  population tree and the second contains the number of 
  #                  alleles to be sampled from that population.
  #
  #         "ploidy": refers to *BEAST's "ploidy" in the xml files
  #                   and modifies the DMV values. When all loci have the same 
  #                   ploidy, it should be left as 1. When ploidy varies, it 
  #                   should be 0.5 for mitoch. and 2 for diploid nuclear.
  #
  #         "stree": all trees must be species trees with associated 
  #                 'dmv' values.
  #
  #         dmv-values: "stree$dmv" = Nu haploid alleles, or 2Nu individuals; 
  #                     hence, multiply by 2 to get the standard 4Nu 
  #                     population scaled mutation rate. BEAST's.dmvparse is 
  #                     equal to Nu, where N is the pop size in alleles. 
  #                     Hence: 2Nu in diploid individuals.

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  nReps = get("P2C2M_flg_nReps", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> ms.exec", fg="red"), sep="")
  }

###########################################
# 1. Preparatory actions, adjusting theta #
###########################################
  n_sp = length(stree$tip.label)

  theta = stree$dmv*2*ploidy
  # relTheta-values are theta-values devided by the root dmv (which is the 
  # last one in the list, hence: tail(x,1)
  relTheta = theta/tail(theta, 1)

  # Branch lengths of the species tree also must be divided by root theta
  stree$edge.length = stree$edge.length/tail(theta, 1)
  # Set up a data frame of branching times, where the row names (which equal
  # the node number) are their own column
  brTms = as.data.frame(sort(ape::branching.times(stree)))
  brTms[,2] = rownames(brTms)
  # nodeNames is a list of node names plus the root
  nodeNames = c(stree$edge[,2], (n_sp+1))

##############################
# 2. Converting associations #
##############################
  # FROM: alpinus,m206395a
  #       alpinus,m206395b
  #       alpinus,m206400a
  #       alpinus,m206400b
  # TO:   alpinus, 4
  assoc = as.matrix(as.data.frame(table(assoc[,1])))

###################################
# 3. Set initial pop sizes for ms #
###################################
  popSizes = ""
  for(i in 1:n_sp) {
    spName = which(stree$tip.label==assoc[i,1])
    # stree$dmv and stree$edge are in the same order
    size_subpop = relTheta[which(stree$edge[,2]==spName)]
   popSizes = paste(popSizes, "-n", i, size_subpop)
  }

##############################
# 4. Set island model for ms #
##############################
  # "length(stree$tip.label)" provides the number of tips in the species tree
  islModel = paste("-I", length(stree$tip.label), paste(assoc[,2], collapse=' '))

############################################################
# 5. Represent the species tree as past demographic events #
############################################################
  demogrEvents = ""
  for (r in 1:nrow(brTms)) {

    # "brTms" are the ages of all nodes, sorted from smallest to largest
    # "nodeNames" is a list of node names plus the root
    brTmName = brTms[r, 2]
    brTmVal = brTms[r, 1]

    # Calculate a scaling factor
    scaledBrTme = relTheta[which(nodeNames==brTmName)]

    # Take those branches that have the shortes length
    children = sort(stree$edge[stree$edge[,1] == brTmName, 2])
    child1 = sort(calchelpers.nodetips(stree, children[1]))[1]
    child2 = sort(calchelpers.nodetips(stree, children[2]))[1]
    tips = sort(c(child1, child2))

    # Check which of the elements in 
    # "stree$tip.label[tips[2]]" are also in "assoc[,1]"
    popI = which(assoc[,1]==stree$tip.label[tips[2]])
    popJ = which(assoc[,1]==stree$tip.label[tips[1]])

    # -ej t i j: Move all lineages in subpopulation i to subpopulation j at time t
    # -en t i x: Set subpop i size to x*N0 at time t and growth rate to zero
    demogrEvents = paste(demogrEvents, "-ej", brTmVal, popI, popJ, 
                                       "-en", brTmVal, popJ, scaledBrTme)
  }

###################################################
# 6. Combining specifications for full ms command #
###################################################
  # Path to Hudson's ms (Hudson 2002)
  # pathToMS = system.file("msdir", "ms", package="P2C2M")
  pathToMS = paste(msdir, "ms", sep="/")
  
  # "nTips" specifies the number of tips in the gene tree
  # "islModel" specifies the number of tips in the species tree and the composition of the number of tips in the gene tree
  # "popSizes" specifies the population sizes in ech subpopulation
  # "demogrEvents" specifies the species tree topology as past demographic events
  fullCmd = paste(pathToMS, nTips, nReps, islModel, popSizes,
                  demogrEvents, "-T")

#################
# 7. Execute ms #
#################
  # The grep command catches the tree spec via the trailing semicolon;
  # the first backslash necessary due to R, the second due to bash shell
  msOutp = system(paste(fullCmd, "| grep \\;"), intern=T)

##############################
# 8. Read in simulated trees #
##############################
  if (nReps == 1) {
    trees = list()
    trees[[1]] = ape::read.tree(text=msOutp)
  }
  else {
    trees = ape::read.tree(text=msOutp)
  }

  # Undoing the adjustment of species tree brlens 
  # as was necessary for ms-simulation
  if (length(trees) > 1) {
    for (t in 1:length(trees)) {
      trees[[t]]$edge.length = trees[[t]]$edge.length*tail(theta, 1)
    }
  }
  else {
    trees$edge.length = trees$edge.length*tail(theta, 1)
  }

  return(trees)
}

calchelpers.nodetips <-
function(tree, node) {
  # Descr:    returns the tip numbers for a given node
  # Deps:     calchelpers.nodetips
  # I/p:      tree = a phylog. tree with numbered internal nodes and numbered tips
  #           node = a node in that tree, specified as an integer

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.nodetips", fg="red"),
        sep="")
  }

    # DEBUGLINES:
    #cat("\ntree$edge\n"); print(tree$edge)
    #cat("\nnode\n"); print(node)

    n_tips = length(tree$tip.label)
    if (node <= n_tips) {node}
    #if (node <= n_tips) {                                              # Potential improvements than line above
    #    outD = node
    #    outD
    #}
    else 
        {
        outD = numeric()
        k = tree$edge[which(tree$edge[,1] == node),2]                   # CURRENT PROBLEM UNDER LCWT: When diff. allele numbers, node nicht in tree$edge[,1]
        for (j in k)
            {
            if (j <= n_tips) {outD = c(outD, j)}
            else {outD = c(outD, calchelpers.nodetips(tree, j))}
            }
        outD        
        }
}

Mode <-
function(x) {
  # function for calculating the mode
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

statshelpers.diffrnce <-
function (post_dist, post_pred_dist) {
  ####################################
  # Function "statshelpers.diffrnce" #
  ####################################
  # Descr:    calculates the difference between post_distirical and post_pred_distulated
  # Deps:     -
  # I/p:      post_dist
  #           post_pred_dist

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.diffrnce", 
        fg="red"), sep="")
  }

  # Note: The difference is "post_distirical-post_pred_distulated". Since "post_distirical" is worse 
  # than "post_pred_distulated" whenever it does not conform to the coalescent model (i.e.
  # has larger values), significant differences will be more positive than 
  # non-significant differences.
  diff = post_dist - post_pred_dist
  # TFL converts from type "list" to type "double"; is important, because
  # is.infinite and is.nan can only work on type "double"
  diff = as.matrix(diff)

  # Removing diff. values that are infinite ("Inf")
  diff = ifelse(is.infinite(diff), NA, diff)
  # Removing diff. values that are "NaN" ("not a number", i.e. 0/0)
  diff = ifelse(is.nan(diff), NA, diff)

  return(diff)
}

################ NDC CALCULATION ###################

calc.ndc <-
function(gTree, sTree, assoc) {
  # Descr:  returns the number of deep coalescences for an entire tree
  # Deps:   calc.parse
  #         calchelpers.gtreeparse
  # I/p:    sTree
  #         gTree
  #         assoc
  # Note:   ndc = "number of deep coalescences"

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.ndc",fg="red"),sep="")
  }

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD
  ndc = c()
  for (node in nodes) {
    tempData = calchelpers.gtreeparse(sTree, gTree, assoc, dmvD, node)
    ndc = c(ndc, tempData$lNd-1)
  }

  return(sum(ndc))
}

calc.parse <-
function(sTree, assoc) {
  # Descr:  parses species tree nodes for metric calculation
  # Deps:   calchelpers.dmvparse
  # I/p:    sTree
  #         assoc

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.parse",fg="red"),sep="")
  }

  # DEBUGLINES:
  #cat("\nassoc\n"); print(assoc)
  #cat("\nsTree$tip.label\n"); print(sTree$tip.label)
  
  spNames = sTree$tip.label
  n_sp = length(spNames)
  #nBr = (2*n_sp)-1
  tiplist = list()
  for(i in 1:n_sp) {
    tiplist[[spNames[i]]] = assoc[which(assoc[,1]==spNames[i]),2]
  }
  dmvD = calchelpers.dmvparse(sTree, n_sp)                              # returns demographic info of the species tree
  n_tips_per_sp = lapply(tiplist, length)                               # calculate number of tips per species

  if(any(n_tips_per_sp==1)) {                                           # evaluate number of tips per species tree
    tmp = dmvD[-which(n_tips_per_sp==1),]                               # Remove dmv values for terminals (i.e., where n_tips_per_sp==1)
  } else {
    tmp = dmvD
  }
  nodes = tmp[,"node"]
  names(nodes) = NULL

  outD = list()
  outD$nodes = nodes
  outD$dmvD = dmvD

  return(outD)
}

calchelpers.dmvparse <-
function(sTree, nSp) {
  # Descr:    parsing demographic data (dmv and dmt) of a branch
  # Deps:     -
  # I/p:      sTree
  #           nSp

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.dmvparse", fg="red"),
        sep="")
  }
  
    # DEBUGLINES:
    #cat("\nsTree\n"); print(sTree)
    #cat("\nsTree$edge\n"); print(sTree$edge)
    #cat("\nsTree$edge.length\n"); print(sTree$edge.length)
    #cat("\nsTree$dmv\n"); print(sTree$dmv)

    dmvD = cbind(sTree$edge[,2], sTree$edge.length)                     # generating rows of node number - branch length pairs
    dmvD = rbind(dmvD, c((nSp+1), Inf))                                 # adding another branch of length Inf
    dmvD = cbind(dmvD, sTree$dmv)                                       # adding dmv values
    dmvD = dmvD[order(dmvD[,1]),]                                       # order the matrix by the first column
    stBt = ape::branching.times(sTree)                                  # calc branching times of the species tree (via ape-function)

    # TFL may not be necessary, as stBt may already be sorted
    stBt = stBt[order(as.numeric(names(stBt)))]                         # sort the branching times by their names (which are numbers)

    pre = structure(rep(0, nSp), .Names=c(1:nSp))                       # add x zeros to the beginning of branching times list, where x = number of species in sTree
    stBt = c(pre, stBt)

    dmvD = cbind(dmvD, stBt)                                            # add the column "stBt" to the matrix "dmvD"
    colnames(dmvD) = c("node", "length", "dmv", "sbt")
    rownames(dmvD) = c(1:length(stBt))

return(dmvD)
}

calchelpers.gtreeparse <-
function(sTree, gTree, assoc, dmvD, node) {
# Descr:    parsing branches for likelihood calculations
# Deps:     calchelpers.descend
# I/p:      sTree
#           gTree
#           assoc
#           dmvD
#           node
# Note:     branching times = distance from each node to the tips, 
#           under the assumption that the tree is ultrametric
#           gBt = gene tree branching times (plural!)
#           fBt = first branching time
#           lBt = last branching time
#           fNd = first node
#           lNd = last node

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.gtreeparse", fg="red"),
        sep="")
  }

    # returns the subtree that starts at node <node>
    subtree = calchelpers.descend(sTree, gTree, assoc, node)
    # infers the branching times for a subtree (here the gene tree)
    gBt = sort(ape::branching.times(subtree))
    # generate a matrix with branching times in first column and node 
    # names in second column
    gBt = c(0, gBt)
    gBt = cbind(gBt, length(gBt):1)
    # get the branching times for a node in the species tree 
    fBt = dmvD[node,"sbt"]
    # get the branching times for a node in the species tree + the branch 
    # lengths for that node
    lBt = (dmvD[node,"sbt"] + dmvD[node,"length"])
    # get those node IDs (column 2), whose gene tree branching times (column 1) 
    # are equal the largest of the species tree branching times
    fBtMax = max(gBt[gBt[,1]<=fBt,1])
    fNd = gBt[gBt[,1]==fBtMax,2]
    # get those node IDs (column 2), whose  gene tree branching times (column 1) 
    # that equal the largest of the species tree branching times + branch lengths
    lBtMax = max(gBt[gBt[,1]<=lBt,1])
    lNd = gBt[gBt[,1]==lBtMax,2]

    outd = list("gBt"=gBt, "fBt"=fBt, "lBt"=lBt, "fNd"=fNd, "lNd"=lNd)

return(outd)
}

calchelpers.descend <-
function (sTree, gTree, assoc, cNode) {
##################################
# Function "calchelpers.descend" #
##################################
# Descr:    returns a tree containing descendants of all gene lineages 
#           that pass through a node; here this function is used to 
#           obtain the gene tree that is part of a species tree
# Deps:     calchelpers.nodetips
# I/p:      sTree
#           gTree
#           assoc
#           cNode = coalescence node

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.descend", fg="red"),
        sep="")
  }

    spNames = sTree$tip.label
    if (!is.na(cNode) && cNode != (length(spNames)+1))                  # Unless cNode is 'NA' AND unless the coalescence node is the root, do ...
        {
        subtreeTaxa = spNames[calchelpers.nodetips(sTree, cNode)]
        gTreeDscr = c()
        for(taxon in subtreeTaxa)
            {
            nodesRaw = assoc[which(assoc[,1] == taxon),2]
            if (is.character(nodesRaw[1]) == T)
                {nodes = nodesRaw}
            if (is.character(nodesRaw[1]) == F)
                {nodes = as.numeric(nodesRaw)}
            gTreeDscr = c(gTreeDscr, nodes)
            }
        taxaIndex = match(gTreeDscr, gTree$tip.label)
        # TFL removes all terminal branches of the gene tree 
        # that are not in the subtree; gTree must be in DNAbin-format
        subtree = ape::drop.tip(gTree, gTree$tip.label[-taxaIndex])
        }
    else {subtree = gTree}

return(subtree)
}

frmtMntse <-
function (inNum, mntse) {
  # function to control the significands of a number
return(as.numeric(format(round(inNum, mntse), nsmall=mntse)))
}


################ LCWT CALCULATION ###################

calc.lcwt <-
function(gTree, sTree, assoc, ploidy) {
  # Descr:    calculates the probability of a whole gene tree
  #           by calculating the probability of subtrees defined 
  #           by their MRCA (i.e. their node); done for all 
  #           nodes in a tree, values are summed up thereafter
  # Deps:     calc.parse
  #           calchelpers.brprob
  # I/p:      sTree
  #           gTree
  #           assoc
  #           ploidy
  # Note:     gtp = "gene tree probability" (a form of coalescent likelihood)

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.lcwt",fg="red"),sep="")
  }

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD
  
  # DEBUGLINES:
  #cat("\nnodes\n"); print(nodes)
  #cat("\ndmvD\n"); print(dmvD)

  lnP = c()
  for(node in nodes) {
    lnP = c(lnP, log(calchelpers.brprob(sTree, gTree, assoc, 
                                        ploidy, dmvD, node)))
  }

  return(sum(lnP))
}

calchelpers.brprob <-
function(sTree, gTree, assoc, ploidy, dmvD, node) {
# Descr:    organising likelihood calculations for particular branches
#           Calculate prob of gene tree part, which is a subset 
#           of the species tree.
#           More general: Calculate the prob of subtrees defined 
#           by their MRCA (i.e. their node). 
# Deps:     calchelpers.gtreeparse
#           calchelpers.probcalc
# I/p:      sTree
#           gTree
#           assoc
#           ploidy
#           dmvD = demographic data
#           node
# Note:     branching times = distance from each node to the tips, 
#            under the assumption that the tree is ultrametric

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.brprob", fg="red"),
        sep="")
  }

    ind = calchelpers.gtreeparse(sTree, gTree, assoc, dmvD, node)

    gBt = ind$gBt
    fBt = ind$fBt
    lBt = ind$lBt
    fNd = ind$fNd
    lNd = ind$lNds

    if(fBt!=0) {prob=1}
    else
        {
        # Retain those gTree branching times (i.e. column 1) whose node names 
        # (i.e. column 2) fall within "fBt" and "lBt"
        gBt = gBt[gBt[,1]>=fBt & gBt[,1]<=lBt,]
        # Add lNd row at bottom of gBt-matrix
        gBt = rbind(gBt, c(lBt, lNd))
        # Perform actual prob calculation
        prob = calchelpers.probcalc(gBt, dmvD, ploidy, node, fBt, lBt, fNd, lNd)
        }

return(prob)
}

calchelpers.probcalc <-
function(gBt, dmvD, ploidy, node, fBt, lBt, fNd, lNd) {
# Descr:    calculating likelihoods for particular branches
# Deps:     -
# I/p:      gBt
#           dmvD
#           ploidy
#           node
#           fBt
#           lBt
#           fNd
#           lNd
# Note:     The input variable "gBt" contains two columns: 
#           branching time differences and node ID
#           branching times = distance from each node to the tips, 
#           under the assumption that the tree is ultrametric

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.probcalc", fg="red"),
        sep="")
  }

# 1. Calculating differences in branching times
    value = length(gBt[,1])
    # waitTms = branching time differences between two branches
    waitTms = gBt[2:value,1] - gBt[1:value-1,1]
    # append node IDs again (so that "waitTms" mimicks "gBt"), 
    # with the exception of last node ID in list
    waitTms = cbind(waitTms, gBt[1:value-1,2])
# 2. Calculating lists of values
    dmv = (2*dmvD[node,"dmv"])*ploidy
    # lambda is a list of values, each calculated according to the following formula
    lambda = (waitTms[,2]*(waitTms[,2]-1))/dmv
    # exponent is a list of values; Euler-constant to the power of (-lambda*waitTms[,1])
    exponent = exp(-lambda*waitTms[,1])
# 3. Calculating products for each list
    # Calculate product of exponent list, given that lambda is never 0
    exponent = prod(exponent[lambda!=0])
    # Calculate product of lambda list, however disregarding last list element
    lambda = prod(lambda[1:(length(lambda)-1)])
# 4. Calculating overall probability
    prob = lambda * exponent

return(prob)
}

################ P2C2M ANALYZE ###################

p2c2m.analyze <-
function (inFn, inData, prmFHndl) {
  # Descr:  core of script, coordinates the metric calculation
  # Deps:   (various)
  # I/p:    inFn
  #         inData
  #         prmFHndl

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  descrStats = get("P2C2M_flg_dscrStat", envir=P2C2M_globalVars)
  mpiBool = get("P2C2M_flg_mpiBool", envir=P2C2M_globalVars)
  nReps = get("P2C2M_flg_nReps", envir=P2C2M_globalVars)
  singleAllele = get("P2C2M_flg_snglAl", envir=P2C2M_globalVars)
  srtBool = get("P2C2M_flg_srtBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> core.analyze", fg="red"), sep="")
  }

##########################
# 1. Preparatory actions #
##########################
  n_loci = length(inData$loci$dta)
  n_sTrees = length(inData$stre$dta)
  # generate empty matrices of dimensions 'matrix'
  post_dist = post_pred_dist = matrix(nrow=n_sTrees, ncol=n_loci)
  # label the columns of the empty matrices by the locus names
  colnames(post_dist) = colnames(post_pred_dist) = inData$loci$dta
  # Initilizing outdata 
  outD = list()

##########################################################
# 2. Metrics for the trees of the posterior distribution #
##########################################################
  loghelpers.prntmngr(paste("Metrics for the trees of the", 
                            "posterior distribution"), uprFlg=T)

 ##########################################################
 # 2.a. Looping through genes, calc. descript. statistics #
 ##########################################################
  # Loop through all gene tree names
  for (j in 1:n_loci) {
    loghelpers.prntmngr(paste("Locus '", inData$loci$dta[j], "'",
                              sep=""))
    loghelpers.prntmngr("Analyzing gene trees:\n", nwlFlg=F)
    # Loop through all species trees
    for (i in 1:n_sTrees) {
      # Display tree number on screen
      loghelpers.prntmngr(paste(" ", i, sep=""), nwlFlg=F)
      # Save tree number as replicate ID
      assign("P2C2M_flg_repID", paste("post_dist", i, sep="."), 
             envir = P2C2M_globalVars)

    # DEBUGLINES:
    #cat("\ninData$gtre$dta[[j]][[i]]\n"); print(inData$gtre$dta[[j]][[i]])
    #cat("\ninData$ptre$dta$tree[[i]]\n"); print(inData$ptre$dta$tree[[i]])
    #cat("\ninData$ptre$dta$names\n"); print(inData$ptre$dta$names)
    #cat("\ninData$stre$dta[[i]]\n"); print(inData$stre$dta[[i]])
    #cat("\ninData$asoc$dta\n"); print(inData$asoc$dta)
    #cat("\ninData$pldy$dta[j]\n"); print(inData$pldy$dta[j])

      post_dist[i,j] = corehelpers.metrics(inData$gtre$dta[[j]][[i]],
                                     inData$ptre$dta$tree[[i]],
                                     inData$ptre$dta$names,
                                     inData$stre$dta[[i]],
                                     inData$asoc$dta[[j]],
                                     inData$pldy$dta[j],
                                     descrStats,
                                     singleAllele)
    }
    loghelpers.prntmngr("\n", nwlFlg=F)
  }

 ###########################
 # 2.b. Assembling results #
 ###########################
  for (s in descrStats) {
    # Populating list "outD" with the results of the trees of the posterior distribution
    outD[[s]]$post_dist$unsort = as.data.frame(rtnlc(post_dist, which(s==descrStats)))
    # Assigning colnames
    colnames(outD[[s]]$post_dist$unsort) = inData$loci$dta
    # Generating a sorted result clone
    outD[[s]]$post_dist$sorted = apply(outD[[s]]$post_dist$unsort, MARGIN=2, sort)
  }

#####################################################################
# 3. Metrics for the trees of the posterior predictive distribution #
#####################################################################
  loghelpers.prntmngr(paste("Metrics for the trees of the", 
                            "post. predictive distribution"), uprFlg=T)

 ##################################
 # 3.a. Initiate MPI (if applic.) #
 ##################################
  # Determining number of processors (i.e. CPUs) available
  if (Sys.info()['sysname'] == "Linux") {
    nproc = as.numeric(system("nproc", intern=T))
  }
  # MacOS option
  if (Sys.info()['sysname'] == "Darwin") {
    nproc = as.numeric(system("sysctl -n hw.ncpu", intern=T))
  }

  if (mpiBool && nproc >= 3) {
    loghelpers.prntmngr(paste("\n", "Initialize parallelization", "|", 
                              "N of processors:", toString(nproc)))
    # Specifying one master and nproc-1 slaves
    #Rmpi:: mpi.spawn.Rslaves(nslaves=nproc-1, needlog=F)
    Rmpi:: mpi.spawn.Rslaves(nslaves=nproc-1)
    # Passing separate functions to MPI environment
    mpiEnvir = c(calc.lcwt,
               #  calc.gsi,
                 calc.ndc,
                 calc.parse,
               #  calc.coal,
                 calchelpers.brprob,
                 calchelpers.descend,
                 calchelpers.dmvparse,
                 calchelpers.gtreeparse,
                 calchelpers.nodetips,
                 calchelpers.probcalc,
                 corehelpers.metrics,
                 # LEGACY: corehelpers.repl, 
                 frmtMntse,
                 # LEGACY: apTreeshape::as.treeshape.phylo,
                 ape::branching.times,
                 ape::drop.tip)
             #    genealogicalSorting::gsi,
             #    phybase::loglikeSP)
    lapply(mpiEnvir, Rmpi::mpi.bcast.Robj2slave)
  }

 ##############################
 # 3.b. Looping through genes #
 ##############################
  simTrees_across_loci = list()
  for (j in 1:n_loci) {
    loghelpers.prntmngr(paste("Locus '", inData$loci$dta[j], "'",
                              sep=""))
    loghelpers.prntmngr("Analyzing the posterior predictive trees:\n",
                        nwlFlg=F)
    n_tips = length(inData$gtre$dta[[j]][[1]]$tip.label)
    simTrees_across_sTrees = list()
    for (i in 1:n_sTrees) {
      # Display tree number on screen
      loghelpers.prntmngr(paste(" ", i, sep=""), nwlFlg=F)
      # Save tree number as replicate ID
      assign("P2C2M_flg_repID", paste("post_pred_dist", j, i, sep="."), 
             envir = P2C2M_globalVars)

 ############################
 # 3.c. Simulation of trees #
 ############################
      # print(inData$assoc)
      # speciesName   alleleName
      # [1,] "A"      "a1"
      # [2,] "A"      "a10"
      # [3,] "A"      "a11"
      # [4,] "A"      "a12"
      # [5,] "A"      "a2"

      # A constituent-frequency table is required for the simulations in ms,
      # because ms saves the alleles as numbers such as "1", "2", "3", ...
      CFTable = readhelpers.makeCFtable(inData$asoc$dta[[j]])

      # print(CFTable)
      #       aList  V2
      # Var1  "A"    "1"
      # Var1  "A"    "2"
      # Var1  "A"    "3"

      # Simulate gene trees
      simTrees = ms.exec(CFTable, inData$stre$dta[[i]], n_tips, 
                         inData$pldy$dta[j], prmFHndl)
      class(simTrees) = "multiPhylo"

      ## Possible improvement in TFL: why not "length(simTrees)" instead of 
      ## "as.numeric(nReps)"
      simReps = matrix(nrow=as.numeric(nReps), ncol=1)

 ##################################################################
 # 3.d. Calculation of descript. statistics depend. on MPI status #
 ##################################################################
      if (mpiBool && nproc >= 3) {
        # Applying function "corehelpers.metrics" in parallel
        #  The second element in mpi.parLapply has to be the function that
        #  the data is to be applied
        simReps = Rmpi::mpi.parLapply(simTrees,
                                      corehelpers.metrics,
                                      inData$ptre$dta$tree[[i]],
                                      inData$ptre$dta$names,
                                      inData$stre$dta[[i]],
                                      CFTable,
                                      inData$pldy$dta[j],
                                      descrStats,
                                      singleAllele)
        simReps = as.matrix(as.character(simReps))
      }

      else {
        for (k in 1:nReps) {
        simReps[k,1] = corehelpers.metrics(simTrees[[k]],
                                           inData$ptre$dta$tree[[i]],
                                           inData$ptre$dta$names,
                                           inData$stre$dta[[i]],
                                           CFTable,
                                           inData$pldy$dta[j],
                                           descrStats,
                                           singleAllele)
        }
      }

 #########################################
 # 3.e. Parsing of average metric values #
 #########################################
      # Extracting the listcolumns of the matrix simReps via GSO.rtnlc
      sumD = list()
      for (s in descrStats) {
        sumD[[s]] = Mode(rtnlc(simReps, which(s==descrStats)))
      }
      sumD = unlist(sumD)
      names(sumD) = NULL
      
      # Save the replicate values to variable "post_pred_dist"
      post_pred_dist[i,j] = toString(sumD)

 #############################################
 # 3.f. Preparing simulated trees for saving #
 #############################################

      # Loop over simulated replicates
      ## Possible improvement: Use apply function instead of loop
      # Generate set of print-ready simtrees
      simTrees_wLabels = list()
      for (k in 1:nReps) {
        # LEGA!CY: # Generating a simTree version where tips have been replaced with the correct allele names
        # LEGA!CY: assoc_ids = unlist(lapply(inData$asoc$dta[[j]][,2], substring, first=2))
        # LEGA!CY: # TFL removes leading zeros
        # LEGA!CY: assoc_ids = gsub("^[0]+", "", assoc_ids)
        # LEGA!CY: # TFL removes all non-numeric characters
        # LEGA!CY: oldLbls = gsub("[^0-9]", "", simTrees[[k]]$tip.label)

        # LEGA!CY: # DEBUGLINES: 
        # LEGA!CY: #cat("\nassoc_ids\n"); print(assoc_ids)
        # LEGA!CY: #cat("\nsimTrees[[k]]$tip.label\n"); print(simTrees[[k]]$tip.label)
        # LEGA!CY: #cat("\noldLbls\n"); print(oldLbls)
        
        # LEGA!CY: match_indcs = match(oldLbls, assoc_ids)
        # LEGA!CY: # Error handling
        # LEGA!CY: if (NA %in% match_indcs | "NA" %in% match_indcs) {
        # LEGA!CY:   stop(cat("\nERROR: Cannot label the terminals of the posterior predictive trees correctly.\n"))
        # LEGA!CY: }
        # LEGA!CY: nwLbls = inData$asoc$dta[[j]][,2][match_indcs]
        # LEGA!CY: simTrees_wLabels[[k]] = simTrees[[k]]
        # LEGA!CY: simTrees_wLabels[[k]]$tip.label = mapply(toString, nwLbls)
        simLbls = gsub("[^0-9]", "", simTrees[[k]]$tip.label)
        nwLbls = inData$asoc$dta[[j]][as.integer(simLbls),2]
        simTrees_wLabels[[k]] = simTrees[[k]]
        simTrees_wLabels[[k]]$tip.label = mapply(toString, nwLbls)       
      class(simTrees_wLabels) = "multiPhylo"
      }
      simTrees_across_sTrees[[paste("sTree_", i, sep="")]] = simTrees_wLabels


    # End of looping through sTrees
    }
    simTrees_across_loci[[inData$loci$dta[[j]]]] = simTrees_across_sTrees
    loghelpers.prntmngr("\n", nwlFlg=F)

  # End of looping through loci
  }

 ##################
 # 3.g. Close MPI #
 ##################
  if (mpiBool && nproc >= 3) {
    # Closing slaves
    loghelpers.prntmngr(paste("\n", "Close parallelization", "\n",
                              sep=""))
    Rmpi::mpi.close.Rslaves(dellog=T)
  }
  
 #############################
 # 3.g. Save simulated trees #
 #############################
  outD$simTrees = simTrees_across_loci 

 ###########################
 # 3.h. Assembling results #
 ###########################
  for (s in descrStats) {
    # Populating list "outD" with the simulated results
    # Both "as.matrix" and "as.data.frame" are critical below
    outD[[s]]$post_pred_dist$unsort = as.matrix(as.data.frame(rtnlc(post_pred_dist, which(s==descrStats))))
    # Assigning colnames
    colnames(outD[[s]]$post_pred_dist$unsort) = inData$loci$dta
    # Generating a sorted result clone
    outD[[s]]$post_pred_dist$sorted = apply(outD[[s]]$post_pred_dist$unsort, MARGIN=2, sort)
  }

#####################################
# 4. Calculating metric differences #
#####################################
  loghelpers.prntmngr(paste("Calculating differences btw. the",
                            "posterior and the posterior predictive", 
                            "distribution"))
  
  for (s in descrStats) {
    # Populating list "outD" with the differences
    if (srtBool) {
      outD[[s]]$dif = statshelpers.diffrnce(outD[[s]]$post_dist$sorted, 
                                            outD[[s]]$post_pred_dist$sorted)
      # Remove the outD set not chosen
      outD[[s]]$post_dist$unsort = NULL
      outD[[s]]$post_pred_dist$unsort = NULL
    }
    else {
      outD[[s]]$dif = statshelpers.diffrnce(outD[[s]]$post_dist$unsort, 
                                            outD[[s]]$post_pred_dist$unsort)
      # Remove the outD set not chosen
      outD[[s]]$post_dist$sorted = NULL
      outD[[s]]$post_pred_dist$sorted = NULL
    }
  }

  return(outD)
}

# resultsP2C2M$inData is the format that P2C2M parses Beast data using the Python scripts
# make inData format of P2C2M from the inData format of starbeastPPS2


statistics<-p2c2m.analyze("/Users/arley/Desktop/pruebaSBPPS/beast1/beast1_8.xml", inLwieg)


#################### STATS MAIN #######################

stats.main <-
function (path, xml.file, loci, resultData, prmFile) {
  # Descr:    coordinates executing of modules
  # Deps:     (various)
  # I/p:      path = absolute path to Working directory
  #           xml.file = name of infile
  #           loci = names of loci
  #           resultData
  #           prmFile

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.main", fg="red"), "\n",
        sep="")
  }

##########################
# 1. Summarizing results #
##########################
  loghelpers.prntmngr("Summarizing results", uprFlg=T)
                                               
  # Setting outdata
  outData = list()
  # Setting alpha values
  alphaValues = c(0.1, 0.05, 0.01)
  for (val in alphaValues) {
    # Coordinating the calculation of the statistics
    results = stats.coord(resultData, loci, val)
    valStr = paste("alpha", as.character(val), sep="")
    outData[[valStr]] = results
  }

#####################
# 2. Writing legend #
#####################
  legend = "Differences between the posterior and the posterior predictive distributions. Each cell contains the following information in said order: mean, standard deviation, significance level. Codes in square brackets indicate the number of tails. Alpha values are automatically adjusted for the number of tails."
  outData$legend = legend
return(outData)
}

stats.coord <-
function(ind, loci, alpha) {
  # Descr:  generate significance tables
  # Deps:   stats.outmaker
  # I/p:    ind
  #         loci
  #         alpha
  # Note:   CV = coefficient of variance
  #         avgs = arithmetic means
  #         qntls = quantiles

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  slctn = get("P2C2M_flg_dscrStat", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.coord", fg="red"), "\n",
        sep="")
  }

##############################
# 1. Setting number of tails #
##############################

  # T-tests for all descriptive statistics are two-tailed, because there is no
  # a priori reason in which direction they should differ.
  tailL = list()
  tailL$LCWT = "2"
  tailL$COAL = "2"
  tailL$NDC = "2"
  tailL$GSI = "2"
  ## T-tests for NDC should be left one-tailed, because trees not compliant 
  ## with the coalescent model have a higher number of deep coalescences.
  #NDCtail = "1l"
  ## T-tests for GSI should be right one-tailed, because trees not compliant 
  ## with the coalescent model have lower values.
  #GSItail = "1r"

#############################
# 2. Inferring significance #
#############################
  perGene = acrGenes = list()
  for (s in slctn) {
    perGene[[s]] = stats.perGene(ind[[s]]$dif, alpha, tailL[[s]])
    acrGenes[[s]] = stats.acrGenes(ind[[s]]$dif, alpha, tailL[[s]])
  }

#######################
# 3. Combining output #
#######################
  outList = list()
  # perGene output
  outList$perGene = sapply(perGene, cbind)
  rownames(outList$perGene) = c(loci)
  # acrossGene output
  outList$acrGenes = sapply(acrGenes, cbind)
  rownames(outList$acrGenes) = c("Sum", "Mean", "Median", "Mode", "CV")
  # Naming rows of output
  names = c()
  for (stat in slctn) {
    names = c(names, paste(stat, "[", tailL[[stat]], "]", sep=""))
  }
  colnames(outList$perGene) = colnames(outList$acrGenes) = names
  
  return(outList)
}

stats.perGene <-
function (diff, alpha, tail) {
  # Descr:    calculates statistics per gene
  # Deps:     statshelpers.qntls
  # I/p:      diff
  #           alpha
  #           tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.perGene", fg="red"), 
        sep="")
  }

  # Structure of "diff" at this point:
  #        [gene1]   [gene2]  [gene3]
  # [1,]  -2.72395  -9.48870 45.79565
  # [2,]  14.97560  38.06380 88.60285
  # [3,]   8.76280  21.09895 50.51950
  qntls = statshelpers.qntls(diff, alpha, tail)
  sigSgns = statshelpers.sigsgn(qntls, tail)
  # Calculation of means
  # Note: "MARGIN=2" means "across a column"
  subCol1 = c(apply(diff, MARGIN=2, mean))
  # Calculation of stdv
  subCol2 = c(apply(diff, MARGIN=2, sd))
  # Assignment of signif. values
  subCol3 = c(sigSgns)
  # "\u00B1" is unicode sign of "plusminus"
  # Output format:   mean(+-1 SD)[:space:]signif.level
  #                  Mean and SD are rounded to two decimal places.
  outData = paste(round(subCol1, 2), 
               paste("(", "\u00B1", round(subCol2, 2), ")", sep=""), 
               subCol3)

  return(outData)
}

statshelpers.qntls <-
function (ind, alpha, tail) {
  # Descr:  generates a distribution of quantiles;
  #         quantile levels are hardcoded in variable "qntlLevels"
  # Deps:   (none)
  #         ind = a data frame
  #         alpha
  #         tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.qntls",
        fg="red"), sep="")
  }

  #################################
  # Set quantile levels via alpha #
  #################################
  if (tail=="1l" | tail=="1r") {
    qntlLevels = c(alpha, 1-alpha)
  }
  # Adjusting alpha value for two-tailed test
  if (tail=="2") {
    qntlLevels = c((alpha/2), 1-(alpha/2))
  }

  qntls = t(apply(ind, MARGIN=2, quantile, qntlLevels, na.rm=T))
  #qntls = rbind(qntls, quantile(acrossGenes, qntlLevels, na.rm=T))

  return(qntls)
}

statshelpers.sigsgn <-
function (qntls, tail) {
  ##################################
  # Function "statshelpers.sigsgn" #
  ##################################
  # Descr:    applies significance signs
  # Deps:     -
  # I/p:      qntls = a data frame of two columns

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.sigsgn", 
        fg="red"), sep="")
  }

  sigSgns = rep(0, length(qntls[,1]))

  # Schemes for one-tailed test
  if (tail=="1l") {
    # Read: "Whichever elements of quants are smaller than zero,
    #        are considered significant and, hence, receive a star 
    #        in the first column."
    sigSgns[which(qntls[,1] > 0)] = "*"
    sigSgns[which(qntls[,1] < 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,1] == 0)] = "n.s."
  }
  if (tail=="1r") {
    sigSgns[which(qntls[,2] < 0)] = "*"
    sigSgns[which(qntls[,2] > 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,2] == 0)] = "n.s."
  }

  # Scheme for two-tailed test
  if (tail=="2") {
    sigSgns[which(qntls[,1] > 0 | qntls[,2] < 0)] = "*"
    sigSgns[which(qntls[,1] < 0 & qntls[,2] > 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,1] == 0 & qntls[,2] > 0)] = "n.s."
    sigSgns[which(qntls[,1] < 0 & qntls[,2] == 0)] = "n.s."
  }

  return(sigSgns)
}

stats.acrGenes <-
function (diff, alpha, tail) {
  # Descr:    calculates statistics across genes
  # Deps:     statshelpers.qntls
  #           statshelpers.cv
  # I/p:      diff
  #           alpha
  #           tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.acrGenes", fg="red"),
        sep="")
  }

  # Perform calculations per row (i.e per MCMC generation)
  acrossG_Sum = rowSums(diff)
  acrossG_Mean = rowMeans(diff)
  acrossG_Median = rowMedians(diff)
  acrossG_Mode = rowModes(diff)
  acrossG_CV = statshelpers.cv(diff)

  # "acrossG" is a matrix with four columns hereafter
  acrossG = cbind(acrossG_Sum, 
                  acrossG_Mean, 
                  acrossG_Median,
                  acrossG_Mode,
                  acrossG_CV)

  qntls = statshelpers.qntls(acrossG, alpha, tail)
  sigSgns = statshelpers.sigsgn(qntls, tail)

  subCol1 = c(apply(acrossG, MARGIN=2, mean))
  subCol2 = c(apply(acrossG, MARGIN=2, sd))
  subCol3 = c(sigSgns)
  outData = paste(round(subCol1, 2), 
               paste("(", "\u00B1", round(subCol2, 2), ")", sep=""), 
               subCol3)

  return(outData)
}

rowMedians <-
function(x) {
  # function for calculating the median of all row values
  # analogous to rowSums and rowMeans
  apply(x, MARGIN=2, median, na.rm=T)
}

rowModes <-
function(x) {
  # function for calculating the mode of all row values
  # analogous to rowMeans and rowMedians
  apply(x, MARGIN=2, Mode)
}

statshelpers.cv <-
function(diff) {
  # Descr:    Calculates the coefficient of variance
  # Deps:     (none)
  # I/p:      inD
  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.cv", fg="red"), sep="")
  }

  # Structure of "diff" at this point:
  #        [gene1]   [gene2]  [gene3]
  # [1,]  -2.72395  -9.48870 45.79565
  # [2,]  14.97560  38.06380 88.60285
  # [3,]   8.76280  21.09895 50.51950

#  inD = as.data.frame(diff)
#  inD = ifelse(is.nan(diff), NA, inD)
  Stdv = apply(diff, MARGIN=1, sd, na.rm=T)
  Mean = rowMeans(diff)
  cv = Stdv/Mean

  # TFL is CRITICAL, because there are occasional "Inf" in the matrix for 
  # descriptive statistics "NDC"
  is.na(cv) <- do.call(cbind, lapply(cv, is.infinite))
  return(cv)
}

mainstats <- stats.main("/Users/arley/Desktop/pruebaSBPPS/beast1", "beast1_8.xml", inLwieg$loci$dta, statistics)

results <- list()
results[["stats"]] <- mainstats
results[["NDC"]] <- statistics$NDC
results[["LCWT"]] <- statistics$LCWT
return (results)
    
}