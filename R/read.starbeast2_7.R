read.starbeast2_7 <- function(inFn, apWd) {
    
    ##### edit xml before running starbeast2 in order to write branch rate info in gene tree files
    ##### add the following: substitutions="true" branchratemodel="@StrictClock.c:26" (specify branchRateModel id in the dots)
    ##### read in a beast xml file in order to extract the taxon-allele bindings and
    ##### locus names.

# source("read.gtree.samples2.R")

# source("read.sptree.samples2.R")

xmlfile2 <- paste(unlist(c(apWd, "/", inFn)), collapse="")

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

loci2 <- D2[grep("<distribution id=\"geneTree", D2)]

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

orderploidy  <- gsub(".*<distribution.*geneTree.t:", "", loci2)
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


neworder<-gsub("Tree.t:","",gsub(" spec.*","",gsub("\"","", gsub(".*id=\"","",D2[grep("tree id=", D2)]))))
rownames(genes.to.partitions2)<-genes.to.partitions2[,1]
genes.to.partitions2 <- genes.to.partitions2[neworder,]
#ploidy2 <- ploidy2[neworder]
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

cat("\nbindings, loci, and alignments read\n\n")
    
    ### read in gene trees from combinedfiledirectory. assumes custom perl script used
    ### to thin and burn in beast results.

# library(ape)
    
    gene.trees2 <- list()
    for (i in 1:length(loci2)) {
         treefile2 <- c(apWd, "/", loci2[i], ".trees")
         gene.trees2[[loci2[i]]] <- read.gtree.samples2(paste(unlist(treefile2), collapse = ""))
         cat("locus", loci2[i], "is done\n\n")
    }
    cat("gene trees done\n\n")   
 
    ### read in species trees from combinedfiledirectory. assumes custom perl script
    ### used to thin and burn in beast results.
    
    species.trees2 <- read.sptree.samples2(paste(unlist(c(apWd, "/", 
        "species.trees")), collapse = ""))
        cat("species trees done\n\n")
 

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
    
 #   print(empirical)
    
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
  #  class(output) <- "starbeast.data"
  #  return(output)


#################### OUTPUT TO P2C2M ###########################

inData = list()

inData$asoc$nme = "ordered allele-species associations"
inData$loci$nme = "loci"
inData$gtre$nme = "gene trees"
inData$pldy$nme = "ploidy level"
inData$ptre$nme = "phybase species trees"
inData$stre$nme = "regular species trees"


for (i in 1:length(loci2)) {
	inData$asoc$dta[[i]]<-output$associations$empirical[[i]]
	colnames(inData$asoc$dta[[i]])<-c("speciesName","alleleName")
}

inData$loci$dta<-as.list(names(output$genes))

for (i in 1:length(loci2)) {
	inData$gtre$dta[[i]] <- output$gene.trees[[i]]
}

inData$pldy$dta = output$ploidy

# if using COAL, load tree with phybase using the function readtree.phybase
# inData$ptre$dta = pTrees	

for (i in 1:length(output$species.trees)) {
inData$stre$dta[[i]] <- output$species.trees[[i]]
}

return(inData)

}

