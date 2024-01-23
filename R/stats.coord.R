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
  errBool = get("P2C2M_flg_errBool", envir=P2C2M_globalVars)
  errReps = get("P2C2M_flg_errReps", envir=P2C2M_globalVars)
  
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
  perGene.error = acrGenes.error = list()
  perGene.ersig = acrGenes.ersig = list()
  perGene.ersig.mat = acrGenes.ersig.mat = list()
  
  for (s in slctn) {
    perGene[[s]] = stats.perGene(ind[[s]]$dif, alpha, tailL[[s]])
    acrGenes[[s]] = stats.acrGenes(ind[[s]]$dif, alpha, tailL[[s]])

    if (errBool) {
      for (r in 1:errReps) {
        perGene.error[[s]][[r]] = stats.perGene(ind[[s]]$err[[r]], alpha, tailL[[s]])
        perGene.ersig[[s]][[r]] <- grepl("n.s.",perGene.error[[s]][[r]])
        acrGenes.error[[s]][[r]] = stats.acrGenes(ind[[s]]$err[[r]], alpha, tailL[[s]])
        acrGenes.ersig[[s]][[r]] <- grepl("n.s.",acrGenes.error[[s]][[r]])
      }       
      perGene.ersig.mat[[s]] = round((errReps-rowSums(sapply(perGene.ersig[[s]], rbind)*1))/errReps, digits=3)
      acrGenes.ersig.mat[[s]] = round((errReps-rowSums(sapply(acrGenes.ersig[[s]], rbind)*1))/errReps, digits=3)
    }
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
  
  if (errBool) {
    outList$perGene.error = sapply(perGene.ersig.mat, rbind)
    outList$acrGenes.error = sapply(acrGenes.ersig.mat, rbind)
    rownames(outList$perGene.error) = c(loci)
    rownames(outList$acrGenes.error) = c("Sum", "Mean", "Median", "Mode", "CV")
    names = c()
    for (stat in slctn) {
      names = c(names, paste(stat, "[", tailL[[stat]], "]", sep=""))
    }
    colnames(outList$perGene.error) = colnames(outList$acrGenes.error) = names
 }

  return(outList)
}
