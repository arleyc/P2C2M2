\name{p2c2m2.complete}
\alias{p2c2m2.complete}
\title{Execute the complete P2C2M2 pipeline via a single command}

\description{
  This function executes the complete \pkg{P2C2M2} pipeline from 
  beginning to end.
}

\usage{
p2c2m2.complete(path = "/home/user/Desktop/", xml.file = "beast.xml", 
    descr.stats = "LCWT,NDC", beast.vers = "2.7", 
    single.allele = c("O"), num.reps = 1, error.rate = FALSE, 
    error.reps = 1000, use.sorted = FALSE, use.mpi = FALSE, 
    save.metadata = FALSE, verbose = FALSE, dbg = FALSE)
}

\arguments{
  \item{path}{the absolute file path to the input directory, specified 
    as a double-quoted string; if \code{"/home/user/Desktop/"} 
    (the default), then the desktop itself is considered the input 
    directory.}

  \item{xml.file}{the name of the BEAUTi-generated and XML-formatted 
    input file, specified as a double-quoted string. The default is 
    \code{"beast.xml"}.}

  \item{descr.stats}{the name(s) of the summary statistic(s) to be 
    applied, specified as a double-quoted string. If multiple 
    statistics are specified, they must be separated by commas. A total 
    of four summary statistics is currently available: 
    \code{"COAL"} and \code{"LCWT"} (both Rannala & Yang 2013), 
    \code{"GSI"} (Cummings et al. 2008), \code{"NDC"} (Maddison 1997). 
    The default is \code{"LCWT,NDC"}.}

  \item{beast.vers}{the version of starBEAST2 (...) 
    used to perform the species tree inference, specified as a 
    double-quoted string. Data parsers are located in the subdirectory 
    exec/. Currently, the following parser is available: \code{"2.7"}.}

  \item{single.allele}{the name of a species that is represented by 
    only a single allele, specified as a variable of mode vector. This 
    setting is useful when defining an outgroup, because the species so 
    defined does not contribute towards the calculation of the summary 
    statistic 'GSI'. The default is \code{c("O")}.}

  \item{num.reps}{the number of simulation replicates to be conducted, 
    specified as an integer. The default is \code{1} (i.e., no
    replication).}

  \item{error.rate}{a logical specifying if error rates are estimated
    by sampling from the posterior predictive distributions. The default     is \code{FALSE}.}
  
  \item{error.reps}{the number of pseudo-observed datasets or "pods" to
    be sampled from the posterior predictive distributions, specified as 
    an integer. The default is \code{1000}.}

  \item{use.sorted}{a logical specifying if the summary statistics 
    generated from the posterior and from the posterior predictive 
    distribution are to be ranked by magnitude prior to the calculation 
    of the differences and the formation of the test distribution. 
    The default is \code{FALSE}. This argument is only EXPERIMENTAL and 
    should not be selected by regular users.}

  \item{use.mpi}{a logical specifying if \pkg{P2C2M2} utilizes multiple 
    computer CPUs (if such exist on the system) in order to speed up 
    the calculations. Computations are then executed as parallel 
    processes. The default is \code{FALSE}.}
    
  \item{save.metadata}{a logical specifying if \pkg{P2C2M2} saves 
    the metadata of the analysis to the output variable. The default 
    is \code{FALSE}.}

  \item{verbose}{a logical specifying if \pkg{P2C2M2} prints status 
    information to the screen. The default is \code{FALSE}.}

  \item{dbg}{a logical specifying if \pkg{P2C2M2} is to be run in a 
    debug mode. If \code{TRUE}, then only the first 5 percent of 
    input trees are analyzed and information useful for debugging is 
    printed to the screen. Argument \code{dbg = TRUE} must be set in 
    combination with argument \code{verbose = TRUE}. The default is 
    \code{FALSE}. This argument is intended for developers and should 
    not be selected by regular users.}
}

\value{
  The results of a \pkg{P2C2M2} run comprise test statistics, measures 
  of data dispersion and deviations marked at several quantile levels 
  (analogous to P-values under different alpha-levels in a parametric 
  simulation) for each gene under study and of the sum of all genes. If
  error rates are estimated (using the option error.rate = TRUE), the
  proportion of pods that significantly reject the null, multi-species
  coalescent model are calculated for each gene and across genes.
}

\author{
  Arley Camargo, Joaquin Villamil

  Maintainer: Arley Camargo \email{arley.camargo@gmail.com}
}

\references{
  Cummings, M.P., Neel, M.C. and Shaw, K.L. (2008) A genealogical approach to quantifying lineage divergence. \emph{Evolution}, \bold{62}, 2411--2422.
  
  Gruenstaeudl, M., Reid, N.M., Wheeler, G.R. and Carstens, B.C., 2016. Posterior Predictive Checks of Coalescent Models: P2C2M, an R package. \emph{Molecular Ecology Resources}, \bold{6},193--205. 
doi: 10.1111/1755-0998.12435
  
  Heled, J. and Drummond, A.J. (2010) Bayesian inference of species trees from multilocus data. \emph{Molecular Biology And Evolution}, \bold{27}, 570--580.
  
  Maddison, W.P. (1997) Gene trees in species trees. \emph{Systematic Biology}, \bold{46}, 523--536.
  
  Rannala, B. and Yang, Z. (2003) Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences From Multiple Loci. \emph{Genetics}, \bold{164}, 1645--1656.
}

\examples{
## Example of the minimal data requirements to run P2C2M2

# The absolute path to the input directory is set
inPath <- system.file("extdata", "thomomys/", package="P2C2M2")

# The name of the xml-file generated by BEAUTi and located in 
# "inPath" is set
inFile <- "thomomys.xml"

# Posterior predictive simulations with a setting of 2 simulation 
# replicates are preformed
thomomys <- p2c2m2.complete(inPath, inFile, num.reps=2, 
  error.rate = TRUE, save.metadata=TRUE)
}

