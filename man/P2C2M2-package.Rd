\name{P2C2M2-package}
\alias{P2C2M2-package}
\alias{P2C2M2}
\docType{package}

\title{
  Posterior Predictive Checks of Coalescent Models in BEAST2
}

\description{
  \pkg{P2C2M2} provides functions to read default output from BEAST2 
  (Bouckaert et al. 2019) and starBEAST2 (Ogilvie et al. 2017)
  and conduct posterior predictive checks of coalescent models (Reid
  et al. 2014) based on the functions originally provided in the P2C2M
  package (Gruenstaeudl et al. 2016) to parse tree and xml files from
  BEAST and *BEAST. It also implements estimation of type-I (false
  -positive) error rates using pseudo-observed datasets ("pods")
  sampled from the posterior predictive distribution.
}

\note{

  \bold{Installation Instructions}

  To use \pkg{P2C2M2}, the default version of Python must be set to 
  Python 2.7. Users of unix-like operating systems can insure that 
  this requirement is fulfilled by setting the following alias:

  \code{echo 'alias python=python2.7' >> ~/.bashrc}

  Mandatory and optional dependencies of \pkg{P2C2M2} can be installed 
  automatically via two installation scripts that are co-supplied with 
  the package. These scripts were designed for unix-like operating 
  systems and are located in folder /exec. To use these installation 
  scripts, a correct configuration of python2-setuptools is required. 
  Users of unix-like operating systems can insure a correct 
  configuration by setting the following alias:

  \code{echo 'alias python-config=python2-config' >> ~/.bashrc}

  To execute the R installer, please run the following commands in R:

  \code{source('/path_to_P2C2M2/exec/P2C2M.installRlibs.R')}; 
  \code{p2c2m.install()}
  
  To execute the Python installer, please run the following command in 
  a terminal:

  \code{python /path_to_P2C2M2/exec/P2C2M.installPylibs.py}


  \emph{Special Note for MacOS}

  Users of the MacOS operating system need to install the dependencies 
  manually. Prior to their installation, please confirm that file 
  '/usr/bin/python2-config' exists in your file system and that it 
  points to the Python 2.7 executable. Please refer to \url{http://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html} 
  on how to install R packages manually. For the manual installation of 
  Python libraries, please refer to 
  \url{http://docs.python.org/2/using/mac.html}


  \bold{Study Design Requirements}

  In the user-supplied data set, every species should be represented by 
  at least two alleles. Species that are represented by only a single 
  allele, by contrast, must be specified via option 
  \code{"single.allele"} and thereby are not included in the 
  calculation of the summary statistic 'GSI'; misspecifications 
  causes \pkg{P2C2M2} to print the error message \emph{'Error: given 
  group represents one or fewer taxa. Cannot compute index.'}).


  \bold{Input File Requirements}

  In order to execute \pkg{P2C2M2}, a user must provide a directory with 
  three different types of input files: (a) a file that contains 
  species trees, (b) a file that contains gene trees for each gene 
  under study, and (c) an XML-formatted file generated by BEAUTi2, the 
  input generator of BEAST2 (Bouckaert et al. 2019). A species tree 
  file contains a draw of s generations from the posterior distribution 
  of species trees. Each gene tree file contains an equally large draw 
  from the respective posterior distribution of ultrametric 
  genealogies. Please note that the generations recorded in the species 
  tree file must match those in the gene tree files exactly. The input 
  file generated by BEAUTi2 is formatted in XML markup language and 
  represents the starting point for a species tree inference in
  starBEAST2. Here, it provides information on allele and species names,
  the association between alleles and species, and ploidy levels to 
  \pkg{P2C2M2}.


  \bold{File Name Requirements}

  The following requirements for input file names are in place: The 
  species tree file must be named 'species.trees'. Each gene tree file 
  must be named 'g.trees', where the letter g is substituted with the 
  actual name of the gene. The name of the xml-formatted input file is 
  not constrained and at the discretion of the user. Please be aware 
  that \pkg{P2C2M2} uses the name of the xml-formatted input file name 
  to label all subsequent output of the package.
}

\author{
  Arley Camargo, Joaquin Villamil

  Maintainer: Arley Camargo \email{arley.camargo@gmail.com}
}

\references{
  
  Bouckaert, R., Vaughan, T.G., Barido-Sottani, J., Duchêne, S., Fourment, M., Gavryushkina, A., Heled, J., Jones, G., Kühnert, D., De Maio, N., Matschiner, M., Mendes, F.K., Müller, N.F., Ogilvie, H.A., du Plessis, L., Popinga, A., Rambaut, A., Rasmussen, D., Siveroni, I., Suchard, M.A., Wu, C.H., Xie, D., Zhang, C., Stadler, T., Drummond, A.J. (2019) BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. \emph{PLoS Computational Biology}, \bold{15}, e1006650.
  
  Drummond, A.J. and Rambaut, A. (2007) BEAST: Bayesian evolutionary analysis by sampling trees. \emph{BMC Evolutionary Biology}, \bold{7}, 214.

  Gruenstaeudl, M., Reid, N.M., Wheeler, G.R. and Carstens, B.C. (2016) Posterior predictive checks of coalescent nodels: P2C2M, an R package. \emph{Molecular Ecology Resources}, \bold{6},193--205.

  Heled, J. and Drummond, A.J. (2010) Bayesian inference of species trees from multilocus data. \emph{Molecular Biology and Evolution}, \bold{27}, 570--580.
  
  Ogilvie, H.A., Bouckaert, R.R., Drummond, A.J. (2017) StarBEAST2 brings faster species tree inference and accurate estimates of substitution rates. \emph{Molecular Biology and Evolution}, \bold{34}, 2101-–2114.

  Reid, N.M., Brown, J.M., Satler, J.D., Pelletier, T.A., McVay, J.D., Hird, S.M. and Carstens, B.C. (2014) Poor fit to the multi-species coalescent model is widely detectable in empirical data. \emph{Systematic Biology}, \bold{63}, 322--333.
}

