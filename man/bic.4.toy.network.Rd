\name{bic.4.toy.network}
\alias{bic.4.toy.network}
\docType{data}
\title{List containing the results from the application of the function \code{bic.generation} for a toy dataset based on the \code{breastCancerVDX} dataset}
\description{
This dataset contains a list with the results from the application of the function \code{bic.generation}  for a toy dataset based on the \code{breastCancerVDX} dataset. It is based on the results from the application of the function \code{data.generation}. See the package Vignette for a detailed usage.}
\usage{data(bic.4.toy.network)
}
\format{
  A list containing the information in the following way:
  \itemize{
   \item \code{BIC.Gene.0.Pred} A vector with the BIC calculations were genes are set as objective nodes in the network with no source gene.
   \item \code{BIC.Gene.1.Pred} A matrix with the BIC calculations were genes are set as objective nodes in the network with one source gene.
   \item \code{BIC.Set.0.Pred} A vector with the BIC calculations were gene sets are set as objective nodes in the network with no source gene sets.
   \item \code{BIC.Set.1.Pred} A matrix with the BIC calculations were gene sets are set as objective nodes in the network with one source gene set.
  }
}
\details{
This data set is calculated for the particular case were it is of interest to obtain the Functional Network  for a toy dataset based on the \code{breastCancerVDX} dataset and the annotation system chosen is the MF from GO. The phenotype chose are the samples with Estrogen receptor negative (ER-) and the minimum number of genes on each MF is at least 5.
} 
%\itemize{
%}
%\source{
%url{http://www.broadinstitute.org/gsea/msigdb/collections.jsp#C5}
%}
\references{
Quiroz-Zarate A, Haibe-Kains, B and   Quackenbush J (2013). "Manuscript in preparation".
}
\author{
A. Quiroz-Zarate.
\email{aquiroz@jimmy.harvard.edu}
}
\examples{
# See Package Vignette for a detailed example of the usage of this dataset
# vigenette("FunctionalNetworks")
}
\keyword{Bayesian Information Criteria}
\keyword{Gene sets}