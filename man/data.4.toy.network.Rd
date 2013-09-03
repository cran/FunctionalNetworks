\name{data.4.toy.network}
\alias{data.4.toy.network}
\docType{data}
\title{List containing the results from the application of the function \code{data.generation} for a toy dataset based on the \code{breastCancerVDX} dataset}
\description{
This dataset contains a list with the results from the application of the function \code{data.generation} for a toy dataset based on the \code{breastCancerVDX} dataset. It is based on the Molecular Functions ontology from GO by the MSigBD. Only MF with at least 5 genes were considered. See the package Vignette for a detailed usage.}
\usage{data(data.4.toy.network)
}
\format{
  A list containing the information in the following way:
  \itemize{
   \item \code{gene.data} A matrix with gene expression data. Rows are genes (HUGO ids) and columns are samples.
   \item \code{set.data} A matrix with expression data for the gene sets. Rows are gene sets and columns are samples.
   \item \code{affy.loc} A vector that contains the row numbers corresponding to the genes that are reported in the annotation system chosen
   \item \code{gene2set.mat} A matrix that where each row represents a gene and each column contains the gene set to which this gene is reported based on the annotation system chosen.
   \item \code{set2gene.mat} A matrix that where each row represents a gene set and each column contains the genes reported to be in the set based on the annotation system chosen.
   \item \code{Set.obj} A vector with the gene sets set as objective nodes in the network. 
   \item \code{Set.src} A vector with the gene sets set as source nodes in the network. 
   \item \code{G.obj} A vector with the genes set as objective nodes in the network. 
   \item \code{G.src} A vector with the genes set as source nodes in the network.
  }
}
\details{
This data set is calculated for the particular case were it is of interest to obtain the Functional Network for a toy dataset based on the \code{breastCancerVDX} dataset and the annotation system chosen is the MF from GO. The phenotype chose are the samples with Estrogen receptor negative (ER-).
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
\keyword{gmt file}
\keyword{Network}
\keyword{Gene expression data}