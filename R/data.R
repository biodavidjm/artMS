#' Random data set
#'
#' Dataset randomly generated for testing purposes
#'
#' @format A data frame with 100 rows and 10 variables:
#' \describe{
#' Dataset generated using this code
#' 
#' `data.frame(replicate(10,sample(0:1,100,rep=TRUE)))`
#' }
"randomDF"

#' MaxQuant Phosphorylation evidence file
#' 
#' Phosphorylation data set generated on Human cell lines
#' 
#' Two different cell lines analyzed:  `Cal33` and `HSC6`
#' 
#' 4 replicates per cell line
#' 
#' The information about the experimental design of the evidence file can 
#' be found in the data.frame also available in this package named `ph_keys`
#' 
#' @format tab delimited file
#' \describe{
#' \item{ph_evidence}{The MaxQuant evidence file. Check this link to find out 
#' more about it:
#' \url{http://www.coxdocs.org/doku.php?id=maxquant:start}}
#' }
"ph_evidence"

#' Keys file: Experimental design describing the `ph_evidence` file
#' 
#' Data frame with the information about the experimental design.
#' It provides the metadata of the `ph_evidence` data.frame also available in 
#' this package (see `vignettes` to find out more)
#' 
#' Phosphorylation data set generated on Human cell lines
#' 
#' Two different cell lines analyzed:  Cal33 and HSC6
#' 
#' 4 replicates per cell line
#' 
#' @format Tab delimited file with the following columns:
#' \describe{
#' \item{Raw.file}{Raw file processed. Each one should be either a
#' biological replicate or a technical replicate}
#' 
#' \item{IsotopeLabelType}{Type of labeling. `L` is used for label free 
#' experiments}
#' 
#' \item{Condition}{Label for conditions. Only alpha-numeric characters and 
#' `underscore (_)` are allowed}
#' 
#' \item{BioReplicate}{Label for the Biological replicates. The consensus is
#' to use the same label as the Condition, but adding a `dash (-)` corresponding
#' to the number of biological replicate. For example, `Cal-1`, `Cal-2`, 
#' `Cal-3`}
#' }
"ph_keys"


