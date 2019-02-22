#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectLymphocyte-class
#'
#' @name TrenaProjectLymphocyte-class
#' @rdname TrenaProjectLymphocyte-class
#' @aliases TrenaProjectLymphocyte
#' @exportClass TrenaProjectLymphocyte
#'

.TrenaProjectLymphocyte <- setClass("TrenaProjectLymphocyte",
                                  contains="TrenaProject")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectLymphocyte
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectLymphocyte-class
#'
#' @export
#'
#' @return An object of the TrenaProjectLymphocyte class
#'

TrenaProjectLymphocyte <- function(quiet=TRUE)

{
   genomeName <- "hg38"

   directory <- system.file(package="TrenaProjectLymphocyte", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()
   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE)
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   footprintDatabaseNames <- NA_character_;
   expressionDirectory <- system.file(package="TrenaProjectLymphocyte", "extdata", "expression")
   variantsDirectory <- system.file(package="TrenaProjectLymphocyte", "extdata", "variants")
   footprintDatabaseHost <- NA_character_;

   covariatesFile <- NA_character_;

   stopifnot(file.exists(expressionDirectory))

   .TrenaProjectLymphocyte(TrenaProject(supportedGenes=geneSets[[1]],
                                        genomeName=genomeName,
                                        footprintDatabaseHost=footprintDatabaseHost,
                                        footprintDatabaseNames=footprintDatabaseNames,
                                        expressionDirectory=expressionDirectory,
                                        variantsDirectory=variantsDirectory,
                                        covariatesFile=covariatesFile,
                                        quiet=quiet
                                        ))

} # TrenaProjectLymphocyte, the constructor
#----------------------------------------------------------------------------------------------------
