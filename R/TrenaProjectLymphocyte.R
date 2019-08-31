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
                                  contains="TrenaProjectHG38")

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
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE, comment.char="#")
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   footprintDatabaseNames <- c("lymphoblast_hint_16",  "lymphoblast_hint_20", "lymphoblast_wellington_16", "lymphoblast_wellington_20")
   dataDirectory <- system.file(package="TrenaProjectLymphocyte", "extdata")
   footprintDatabaseHost <- "khaleesi.systemsbiology.net"
   footprintDatabasePort <- 5432

   covariatesFile <- NA_character_;

   stopifnot(file.exists(dataDirectory))

     # TODO: this should be hidden by, provided by TrenaProjectHG38:
   geneInfoTable.path <- system.file(package="TrenaProjectHG38", "extdata", "geneInfoTable.RData")

   .TrenaProjectLymphocyte(TrenaProjectHG38(projectName="TrenaProject Lymphocyte",
                                            supportedGenes=geneSets[[2]],
                                            footprintDatabaseHost=footprintDatabaseHost,
                                            footprintDatabasePort=footprintDatabasePort,
                                            footprintDatabaseNames=footprintDatabaseNames,
                                            packageDataDirectory=dataDirectory,
                                            quiet=quiet
                                            ))

} # TrenaProjectLymphocyte, the constructor
#----------------------------------------------------------------------------------------------------
