
biocGet <- function (pkgs){
   library(BiocManager)
   BiocManager::install(pkgs)
   }

pkgs.needed <- c("devtools",
                 "motifmatchr",
                 "TFBSTools",
                 "universalmotif",
                 "ggseqlogo",
                 "yaml",
                 "later",
                 "org.Hs.eg.db",
                 "AnnotationDbi",
                 "Biostrings",
                 "XVector",
                 "glmnet",
                 "foreach",
                 "Matrix",
                 "GenomicRanges",
                 "GenomeInfoDb",
                 "IRanges",
                 "S4Vectors",
                 "BiocGenerics",
                 "DT",
                 "V8",
                 "shiny",
                 "shinydashboard",
                 "shinyjs",
                 "colourpicker",
                 "RColorBrewer")

biocGet(pkgs.needed)

library(devtools)
install_github("paul-shannon/igvShiny")
install_github("paul-shannon/cyjShiny")
