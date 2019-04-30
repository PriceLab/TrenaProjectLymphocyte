biocGet <- function (pkgs){
   library(BiocManager)
   BiocManager::install(pkgs)
   }

pkgs.needed <- c("motifmatchr",
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
                 "shinydashboard",
                 "shinyjs",
                 "shiny",
                 "colourpicker",
                 "RColorBrewer")

biocGet(pkgs.needed)

