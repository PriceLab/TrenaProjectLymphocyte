library(trenaSGM)
library(TrenaProjectLymphocyte)
library(org.Hs.eg.db)

genome <- "hg38"
targetGene <- "LAG3"

tp <- TrenaProjectLymphocyte();

setTargetGene(tp, targetGene)
tbl.info <- getTranscriptsTable(tp)

chromosome <- tbl.info$chrom
tss <- tbl.info$tss
   # strand-aware start and end: atf1 is on the + strand
start <- tss - 5000
end   <- tss + 5000

tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
mtx <- getExpressionMatrix(tp, "GTEx.lymphocyte.geneSymbols.matrix.asinh")

fpdbs <- c("lymphoblast_hint_16", "lymphoblast_hint_20", "lymphoblast_wellington_16", "lymphoblast_wellington_20")[1:2]

build.spec <- list(title="unit test on LAG3",
                   type="footprint.database",
                   regions=tbl.regions,
                   geneSymbol=targetGene,
                   tss=tss,
                   matrix=mtx,
                   db.host="khaleesi.systemsbiology.net",
                   db.port=5432,
                   databases=fpdbs,
                   annotationDbFile=dbfile(org.Hs.eg.db),
                   motifDiscovery="builtinFimo",
                   tfPool=allKnownTFs(),
                   tfMapping="MotifDB",
                   tfPrefilterCorrelation=0.1,
                   orderModelByColumn="rfScore",
                   solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
suppressWarnings(x <- build(fpBuilder))
lapply(x, dim)
tbls.trimmed <- trimModel(x$model, x$regulatoryRegions, 3)
lapply(tbls.trimmed, dim)
tbl.model <- tbls.trimmed$model
tbl.reg   <- tbls.trimmed$regulatoryRegions

