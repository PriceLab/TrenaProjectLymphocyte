library(RUnit)
library(trenaSGM)
library(TrenaProjectLymphocyte)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "C")
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_OR4G11P()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_OR4G11P <- function()
{
   printf("--- test_OR4611P")

   proj.L <- TrenaProjectLymphocyte()
   targetGene <- "ENSG00000240361"  # "OR4G11P"
   setTargetGene(proj.L, targetGene)
   xt <- getTranscriptsTable(proj.L)

   tbl.regions <- data.frame(chrom=xt$chrom, start=xt$end-1000, end=xt$end+2000, stringsAsFactors=FALSE)
   mtx.L <- getExpressionMatrix(proj.L, "GTEX.wholeBlood.rna-seq-geneSymbols.filtered")

   build.spec <- list(title="fp.2000up.1000down",
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=xt$tss,
                      matrix=mtx.L,
                      db.host="khaleesi.systemsbiology.net",
                      db.port=5432,
                      databases="lymphoblast_hint_20",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(getGenome(proj.L), targetGene,  build.spec, quiet=TRUE)
   x <-build(fpBuilder)
   lapply(x, dim)
   which(x$model$pearsonCoeff < 0)  # just #66
   tfs <- allKnownTFs()
   length(tfs)
   tfs.in.mtx <- intersect(tfs, rownames(mtx.L))
   length(tfs.in.mtx)  # 885
   cors <- lapply(tfs.in.mtx, function(tf) cor(mtx.L[tf,], mtx[targetGene,]))
   names(cors) <- tfs.in.mtx
   fivenum(as.numeric(cors)) # [1] -0.2420221  0.1791781  0.3808180  0.5242496  0.7565293
   length(which(cors < 0))   #37

} # test_OR4611P
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
