library(TrenaProjectLymphocyte)
library(RUnit)
library(org.Hs.eg.db)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")) {
   message(sprintf("--- creating instance of TrenaProjectLymphocyte"))
   tp <- TrenaProjectLymphocyte();
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_supportedGenes()
   test_variants()
   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()
   test_buildSingleGeneModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectLymphocyte", "TrenaProjectHG38", "TrenaProject") %in% is(tp)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("EOMES", "IL6")
   checkTrue(all(subset.expected %in% getSupportedGenes(tp)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

     # this variant file mentioned here is very large, not kept in the github repo
     # and thus not in any but the "home" machine for this class - currently pshannon's
     # riptide

   if(file.exists(system.file(package="TrenaProjectLymphoctye", "extdata", "variants",
                              "sebastiani.2017.gwas.RData")))
       checkTrue("sebastiani.2017.gwas" %in% getVariantDatasetNames(tp))

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c("lymphoblast_hint_16", "lymphoblast_hint_20", "lymphoblast_wellington_16", "lymphoblast_wellington_20")

   checkTrue(all(expected %in% getFootprintDatabaseNames(tp)))
   checkEquals(getFootprintDatabaseHost(tp), "khaleesi.systemsbiology.net")

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("GTEX.wholeBlood.rna-seq-geneSymbols.22330x407",
                 "GTEX.lymphocyte.rna-seq-geneSymbols.21415x130")

   checkTrue(all(expected %in% getExpressionMatrixNames(tp)))

   mtx <- getExpressionMatrix(tp, expected[1])
   checkEquals(dim(mtx), c(22330, 407))
   checkEquals(head(sort(rownames(mtx)), n=3), c("A1BG", "A1BG-AS1", "A1CF"))
   checkTrue(max(mtx) < 100)

   mtx <- getExpressionMatrix(tp, expected[2])
   checkEquals(dim(mtx), c(21415, 130))
   checkEquals(head(sort(rownames(mtx)), n=3), c("A1BG", "A1BG-AS1", "A1CF"))
   checkTrue(max(mtx) < 100)

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

   setTargetGene(tp, "EOMES")
   checkEquals(getTargetGene(tp), "EOMES")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tp)
   checkTrue(nrow(tbl.transcripts) == 1)
   checkEquals(tbl.transcripts$chr, "chr3")

   checkEquals(tbl.transcripts$start, 27715949)
   checkEquals(tbl.transcripts$end,  27722711)
   checkEquals(tbl.transcripts$tss, 27722498)
   checkEquals(tbl.transcripts$strand, -1)

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tp, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr3:27715949-27722711")

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tp)
   checkEquals(head(colnames(tbl.enhancers)), c("chrom", "start", "end", "gene", "eqtl", "hic"))

   checkTrue(nrow(tbl.enhancers) >= 0)

   message(sprintf("    geneGeneEnhancersRegion"))
   region <- getGeneEnhancersRegion(tp, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chrom, "chr3")  # start and end likely to change over time with genehancer updates

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tp)
   checkTrue(nrow(tbl.dhs) > 200)

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- with(tbl.transcripts, getChipSeq(tp, chrom=chrom, start=start, end=end, tfs="BCLAF1"))
   checkEquals(nrow(tbl.chipSeq), 1)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel <- function()
{
   printf("--- test_buildSingleGeneModel")

   genome <- "hg38"
   targetGene <- "LAG3"
   setTargetGene(tp, targetGene)
   tbl.info <- getTranscriptsTable(tp)

   chromosome <- tbl.info$chrom
   tss <- tbl.info$tss
      # strand-aware start and end: atf1 is on the + strand
   start <- tss - 5000
   end   <- tss + 5000

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
   matrix.name <-"GTEX.wholeBlood.rna-seq-geneSymbols.22330x407"
   checkTrue(matrix.name %in% getExpressionMatrixNames(tp))

   mtx <- getExpressionMatrix(tp, matrix.name)

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
   checkEquals(sort(names(x)), c("model", "regulatoryRegions"))
   checkTrue(nrow(x$model) > 5)
   some.expected.tfs <- c("RUNX3", "EOMES", "TBX21")
   checkTrue(all(some.expected.tfs %in% head(x$model$gene, n=10)))

} # test_buildSingleGeneModel
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel_IRF4 <- function()
{
   printf("--- test_buildSingleGeneModel_IRF4")

   genome <- "hg38"
   targetGene <- "IRF4"
   setTargetGene(tp, targetGene)
   tbl.info <- getTranscriptsTable(tp)

   chromosome <- tbl.info$chrom
   tss <- tbl.info$tss
      # strand-aware start and end: atf1 is on the + strand
   start <- tss - 5000
   end   <- tss + 5000

   mtx <- getExpressionMatrix(tp, "GTEX.lymphocyte.rna-seq-geneSymbols.21415x130")

   # for a very small model + regions
   # mtx <- getExpressionMatrix(tp, "GTEX.wholeBlood.rna-seq-geneSymbols.22330x407")
   # start <- tss - 2000
   # end <- tss + 200
   # then, with results in hand:
   #   tbl.model <- head(x$model, n=5)
   #   tbl.regRaw <- subset(x$regulatoryRegions, geneSymbol %in% tbl.model$gene)
   #   dups <- which(duplicated(tbl.regRaw$loc))
   #   tbl.reg  <- tbl.regRaw[-dups,]
   #   save(tbl.model, tbl.reg, file="../extdata/model.and.regRegions.irf4.top5.RData")

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   fpdbs <- c("lymphoblast_hint_20")

   build.spec <- list(title="unit test on IRF4",
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
                      motifSpeciesRestriction="hsapiens",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.1,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   suppressWarnings(x <- build(fpBuilder))
   checkEquals(sort(names(x)), c("model", "regulatoryRegions"))
   checkTrue(nrow(x$model) > 5)

} # test_buildSingleGeneModel
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
