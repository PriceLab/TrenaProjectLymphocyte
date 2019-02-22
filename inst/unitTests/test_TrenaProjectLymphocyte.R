library(TrenaProjectLymphocyte)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tProj")) {
   message(sprintf("--- creating instance of TrenaProjectLymphocyte"))
   tProj <- TrenaProjectLymphocyte();
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

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectLymphocyte", "TrenaProject") %in% is(tProj)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("Abca1")
   checkTrue(all(subset.expected %in% getSupportedGenes(tProj)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(getVariantDatasetNames(tProj), character(0))

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c("extraembryonic_structure_wellington_16", "extraembryonic_structure_wellington_20",
                 "extraembryonic_structure_hint_16", "extraembryonic_structure_hint_20")
   checkTrue(is.na(getFootprintDatabaseNames(tProj)))
   checkTrue(is.na(getFootprintDatabaseHost(tProj)))

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("thioglycollate-elicited-peritoneal-macrophages")
   checkTrue(all(expected %in% getExpressionMatrixNames(tProj)))

   mtx <- getExpressionMatrix(tProj, expected[1])
   checkEquals(dim(mtx), c(6890, 14))

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

   setTargetGene(tProj, "Abca1")
   checkEquals(getTargetGene(tProj), "Abca1")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tProj)
   checkTrue(nrow(tbl.transcripts) >= 1)
   checkEquals(tbl.transcripts$chr, "chr4")
   checkEquals(tbl.transcripts$start, 53030787)
   checkEquals(tbl.transcripts$end , 53159895)
   checkEquals(tbl.transcripts$tss, 53159895)
   checkEquals(tbl.transcripts$strand, -1)

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tProj, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr4:53030787-53159895")

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tProj)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) >= 0)

   message(sprintf("    geneGeneEnhancersRegion"))
   region <- getGeneEnhancersRegion(tProj, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr4:52901678-53289004")

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tProj)
   checkEquals(nrow(tbl.dhs), 0)

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- getChipSeq(tProj, chrom=chromosome, start=start, end=end, tfs=NA)
   checkEquals(nrow(tbl.chipSeq), 0)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
