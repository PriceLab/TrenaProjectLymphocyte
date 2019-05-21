system("scp pshannon@whovian:/local/Cory/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt .")
system("scp pshannon@whovian:/local/Cory/gtex/GeneReads.rds .")
mtx.gtex <- readRDS("GeneReads.rds")

tbl.md <- read.table("GTEx_v7_Annotations_SampleAttributesDS.txt", sep = '\t', header=TRUE,
                     quote="", as.is=TRUE, row.names=1)

dim(tbl.md) # 15598 63

matrix.colnames <- colnames(mtx.gtex)
new.colnames <- gsub(".", "-", matrix.colnames, fixed=TRUE)
length(intersect(new.colnames, rownames(tbl.md)))

tbl.tissues <- as.data.frame(table(tbl.md$SMTSD), stringsAsFactors=FALSE)
colnames(tbl.tissues) <- c("tissue", "sampleCount")
descending.order <- order(tbl.tissues$sampleCount, decreasing=TRUE)
tbl.tissues <- tbl.tissues[descending.order,]

head(tbl.tissues)   # see that Whole_Blood has 2412 (alleged) samples

id.oi <- rownames(subset(tbl.md, SMTSD=="Whole_Blood"))
length(id.oi)   # 2412
id.oi <- intersect(id.oi, new.colnames)    # but not all are actually in the expression matrix
length(id.oi)   # 407
id.oi <- gsub("-", ".", id.oi, fixed=TRUE) # now they match the spelling (with dots, not dashses) found in mtx.gtex
length(intersect(id.oi, colnames(mtx.gtex)))

# extract from the big all-tissue matrix
mtx <- mtx.gtex[, id.oi]
dim(mtx)   # 56202 407
new.rownames <- sub("\\.[0-9]*", "", rownames(mtx))
rownames(mtx) <- new.rownames

irf4 <- "ENSG00000137265"
tet2 <- "ENSG00000168769"

mtx[irf4,]
mtx[tet2,]
mtx[1:10, 1:10]
fivenum(mtx)
mtx <- asinh(mtx)
fivenum(mtx)

maxes <- apply(mtx, 1, max)
keepers <- which(maxes > 3)
length(keepers)

irf4 %in% names(keepers)
tet2 %in% names(keepers)
# length(intersect(names(keepers), allKnownTFs()))

mtx <- mtx[keepers,]
dim(mtx)   # 26331
save(mtx, file="../../../inst/extdata/expression/GTEX.wholeBlood.rna-seq.filtered.ensg.RData")

load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))
dim(tbl.geneInfo)

matches <- match(rownames(mtx), tbl.geneInfo$ensg)
length(matches) # 26331
syms <- tbl.geneInfo$geneSymbol[matches]
length(intersect(allKnownTFs(), syms))
"IRF4" %in% syms
"TET2" %in% syms

rownames(mtx) <- syms
length(which(is.na(rownames(mtx))))
deleters <- which(is.na(rownames(mtx)))
length(deleters)
mtx <- mtx[-deleters,]
dim(mtx)

all(c("IRF4", "TET2", "LAG3", "PDCD1", "HAVCR2") %in% rownames(mtx))
ensgs.to.delete <- grep("^ENSG", rownames(mtx))
length(ensgs.to.delete)
mtx <- mtx[-ensgs.to.delete,]
dim(mtx)

setdiff(c("IRF4", "TET2", "LAG3", "PDCD1", "HAVCR2"), rownames(mtx))
mtx[1:10, 1:10]
save(mtx, file="../../../inst/extdata/expression/GTEX.wholeBlood.rna-seq-geneSymbols.22330.RData")






