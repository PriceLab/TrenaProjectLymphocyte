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
new.rownames <- sub("\\.[0-9]", "", rownames(mtx))
rownames(mtx) <- new.rownames

mtx[1:10, 1:10]

save(mtx, file="../../../inst/extdata/expression/GTEX.wholeBlood.rna-seq.RData")

