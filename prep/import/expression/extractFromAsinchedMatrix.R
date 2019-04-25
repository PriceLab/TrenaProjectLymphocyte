library(TrenaProject)
print(load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData")))
library(trenaSGM)  # for allKnownTFs
print(load("mtx.gtex.asinh.RData"))
mtx <- mtx.gtex.a
rownames(mtx) <- sub("\\.[0-9]", "", rownames(mtx))
dim(mtx)

tbl.md <- read.table("GTEx_v7_Annotations_SampleAttributesDS.txt", sep = '\t', header=TRUE, quote="", as.is=TRUE, row.names=1)
dim(tbl.md) # 15598 63
skin.samples <- rownames(tbl.md[grep("skin", tbl.md$SMTS, ignore.case=TRUE), ])
length(skin.samples)  # 1362
skin.samples <- gsub("-", ".", skin.samples, fixed=TRUE)
keepers <- intersect(colnames(mtx), skin.samples)
length(keepers) # 1203
mtx.skin <- mtx[, keepers]
dim(mtx.skin)

blood.samples <- rownames(tbl.md[grep("^Blood$", tbl.md$SMTS, ignore.case=TRUE), ])
length(blood.samples)
blood.samples <- gsub("-", ".", blood.samples, fixed=TRUE)
keepers <- intersect(colnames(mtx), blood.samples)
length(keepers) # 537
mtx.blood <- mtx[, keepers]
dim(mtx.blood)   #56202 537

wholeBlood.samples <- rownames(tbl.md[grep("^Whole_Blood", tbl.md$SMTSD, ignore.case=TRUE), ])

length(wholeBlood.samples)
wholeBlood.samples <- gsub("-", ".", wholeBlood.samples, fixed=TRUE)
keepers <- intersect(colnames(mtx), wholeBlood.samples)
length(keepers) # 407
mtx.wholeBlood <- mtx[, keepers]
dim(mtx.wholeBlood)   #56202 537


#target.gene <- "COL1A1"
target.gene <- "BEND3"
target.ensg <- subset(tbl.geneInfo, geneSymbol==target.gene)$ensg
tfs.ensg <- allKnownTFs(identifierType="ensemblGeneID")
tfs.ensg.in.mtx <- intersect(tfs.ensg, rownames(mtx.skin))
length(tfs.ensg.in.mtx)  # 1024

mtx.target <- mtx.skin
mtx.target <- mtx.blood
mtx.target <- mtx.wholeBlood

cors <- lapply(tfs.ensg.in.mtx, function(tf) cor(mtx.target[tf,], mtx.target[target.ensg,]))
hist(as.numeric(cors), main=sprintf("all blood: %s", target.gene))
hist(mtx, main="mtx asinh")



matrix.colnames <- colnames(mtx)
new.colnames <- gsub(".", "-", matrix.colnames, fixed=TRUE)
length(intersect(new.colnames, rownames(tbl.md)))

tbl.tissues <- as.data.frame(table(tbl.md$SMTSD), stringsAsFactors=FALSE)
colnames(tbl.tissues) <- c("tissue", "sampleCount")
descending.order <- order(tbl.tissues$sampleCount, decreasing=TRUE)
tbl.tissues <- tbl.tissues[descending.order,]

head(tbl.tissues)   # see that Whole_Blood has 2412 (alleged) samples

wholeBloodSamples <- grep("^Whole_Blood", tbl.md$SMTSD, ignore.case=TRUE)
length(wholeBloodSamples)         # 2412
bloodSamples <- grep("^Blood", tbl.md$SMTS, ignore.case=TRUE)
length(bloodSamples)              # 3602
length(intersect(wholeBloodSamples, bloodSamples))  # 2412

id.oi <- rownames(subset(tbl.md, SMTSD=="Whole_Blood"))
length(id.oi)   # 2412
id.oi <- intersect(id.oi, new.colnames)    # but not all are actually in the expression matrix
length(id.oi)   # 407
id.oi <- gsub("-", ".", id.oi, fixed=TRUE) # now they match the spelling (with dots, not dashses) found in mtx
length(intersect(id.oi, colnames(mtx)))

# extract from the big all-tissue matrix
mtx <- mtx[, id.oi]
dim(mtx)   # 56202 407
new.rownames <- sub("\\.[0-9]", "", rownames(mtx))
rownames(mtx) <- new.rownames

mtx[1:10, 1:10]
fivenum(mtx)
mtx <- asinh(mtx)
fivenum(mtx)

means <- apply(mtx, 1, mean)
as.integer(fivenum(means))  #   0  0  0  3 16
length(which(means > 0.5))  # 26k
keepers <- names(which(means > 0.5))
mtx <- mtx.asinh[keepers,]
dim(mtx)   # 26331
save(mtx, file="../../../inst/extdata/expression/GTEX.wholeBlood.rna-seq.filtered.RData")

load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))
dim(tbl.geneInfo)

matches <- match(rownames(mtx), tbl.geneInfo$ensg)
length(matches) # 26331
syms <- tbl.geneInfo$geneSymbol[matches]
rownames(mtx) <- syms
deleters <- which(is.na(rownames(mtx)))
length(deleters)
mtx <- mtx[-deleters,]
dim(mtx)

mtx[1:10, 1:10]
save(mtx, file="../../../inst/extdata/expression/GTEX.wholeBlood.rna-seq-geneSymbols.filtered.RData")






