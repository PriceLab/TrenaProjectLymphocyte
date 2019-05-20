max <- 22330
starts <- seq(1, max, 20)
ends <- starts + 19
ends[length(ends)] <- max
bashScript <- file("runByChunks.sh")
lines <- sprintf("Rscript runMany.R %d %d", starts, ends)
writeLines(lines, bashScript)
close(bashScript)
