source(file.path("R", "Motif_finder.R"))

results <- find_motifs("GATCGATC", c("GATC", "ATCG"))
stopifnot(nrow(results) == 6)
stopifnot(all(c("motif", "start", "end", "strand", "match") %in% names(results)))
stopifnot(any(results$motif == "GATC" & results$strand == "+" & results$start == 1))
stopifnot(any(results$motif == "GATC" & results$strand == "-" & results$start == 1))

cat("All motif finder tests passed.\n")
