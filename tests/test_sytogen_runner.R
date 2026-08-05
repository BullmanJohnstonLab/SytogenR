source(file.path("R", "parse_motif.R"))
source(file.path("R", "codon_bias.R"))
source(file.path("R", "Motif_finder.R"))
source(file.path("R", "SyToGen.R"))

motif_text <- "motif\nGATC\n"
input_sequence <- "ATGGATCCGCTGATC"
cds_ranges <- data.frame(start = 1L, end = 15L, stringsAsFactors = FALSE)

result <- run_sytogen(
  sequence = input_sequence,
  motif_text = motif_text,
  cds_ranges = cds_ranges,
  include_reverse = TRUE
)

stopifnot(all(c("sequence", "motifs", "motif_hits", "codon_bias", "decision_matrix", "mutated_sequence") %in% names(result)))
stopifnot(nrow(result$motif_hits) >= 1)
stopifnot(nrow(result$decision_matrix) >= 1)
stopifnot(result$mutated_sequence != result$sequence)
stopifnot(!identical(result$mutated_sequence, ""))

runner_result <- sytogen_runner(
  sequence = input_sequence,
  motifs = c("GATC"),
  cds_ranges = cds_ranges,
  include_reverse = TRUE
)

stopifnot(nrow(runner_result$decision_matrix) >= 1)
stopifnot(runner_result$mutated_sequence != runner_result$sequence)

cat("All SyToGen runner tests passed.\n")
