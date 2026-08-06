source(file.path("R", "parse_motif.R"))
source(file.path("R", "codon_bias.R"))
source(file.path("R", "Motif_finder.R"))
source(file.path("R", "SyToGen.R"))

sequence <- "ATGGATCGCTGATC"
motif_df <- data.frame(
  motif = c("GATC", "ATCG"),
  enz_type = c("4", "2"),
  stringsAsFactors = FALSE
)
cds_df <- data.frame(
  codon = c("GAT", "GAC", "GCT", "GCC"),
  frequency = c(0.1, 0.9, 0.2, 0.8),
  stringsAsFactors = FALSE
)
params <- list(
  topology = "linear",
  preserve_gc = TRUE,
  cds_ranges = data.frame(start = 1L, end = 12L, stringsAsFactors = FALSE),
  protected_ranges = data.frame(start = 2L, end = 2L, stringsAsFactors = FALSE)
)

result <- run_sytogen_pipeline(sequence, codon_df = cds_df, motif_df = motif_df, params = params)

stopifnot(all(c("summary", "motif_summary", "decision_matrix", "mutated_sequence", "new_motifs") %in% names(result)))
stopifnot(result$summary$motifs_input >= 1)
stopifnot(nrow(result$decision_matrix) >= 1)
stopifnot(nrow(result$motif_summary) >= 1)
stopifnot(any(result$decision_matrix$chosen))
stopifnot(length(result$new_motifs) == 0)
if (nrow(result$decision_matrix) > 0 && any(result$decision_matrix$edit_position != "")) {
  stopifnot(!any(as.integer(result$decision_matrix$edit_position[result$decision_matrix$chosen]) == 2L))
}
stopifnot(any(result$decision_matrix$skip_reason == "type_iv"))

record_input <- list(
  sequence = sequence,
  id = "record_1",
  features = data.frame(
    type = c("CDS", "misc_feature"),
    start = c(1L, 2L),
    end = c(12L, 2L),
    strand = c("+", "+"),
    gene = c("bla", ""),
    stringsAsFactors = FALSE
  )
)

record_result <- run_sytogen_pipeline(record_input, codon_df = cds_df, motif_df = motif_df, params = list(topology = "linear", preserve_gc = TRUE))

stopifnot(record_result$summary$sequence_id == "record_1")
stopifnot(nrow(record_result$decision_matrix) >= 1)
stopifnot(length(record_result$motif_hits) >= 1)

cat("All SyToGen pipeline tests passed.\n")
