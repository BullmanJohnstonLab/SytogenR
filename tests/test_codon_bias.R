source(file.path("R", "codon_bias.R"))

result <- calculate_codon_bias(c("GCT", "GCC", "GCA", "GCT"))

stopifnot(all(c("codon_counts", "amino_acid_counts", "codon_bias", "preferred_codons") %in% names(result)))
stopifnot(result$codon_counts["GCT"] == 2)
stopifnot(result$codon_counts["GCC"] == 1)
stopifnot(result$codon_counts["GCA"] == 1)
stopifnot(result$preferred_codons[["A"]] == "GCT")

cat("All codon bias tests passed.\n")
