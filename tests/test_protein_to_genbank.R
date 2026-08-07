source(file.path("R", "sequence_io.R"))

# ── 1. Unique hit (best policy, default) ────────────────────────────────────
fna_path <- tempfile(fileext = ".fna")
faa_path <- tempfile(fileext = ".faa")
writeLines(c(">genome1", "ATGGCTGCTTAA"), fna_path)
writeLines(c(">geneA product=example protein", "MAA"), faa_path)

record <- protein_fasta_to_genbank(faa_path, fna_path, topology = "linear")
stopifnot(record$sequence == "ATGGCTGCTTAA")
stopifnot(nrow(record$features) == 1)
stopifnot(record$features$type[1] == "CDS")
stopifnot(record$features$start[1] == 1)
stopifnot(record$features$end[1] == 9)
stopifnot(record$features$strand[1] == "+")
stopifnot(record$features$gene[1] == "geneA")
stopifnot(!record$features$partial[1])
stopifnot(!record$features$multi_hit[1])

# ── 2. Multi-hit: policy="all" returns every occurrence ──────────────────────
# MAAM occurs twice in this genome
fna_multi <- tempfile(fileext = ".fna")
faa_multi <- tempfile(fileext = ".faa")
# ATGGCTGCTATGGCTGCT = MAAMAAA (frame 0) – two non-overlapping MAA matches
writeLines(c(">genome_multi", "ATGGCTGCTATGGCTGCT"), fna_multi)
writeLines(c(">geneB", "MAA"), faa_multi)

record_all <- protein_fasta_to_genbank(faa_multi, fna_multi, multi_hit_policy = "all")
stopifnot(nrow(record_all$features) >= 2)

# ── 3. Multi-hit: policy="best" picks earliest forward hit ───────────────────
record_best <- protein_fasta_to_genbank(faa_multi, fna_multi, multi_hit_policy = "best")
stopifnot(nrow(record_best$features) == 1)
stopifnot(record_best$features$start[1] == 1)
stopifnot(!record_best$features$multi_hit[1])

# ── 4. Multi-hit: policy="strict" flags all hits ─────────────────────────────
record_strict <- suppressWarnings(
  protein_fasta_to_genbank(faa_multi, fna_multi, multi_hit_policy = "strict")
)
stopifnot(nrow(record_strict$features) >= 2)
stopifnot(all(record_strict$features$multi_hit))

# ── 4b. Overlapping multi-hit behavior is preserved ─────────────────────────
fna_overlap <- tempfile(fileext = ".fna")
faa_overlap <- tempfile(fileext = ".faa")
# GCT x 5 encodes AAAAA; query AAAA should match at aa positions 1 and 2.
writeLines(c(">genome_overlap", "GCTGCTGCTGCTGCT"), fna_overlap)
writeLines(c(">geneOverlap", "AAAA"), faa_overlap)

record_overlap_all <- protein_fasta_to_genbank(faa_overlap, fna_overlap, multi_hit_policy = "all")
stopifnot(nrow(record_overlap_all$features) >= 2)
forward_overlap <- record_overlap_all$features[record_overlap_all$features$strand == "+", , drop = FALSE]
stopifnot(nrow(forward_overlap) >= 2)
starts_forward <- sort(as.integer(forward_overlap$start))
ends_forward <- sort(as.integer(forward_overlap$end))
stopifnot(all(c(1L, 4L) %in% starts_forward))
stopifnot(all(c(12L, 15L) %in% ends_forward))

record_overlap_best <- protein_fasta_to_genbank(faa_overlap, fna_overlap, multi_hit_policy = "best")
stopifnot(nrow(record_overlap_best$features) == 1)
stopifnot(as.integer(record_overlap_best$features$start[1]) == 1L)

record_overlap_strict <- suppressWarnings(
  protein_fasta_to_genbank(faa_overlap, fna_overlap, multi_hit_policy = "strict")
)
stopifnot(nrow(record_overlap_strict$features) >= 2)
stopifnot(all(record_overlap_strict$features$multi_hit))

# ── 5. Partial match: allow_partial=TRUE recovers a fragment ─────────────────
fna_partial <- tempfile(fileext = ".fna")
faa_partial <- tempfile(fileext = ".faa")
# genome encodes MAA; protein is MAAK (K not in genome) → partial match = MAA
writeLines(c(">genome_partial", "ATGGCTGCTTAA"), fna_partial)
writeLines(c(">geneC", "MAAK"), faa_partial)

record_partial <- protein_fasta_to_genbank(
  faa_partial, fna_partial, allow_partial = TRUE, min_partial_fraction = 0.5
)
stopifnot(nrow(record_partial$features) == 1)
stopifnot(record_partial$features$type[1] == "CDS_partial")
stopifnot(record_partial$features$partial[1])

# ── 6. No match and no partial fallback returns empty features ───────────────
faa_nomatch <- tempfile(fileext = ".faa")
writeLines(c(">geneD", "CCCCCC"), faa_nomatch)
record_nomatch <- protein_fasta_to_genbank(faa_nomatch, fna_path, allow_partial = FALSE)
stopifnot(nrow(record_nomatch$features) == 0)

# ── 7. convert_sequence_file + read_faa_fna_to_genbank convenience wrappers ──
text_output <- convert_sequence_file(faa_path, output_format = "genbank", gff3_path = fna_path)
stopifnot(is.list(text_output))
stopifnot(text_output$sequence == "ATGGCTGCTTAA")
stopifnot(nrow(text_output$features) == 1)

output_path <- tempfile(fileext = ".gbk")
converted <- read_faa_fna_to_genbank(faa_path, fna_path, output_path = output_path)
stopifnot(file.exists(output_path))
stopifnot(converted$sequence == "ATGGCTGCTTAA")
stopifnot(nrow(converted$features) == 1)

cat("All protein-to-GenBank tests passed.\n")
