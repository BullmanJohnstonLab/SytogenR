source(file.path("R", "sequence_io.R"))

record <- sequence_record(
  sequence = "ATGGATCCTTAA",
  id = "plasmid1",
  description = "example plasmid",
  topology = "circular",
  features = data.frame(
    type = c("CDS", "misc_feature"),
    start = c(1L, 5L),
    end = c(12L, 7L),
    strand = c("+", "-"),
    gene = c("bla", "ori"),
    stringsAsFactors = FALSE
  )
)

fasta_text <- write_fasta_text(record)
gff3_text <- write_gff3_text(record)
genbank_text <- write_genbank_text(record)

stopifnot(startsWith(fasta_text, ">plasmid1"))
stopifnot(grepl("##gff-version 3", gff3_text, fixed = TRUE))
stopifnot(grepl("LOCUS", genbank_text, fixed = TRUE))

fasta_path <- tempfile(fileext = ".fasta")
gff3_path <- tempfile(fileext = ".gff3")
genbank_path <- tempfile(fileext = ".gbk")
write_fasta_file(record, fasta_path)
write_gff3_file(record, gff3_path)
write_genbank_file(record, genbank_path)

parsed_fasta <- read_fasta_file(fasta_path)
parsed_gff3 <- read_gff3_file(gff3_path)
parsed_genbank <- read_genbank_file(genbank_path)

stopifnot(parsed_fasta$sequence == record$sequence)
stopifnot(nrow(parsed_gff3) == 2)
stopifnot(parsed_genbank$sequence == record$sequence)
stopifnot(nrow(parsed_genbank$features) == 2)

roundtrip_record <- read_fasta_gff3_to_genbank(fasta_path, gff3_path)
stopifnot(roundtrip_record$sequence == record$sequence)
stopifnot(nrow(roundtrip_record$features) == 2)

converted_gb <- convert_sequence_file(gff3_path, output_format = "genbank", gff3_path = fasta_path)
stopifnot(startsWith(converted_gb, "LOCUS"))
stopifnot(grepl("/gene=\"bla\"", converted_gb, fixed = TRUE))

converted_fasta <- convert_sequence_file(genbank_path, output_format = "fasta")
stopifnot(startsWith(converted_fasta, ">plasmid1"))

converted_gff3 <- convert_sequence_file(genbank_path, output_format = "gff3")
stopifnot(grepl("##gff-version 3", converted_gff3, fixed = TRUE))

cat("All sequence I/O tests passed.\n")
