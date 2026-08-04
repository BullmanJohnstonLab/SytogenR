source(file.path("R", "parse_motif.R"))

plain_text <- "motif\tenz_type\nATGC\t2\nGATC\t4\n"
plain_df <- parse_motif_text(plain_text)
stopifnot(all(c("motif", "enz_type") %in% names(plain_df)))
stopifnot(nrow(plain_df) == 2)
stopifnot(plain_df$motif[1] == "ATGC")
stopifnot(plain_df$enz_type[1] == "2")

rebase_text <- "<enz_type>2<rec_seq>ATGC<meth_base>C<><enz_type>4<rec_seq>GATC<>"
rebase_df <- parse_motif_text(rebase_text)
stopifnot(nrow(rebase_df) == 2)
stopifnot(rebase_df$motif[1] == "ATGC")
stopifnot(rebase_df$motif[2] == "GATC")

records <- mymotif_records(data.frame(motif = c("atgc", ""), enz_type = c("2", ""), stringsAsFactors = FALSE))
stopifnot(nrow(records) == 1)
stopifnot(records$rec_seq[1] == "ATGC")
stopifnot(records$enz_type[1] == "2")

cat("All MyMotif parser tests passed.\n")
