# Codon Bias
#
# Take one or more DNA sequences or codons and summarize the codon usage bias.
# The output can be used to guide choices for motif-destroying mutations by
# preferring codons that are common in the construct.

standard_genetic_code <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W",
  CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R",
  GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

normalize_codon_input <- function(sequence) {
  if (is.null(sequence)) {
    stop("Sequence input is empty.")
  }

  sequence <- as.character(sequence)
  sequence <- sequence[!is.na(sequence)]
  if (length(sequence) == 0) {
    stop("Sequence input is empty.")
  }

  if (length(sequence) == 1 && nchar(sequence[1]) == 3) {
    codons <- toupper(sequence[1])
    return(codons)
  }

  if (all(nchar(sequence) == 3)) {
    return(toupper(sequence))
  }

  flattened <- paste(sequence, collapse = "")
  flattened <- toupper(gsub("[^A-Z]", "", flattened))
  if (!nzchar(flattened)) {
    stop("Sequence input is empty after normalization.")
  }

  if (nchar(flattened) %% 3 != 0) {
    stop("DNA sequence length must be a multiple of 3.")
  }

  substring(flattened, seq(1, nchar(flattened) - 2, 3), seq(3, nchar(flattened), 3))
}

translate_dna_sequence <- function(sequence) {
  codons <- normalize_codon_input(sequence)
  amino_acids <- vapply(codons, function(codon) {
    if (!codon %in% names(standard_genetic_code)) {
      stop(sprintf("Unsupported codon: %s", codon))
    }
    standard_genetic_code[[codon]]
  }, character(1))

  paste(amino_acids, collapse = "")
}

calculate_codon_bias <- function(sequence) {
  codons <- normalize_codon_input(sequence)
  invalid_codons <- codons[!codons %in% names(standard_genetic_code)]
  if (length(invalid_codons) > 0) {
    stop(sprintf("Unsupported codon(s): %s", paste(unique(invalid_codons), collapse = ", ")))
  }

  codon_counts <- table(codons)
  codon_counts <- codon_counts[order(names(codon_counts))]

  amino_acids <- vapply(names(codon_counts), function(codon) {
    standard_genetic_code[[codon]]
  }, character(1))

  amino_acid_counts <- table(amino_acids)
  amino_acid_counts <- amino_acid_counts[order(names(amino_acid_counts))]

  bias_table <- data.frame(
    codon = names(codon_counts),
    aa = amino_acids,
    count = as.integer(codon_counts),
    stringsAsFactors = FALSE
  )

  bias_table$relative_synonymous_usage <- vapply(seq_len(nrow(bias_table)), function(i) {
    aa <- bias_table$aa[i]
    aa_total <- sum(bias_table$count[bias_table$aa == aa], na.rm = TRUE)
    if (aa_total == 0) {
      return(0)
    }
    bias_table$count[i] / aa_total
  }, numeric(1))

  preferred_codons <- vapply(split(bias_table, bias_table$aa), function(group) {
    group$codon[which.max(group$count)]
  }, character(1))

  list(
    codon_counts = stats::setNames(as.integer(codon_counts), names(codon_counts)),
    amino_acid_counts = stats::setNames(as.integer(amino_acid_counts), names(amino_acid_counts)),
    codon_bias = bias_table,
    preferred_codons = preferred_codons
  )
}

codon_bias_from_fasta <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("FASTA file is empty.")
  }

  sequences <- list()
  current_name <- NA_character_
  current_sequence <- character()

  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.na(current_name) && length(current_sequence) > 0) {
        sequences[[current_name]] <- paste(current_sequence, collapse = "")
      }
      current_name <- sub("^>", "", trimws(line))
      current_sequence <- character()
    } else {
      current_sequence <- c(current_sequence, toupper(gsub("[^A-Z]", "", line)))
    }
  }

  if (!is.na(current_name) && length(current_sequence) > 0) {
    sequences[[current_name]] <- paste(current_sequence, collapse = "")
  }

  if (length(sequences) == 0) {
    stop("No FASTA entries were parsed.")
  }

  lapply(sequences, calculate_codon_bias)
}
