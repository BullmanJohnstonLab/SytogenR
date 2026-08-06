# Motif finder
#
# Using a user-provided DNA construct, locate motifs in the forward and reverse
# complement orientations and return their genomic coordinates for downstream
# mutation or annotation workflows.

iupac_complements <- c(
  A = "T", C = "G", G = "C", T = "A",
  R = "Y", Y = "R", S = "S", W = "W",
  K = "M", M = "K", B = "V", V = "B",
  D = "H", H = "D", N = "N"
)

iupac_regex_values <- c(
  A = "A", C = "C", G = "G", T = "T",
  R = "[AG]", Y = "[CT]", S = "[GC]", W = "[AT]",
  K = "[GT]", M = "[AC]", B = "[CGT]", D = "[AGT]",
  H = "[ACT]", V = "[ACG]", N = "[ACGT]"
)

extract_sequence_interval <- function(sequence, start, end, circular = FALSE) {
  sequence <- as.character(sequence)
  sequence_length <- nchar(sequence)
  if (sequence_length == 0) {
    return("")
  }

  if (!circular) {
    start <- max(1, start)
    end <- min(sequence_length, end)
    if (start > end) {
      return("")
    }
    return(substr(sequence, start, end))
  }

  start <- ((start - 1) %% sequence_length) + 1
  end <- ((end - 1) %% sequence_length) + 1
  if (start <= end) {
    return(substr(sequence, start, end))
  }

  paste0(substr(sequence, start, sequence_length), substr(sequence, 1, end))
}

iupac_to_regex <- function(pattern) {
  pattern <- toupper(as.character(pattern))
  paste(vapply(strsplit(pattern, "", fixed = TRUE)[[1]], function(base) {
    iupac_regex_values[[base]] %||% base
  }, character(1)), collapse = "")
}

reverse_complement_iupac <- function(sequence) {
  sequence <- toupper(as.character(sequence))
  chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  complements <- vapply(chars, function(base) {
    iupac_complements[[base]] %||% "N"
  }, character(1))
  paste(rev(complements), collapse = "")
}

reverse_complement <- function(sequence) {
  if (is.null(sequence)) {
    stop("Sequence input is empty.")
  }

  sequence <- toupper(as.character(sequence))
  sequence <- gsub("[^A-Z]", "", sequence)
  if (!nzchar(sequence)) {
    stop("Sequence input is empty after normalization.")
  }

  map <- c(A = "T", T = "A", C = "G", G = "C")
  chars <- strsplit(sequence, "")[[1]]
  complement <- vapply(chars, function(base) {
    if (base %in% names(map)) {
      map[[base]]
    } else {
      base
    }
  }, character(1))

  paste(rev(complement), collapse = "")
}

find_motifs <- function(sequence, motifs, include_reverse = TRUE, topology = c("linear", "circular")) {
  topology <- match.arg(topology)
  if (is.null(sequence) || is.na(sequence) || !nzchar(trimws(as.character(sequence)))) {
    stop("Sequence input is empty.")
  }

  if (is.null(motifs) || length(motifs) == 0) {
    stop("At least one motif is required.")
  }

  sequence <- toupper(as.character(sequence))
  sequence <- gsub("[^A-Z]", "", sequence)
  if (!nzchar(sequence)) {
    stop("Sequence input is empty after normalization.")
  }

  motif_list <- unique(as.character(motifs))
  motif_list <- motif_list[nzchar(trimws(motif_list))]
  if (length(motif_list) == 0) {
    stop("At least one motif is required.")
  }

  rows <- lapply(motif_list, function(motif) {
    motif <- toupper(as.character(motif))
    motif <- gsub("[^A-Z]", "", motif)
    if (!nzchar(motif)) {
      return(NULL)
    }

    search_results <- list()
    search_results[[1]] <- .find_motif_matches(sequence, motif, "+", topology = topology)
    if (include_reverse) {
      search_results[[2]] <- .find_motif_matches(sequence, motif, "-", topology = topology)
    }

    do.call(rbind, Filter(function(x) !is.null(x), search_results))
  })

  rows <- Filter(function(x) !is.null(x) && nrow(x) > 0, rows)
  if (length(rows) == 0) {
    return(data.frame(motif = character(), start = integer(), end = integer(), strand = character(), match = character(), stringsAsFactors = FALSE))
  }

  do.call(rbind, rows)
}

.find_motif_matches <- function(sequence, motif, strand, topology = c("linear", "circular")) {
  topology <- match.arg(topology)
  search_pattern <- if (strand == "-") reverse_complement_iupac(motif) else motif
  regex <- iupac_to_regex(search_pattern)

  search_sequence <- toupper(as.character(sequence))
  sequence_length <- nchar(search_sequence)
  if (topology == "circular") {
    search_sequence <- paste0(search_sequence, substr(search_sequence, 1, max(0, nchar(motif) - 1)))
  }

  match_positions <- gregexpr(regex, search_sequence, ignore.case = TRUE, perl = TRUE)[[1]]
  if (length(match_positions) == 1 && match_positions[1] == -1) {
    return(NULL)
  }

  matches <- as.integer(match_positions)
  lengths <- attr(match_positions, "match.length")
  if (length(matches) == 0) {
    return(NULL)
  }

  entries <- lapply(seq_along(matches), function(i) {
    start <- matches[i]
    end <- start + lengths[i] - 1
    if (topology == "circular") {
      if (start > sequence_length) {
        return(NULL)
      }
      if (strand == "-") {
        start_in_original <- sequence_length - end + 1
        end_in_original <- sequence_length - start + 1
      } else {
        start_in_original <- start
        end_in_original <- end
      }
    } else if (strand == "-") {
      start_in_original <- sequence_length - end + 1
      end_in_original <- sequence_length - start + 1
    } else {
      start_in_original <- start
      end_in_original <- end
    }

    strand_value <- strand
    match_sequence <- extract_sequence_interval(sequence, start_in_original, end_in_original, circular = topology == "circular")

    data.frame(
      motif = motif,
      start = start_in_original,
      end = end_in_original,
      strand = strand_value,
      match = match_sequence,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, entries)
}

find_motifs_in_fasta <- function(path, motifs, include_reverse = TRUE) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("FASTA file is empty.")
  }

  entries <- list()
  current_name <- NA_character_
  current_sequence <- character()

  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.na(current_name) && length(current_sequence) > 0) {
        entries[[current_name]] <- find_motifs(paste(current_sequence, collapse = ""), motifs, include_reverse = include_reverse)
      }
      current_name <- sub("^>", "", trimws(line))
      current_sequence <- character()
    } else {
      current_sequence <- c(current_sequence, toupper(gsub("[^A-Z]", "", line)))
    }
  }

  if (!is.na(current_name) && length(current_sequence) > 0) {
    entries[[current_name]] <- find_motifs(paste(current_sequence, collapse = ""), motifs, include_reverse = include_reverse)
  }

  entries
}
