# Motif finder
#
# Using a user-provided DNA construct, locate motifs in the forward and reverse
# complement orientations and return their genomic coordinates for downstream
# mutation or annotation workflows.

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

find_motifs <- function(sequence, motifs, include_reverse = TRUE) {
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
    search_results[[1]] <- .find_motif_matches(sequence, motif, "+")
    if (include_reverse) {
      search_results[[2]] <- .find_motif_matches(sequence, motif, "-")
    }

    do.call(rbind, Filter(function(x) !is.null(x), search_results))
  })

  rows <- Filter(function(x) !is.null(x) && nrow(x) > 0, rows)
  if (length(rows) == 0) {
    return(data.frame(motif = character(), start = integer(), end = integer(), strand = character(), match = character(), stringsAsFactors = FALSE))
  }

  do.call(rbind, rows)
}

.find_motif_matches <- function(sequence, motif, strand) {
  if (strand == "-") {
    search_sequence <- reverse_complement(sequence)
    start_offset <- 0
  } else {
    search_sequence <- sequence
    start_offset <- 0
  }

  match_positions <- gregexpr(motif, search_sequence, ignore.case = TRUE, perl = TRUE)[[1]]
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
    if (strand == "-") {
      start_in_original <- nchar(sequence) - end + 1
      end_in_original <- nchar(sequence) - start + 1
      strand_value <- "-"
    } else {
      start_in_original <- start
      end_in_original <- end
      strand_value <- "+"
    }

    data.frame(
      motif = motif,
      start = start_in_original,
      end = end_in_original,
      strand = strand_value,
      match = substr(sequence, start_in_original, end_in_original),
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
