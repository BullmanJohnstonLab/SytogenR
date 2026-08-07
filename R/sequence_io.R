# Sequence I/O helpers for FASTA, GFF3, and GenBank.

sio_or <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    return(y)
  }
  if (is.character(x) && !nzchar(trimws(x))) {
    return(y)
  }
  if (all(is.na(x))) {
    return(y)
  }
  x
}

normalize_dna_sequence <- function(sequence) {
  sequence <- toupper(paste(as.character(sequence), collapse = ""))
  sequence <- gsub("[^A-Z]", "", sequence)
  if (!nzchar(sequence)) {
    stop("Sequence input is empty after normalization.")
  }
  sequence
}

normalize_protein_sequence <- function(sequence) {
  sequence <- toupper(paste(as.character(sequence), collapse = ""))
  sequence <- gsub("[^A-Z*]", "", sequence)
  if (!nzchar(sequence)) {
    stop("Protein sequence input is empty after normalization.")
  }
  sequence
}

read_fasta_records_text <- function(text) {
  lines <- strsplit(as.character(text), "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("FASTA text is empty.")
  }

  header_indices <- which(startsWith(lines, ">"))
  if (length(header_indices) == 0) {
    stop("FASTA text is missing a header line.")
  }

  header_indices <- c(header_indices, length(lines) + 1)
  records <- list()
  for (i in seq_len(length(header_indices) - 1)) {
    start_line <- header_indices[i]
    end_line <- header_indices[i + 1] - 1
    header <- sub("^>", "", trimws(lines[start_line]))
    sequence_lines <- lines[(start_line + 1):end_line]
    sequence_lines <- sequence_lines[!startsWith(sequence_lines, ">")]
    sequence <- paste(sequence_lines, collapse = "")
    records[[length(records) + 1]] <- list(
      sequence = sequence,
      id = if (nzchar(header)) strsplit(header, "\\s+", perl = TRUE)[[1]][1] else paste0("sequence_", i),
      name = if (nzchar(header)) header else paste0("sequence_", i),
      description = header,
      header = header,
      format = "fasta"
    )
  }
  records
}

read_fasta_records_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  read_fasta_records_text(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

read_protein_fasta_file <- function(path) {
  records <- read_fasta_records_file(path)
  lapply(records, function(record) {
    record$sequence <- normalize_protein_sequence(record$sequence)
    record$format <- "protein_fasta"
    record
  })
}

sio_standard_genetic_code <- function() {
  c(
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
}

sio_reverse_complement_dna <- function(sequence) {
  sequence <- normalize_dna_sequence(sequence)
  complements <- chartr("ACGT", "TGCA", sequence)
  paste(rev(strsplit(complements, "", fixed = TRUE)[[1]]), collapse = "")
}

sio_translate_dna <- function(sequence, frame = 0) {
  code <- sio_standard_genetic_code()
  sequence <- normalize_dna_sequence(sequence)
  if (frame < 0 || frame > 2) {
    stop("Frame must be 0, 1, or 2.")
  }
  if (nchar(sequence) - frame < 3) {
    return("")
  }
  sequence <- substr(sequence, frame + 1, nchar(sequence))
  codon_starts <- seq(1, nchar(sequence) - 2, by = 3)
  if (length(codon_starts) == 0) {
    return("")
  }
  codons <- substring(sequence, codon_starts, codon_starts + 2)
  amino_acids <- vapply(codons, function(codon) {
    aa <- code[[codon]]
    if (is.null(aa)) {
      "X"
    } else {
      aa
    }
  }, character(1))
  paste(amino_acids, collapse = "")
}

sio_find_overlapping_matches <- function(pattern, text) {
  pattern <- as.character(pattern)
  text <- as.character(text)
  pattern_length <- nchar(pattern)
  text_length <- nchar(text)

  if (!nzchar(pattern) || pattern_length > text_length) {
    return(integer())
  }

  match_positions <- integer()
  max_start <- text_length - pattern_length + 1L
  for (start in seq_len(max_start)) {
    if (substr(text, start, start + pattern_length - 1L) == pattern) {
      match_positions <- c(match_positions, as.integer(start))
    }
  }

  match_positions
}

sio_collect_protein_hits <- function(protein_sequence, dna_sequence) {
  protein_sequence <- normalize_protein_sequence(protein_sequence)
  dna_sequence <- normalize_dna_sequence(dna_sequence)
  sequence_length <- nchar(dna_sequence)
  candidates <- list()

  for (strand in c("+", "-")) {
    strand_sequence <- if (strand == "+") dna_sequence else sio_reverse_complement_dna(dna_sequence)
    for (frame in 0:2) {
      translated <- sio_translate_dna(strand_sequence, frame = frame)
      if (!nzchar(translated)) next
      match_positions <- sio_find_overlapping_matches(protein_sequence, translated)
      if (length(match_positions) == 0) next
      match_lengths <- rep.int(nchar(protein_sequence), length(match_positions))
      for (i in seq_along(match_positions)) {
        aa_start <- match_positions[i] - 1
        aa_length <- match_lengths[i]
        nt_start_in_strand <- frame + (aa_start * 3) + 1
        nt_end_in_strand <- nt_start_in_strand + (aa_length * 3) - 1
        if (strand == "+") {
          start <- nt_start_in_strand
          end <- nt_end_in_strand
        } else {
          start <- sequence_length - nt_end_in_strand + 1
          end <- sequence_length - nt_start_in_strand + 1
        }
        candidates[[length(candidates) + 1]] <- list(
          start = as.integer(start),
          end = as.integer(end),
          strand = strand,
          frame = frame + 1,
          aa_start = aa_start + 1,
          aa_end = aa_start + aa_length,
          aa_length = as.integer(aa_length),
          matched_protein = protein_sequence,
          partial = FALSE,
          multi_hit = FALSE
        )
      }
    }
  }
  candidates
}

sio_find_partial_protein_match <- function(protein_sequence, dna_sequence, min_fraction = 0.5) {
  protein_sequence <- normalize_protein_sequence(protein_sequence)
  protein_length <- nchar(protein_sequence)
  min_length <- max(3L, ceiling(protein_length * min_fraction))
  dna_sequence <- normalize_dna_sequence(dna_sequence)
  best_candidate <- NULL
  best_aa_length <- 0L

  for (sub_len in rev(seq(min_length, protein_length))) {
    for (sub_start in seq_len(protein_length - sub_len + 1)) {
      sub_protein <- substr(protein_sequence, sub_start, sub_start + sub_len - 1)
      hits <- sio_collect_protein_hits(sub_protein, dna_sequence)
      if (length(hits) == 0) next
      if (as.integer(hits[[1]]$aa_length) > best_aa_length) {
        best_aa_length <- as.integer(hits[[1]]$aa_length)
        best_candidate <- hits[[1]]
        best_candidate$partial <- TRUE
        best_candidate$partial_fraction <- sub_len / protein_length
        best_candidate$partial_query_start <- sub_start
        best_candidate$partial_query_end <- sub_start + sub_len - 1L
      }
      break
    }
    if (!is.null(best_candidate) && best_candidate$aa_length >= min_length) break
  }
  best_candidate
}

sio_disambiguate_hits <- function(candidates, policy = c("best", "strict", "all"), gene_id = NULL) {
  policy <- match.arg(policy)
  if (length(candidates) == 0) {
    return(list())
  }
  if (policy == "all" || length(candidates) == 1) {
    return(candidates)
  }

  if (policy == "strict" && length(candidates) > 1) {
    warning(sprintf(
      "Protein '%s': %d hits — strict policy requires a unique match; all hits returned with multi_hit flag.",
      if (!is.null(gene_id)) gene_id else "unknown", length(candidates)
    ))
    for (i in seq_along(candidates)) candidates[[i]]$multi_hit <- TRUE
    return(candidates)
  }

  # policy == "best": prefer forward strand, then earliest start
  fwd <- Filter(function(c) c$strand == "+", candidates)
  pool <- if (length(fwd) > 0) fwd else candidates
  best_idx <- which.min(sapply(pool, function(c) c$start))
  pool[best_idx]
}

sio_find_protein_match <- function(protein_sequence, dna_sequence, topology = "linear",
                                   policy = c("best", "strict", "all"),
                                   allow_partial = FALSE, min_partial_fraction = 0.5) {
  policy <- match.arg(policy)
  all_hits <- sio_collect_protein_hits(protein_sequence, dna_sequence)

  if (length(all_hits) == 0) {
    if (isTRUE(allow_partial)) {
      return(sio_find_partial_protein_match(protein_sequence, dna_sequence, min_fraction = min_partial_fraction))
    }
    return(NULL)
  }

  result <- sio_disambiguate_hits(all_hits, policy = policy, gene_id = NULL)
  if (length(result) == 0) {
    return(NULL)
  }
  if (policy %in% c("all", "strict")) {
    return(result)
  }
  result[[1]]
}

protein_fasta_to_genbank <- function(faa_path, fna_path, output_path = NULL, topology = "linear",
                                     multi_hit_policy = c("best", "strict", "all"),
                                     allow_partial = FALSE, min_partial_fraction = 0.5) {
  multi_hit_policy <- match.arg(multi_hit_policy)
  protein_records <- read_protein_fasta_file(faa_path)
  genome_record <- read_fasta_file(fna_path)
  genome_sequence <- normalize_dna_sequence(genome_record$sequence)

  feature_rows <- list()
  for (i in seq_along(protein_records)) {
    protein_record <- protein_records[[i]]
    matches <- sio_find_protein_match(
      protein_record$sequence, genome_sequence,
      topology = topology,
      policy = multi_hit_policy,
      allow_partial = allow_partial,
      min_partial_fraction = min_partial_fraction
    )
    if (is.null(matches)) next

    match_list <- if (is.list(matches) && !is.null(matches$start)) list(matches) else matches
    for (match in match_list) {
      if (is.null(match$start) || is.null(match$end)) next
      note_text <- if (isTRUE(match$partial)) {
        sprintf("partial_match=%.0f_pct", (match$partial_fraction %||% 1) * 100)
      } else {
        ""
      }
      multi_hit_note <- if (isTRUE(match$multi_hit)) "multi_hit=true" else ""
      extra_notes <- paste(Filter(nzchar, c(note_text, multi_hit_note)), collapse = ";")
      attrs <- sprintf("ID=%s;Name=%s;product=%s", protein_record$id, protein_record$name, protein_record$description)
      if (nzchar(extra_notes)) attrs <- paste0(attrs, ";", extra_notes)
      feature_rows[[length(feature_rows) + 1]] <- data.frame(
        type = ifelse(isTRUE(match$partial), "CDS_partial", "CDS"),
        start = as.integer(match$start),
        end = as.integer(match$end),
        strand = as.character(match$strand),
        gene = protein_record$id,
        attributes = attrs,
        partial = isTRUE(match$partial),
        multi_hit = isTRUE(match$multi_hit),
        stringsAsFactors = FALSE
      )
    }
  }

  features <- if (length(feature_rows) > 0) do.call(rbind, feature_rows) else data.frame(type = character(), start = integer(), end = integer(), strand = character(), gene = character(), attributes = character(), partial = logical(), multi_hit = logical(), stringsAsFactors = FALSE)
  record <- sequence_record(
    sequence = genome_sequence,
    id = genome_record$id,
    name = genome_record$name,
    description = genome_record$description,
    features = features,
    annotations = list(molecule_type = "DNA"),
    topology = topology
  )

  if (!is.null(output_path)) {
    write_genbank_file(record, output_path)
  }
  record
}

wrap_sequence <- function(sequence, width = 60) {
  if (!nzchar(sequence)) {
    return(character())
  }
  starts <- seq(1, nchar(sequence), by = width)
  vapply(starts, function(start) {
    substr(sequence, start, min(start + width - 1, nchar(sequence)))
  }, character(1))
}

sequence_record <- function(sequence, id = "sequence", name = NULL, description = "", features = NULL, annotations = NULL, topology = "linear") {
  if (is.null(name) || !nzchar(as.character(name))) {
    name <- id
  }
  if (is.null(features)) {
    features <- data.frame(type = character(), start = integer(), end = integer(), strand = character(), gene = character(), stringsAsFactors = FALSE)
  }
  if (is.null(annotations)) {
    annotations <- list()
  }
  annotations$molecule_type <- sio_or(annotations$molecule_type, "DNA")
  annotations$topology <- topology
  list(
    sequence = normalize_dna_sequence(sequence),
    id = as.character(id),
    name = as.character(name),
    description = as.character(description),
    features = features,
    annotations = annotations,
    topology = topology,
    format = "sequence_record"
  )
}

parse_fasta_text <- function(text) {
  lines <- strsplit(as.character(text), "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("FASTA text is empty.")
  }

  headers <- which(startsWith(lines, ">"))
  if (length(headers) == 0) {
    stop("FASTA text is missing a header line.")
  }

  header <- sub("^>", "", trimws(lines[headers[1]]))
  sequence <- normalize_dna_sequence(lines[!startsWith(lines, ">")])
  list(
    sequence = sequence,
    id = if (nzchar(header)) strsplit(header, "\\s+", perl = TRUE)[[1]][1] else "sequence",
    name = if (nzchar(header)) header else "sequence",
    description = header,
    features = data.frame(),
    annotations = list(molecule_type = "DNA"),
    topology = "linear",
    format = "fasta"
  )
}

read_fasta_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  parse_fasta_text(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

write_fasta_text <- function(record, width = 60) {
  sequence <- normalize_dna_sequence(record$sequence)
  header <- sio_or(record$id, "sequence")
  if (!is.null(record$description) && nzchar(trimws(as.character(record$description)))) {
    header <- paste(header, trimws(as.character(record$description)))
  }
  paste0(
    ">", header, "\n",
    paste(wrap_sequence(sequence, width = width), collapse = "\n"),
    "\n"
  )
}

write_fasta_file <- function(record, path, width = 60) {
  writeLines(write_fasta_text(record, width = width), path, sep = "")
  invisible(path)
}

parse_gff3_attributes <- function(attribute_text) {
  attribute_text <- trimws(as.character(attribute_text))
  if (!nzchar(attribute_text) || attribute_text == ".") {
    return(list())
  }
  chunks <- strsplit(attribute_text, ";", fixed = TRUE)[[1]]
  attributes <- list()
  for (chunk in chunks) {
    if (!nzchar(trimws(chunk))) {
      next
    }
    if (grepl("=", chunk, fixed = TRUE)) {
      pair <- strsplit(chunk, "=", fixed = TRUE)[[1]]
      key <- trimws(pair[1])
      value <- paste(trimws(pair[-1]), collapse = "=")
      attributes[[key]] <- utils::URLdecode(value)
    } else {
      attributes[[trimws(chunk)]] <- TRUE
    }
  }
  attributes
}

format_gff3_attributes <- function(attributes) {
  if (is.null(attributes) || length(attributes) == 0) {
    return(".")
  }
  pieces <- vapply(names(attributes), function(name) {
    value <- attributes[[name]]
    if (isTRUE(value)) {
      return(name)
    }
    paste0(name, "=", utils::URLencode(as.character(value), reserved = TRUE))
  }, character(1))
  paste(pieces, collapse = ";")
}

normalize_gff3_features <- function(features) {
  if (is.null(features) || !is.data.frame(features) || nrow(features) == 0) {
    return(data.frame(type = character(), start = integer(), end = integer(), strand = character(), gene = character(), attributes = character(), stringsAsFactors = FALSE))
  }

  if (!"attributes" %in% names(features)) {
    features$attributes <- ""
  }
  if (!"gene" %in% names(features)) {
    features$gene <- ""
  }
  if (!"source" %in% names(features)) {
    features$source <- "mymotifr"
  }
  if (!"phase" %in% names(features)) {
    features$phase <- "."
  }
  if (!"score" %in% names(features)) {
    features$score <- "."
  }

  for (i in seq_len(nrow(features))) {
    if (!nzchar(trimws(as.character(features$gene[i])))) {
      attrs <- parse_gff3_attributes(features$attributes[i])
      if (!is.null(attrs$gene) && nzchar(trimws(as.character(attrs$gene)))) {
        features$gene[i] <- as.character(attrs$gene)
      }
    }
  }

  features
}

parse_gff3_text <- function(text) {
  lines <- strsplit(as.character(text), "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) {
    return(data.frame(seqid = character(), source = character(), type = character(), start = integer(), end = integer(), score = character(), strand = character(), phase = character(), attributes = character(), stringsAsFactors = FALSE))
  }

  rows <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 9) {
      stop("Invalid GFF3 line: expected 9 tab-separated columns.")
    }
    data.frame(
      seqid = fields[1],
      source = fields[2],
      type = fields[3],
      start = as.integer(fields[4]),
      end = as.integer(fields[5]),
      score = fields[6],
      strand = fields[7],
      phase = fields[8],
      attributes = fields[9],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

read_gff3_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  parse_gff3_text(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

write_gff3_text <- function(record, include_fasta = FALSE) {
  features <- record$features
  if (is.null(features) || !is.data.frame(features)) {
    features <- data.frame()
  }

  lines <- c("##gff-version 3")
  seqid <- sio_or(record$id, "sequence")
  sequence_length <- nchar(normalize_dna_sequence(record$sequence))
  lines <- c(lines, sprintf("##sequence-region %s 1 %s", seqid, sequence_length))

  if (nrow(features) > 0) {
    for (i in seq_len(nrow(features))) {
      feature <- features[i, , drop = FALSE]
      attributes <- list()
      if ("gene" %in% names(features) && nzchar(trimws(as.character(feature$gene)))) {
        attributes$ID <- feature$gene
        attributes$Name <- feature$gene
      } else if ("id" %in% names(features) && nzchar(trimws(as.character(feature$id)))) {
        attributes$ID <- feature$id
      } else {
        attributes$ID <- paste0(tolower(as.character(feature$type)), "_", i)
      }
      if ("attributes" %in% names(features) && nzchar(trimws(as.character(feature$attributes)))) {
        extra <- parse_gff3_attributes(feature$attributes)
        for (nm in names(extra)) {
          attributes[[nm]] <- extra[[nm]]
        }
      }
      if ("gene" %in% names(features) && nzchar(trimws(as.character(feature$gene)))) {
        attributes$gene <- feature$gene
      }
      lines <- c(lines, paste(
        seqid,
        ifelse("source" %in% names(features), as.character(feature$source), "mymotifr"),
        as.character(feature$type),
        as.integer(feature$start),
        as.integer(feature$end),
        ".",
        ifelse("strand" %in% names(features), as.character(feature$strand), "+"),
        ifelse("phase" %in% names(features), as.character(feature$phase), "."),
        format_gff3_attributes(attributes),
        sep = "\t"
      ))
    }
  }

  if (isTRUE(include_fasta)) {
    lines <- c(lines, "##FASTA", paste0(">", seqid), paste(wrap_sequence(normalize_dna_sequence(record$sequence), 60), collapse = "\n"))
  }

  paste(lines, collapse = "\n")
}

write_gff3_file <- function(record, path, include_fasta = FALSE) {
  writeLines(write_gff3_text(record, include_fasta = include_fasta), path)
  invisible(path)
}

parse_genbank_location <- function(location_text) {
  location_text <- gsub("\\s+", "", as.character(location_text))
  strand <- "+"
  if (startsWith(location_text, "complement(")) {
    strand <- "-"
    location_text <- sub("^complement\\((.*)\\)$", "\\1", location_text)
  }
  if (startsWith(location_text, "join(")) {
    location_text <- sub("^join\\((.*)\\)$", "\\1", location_text)
  }
  parts <- strsplit(location_text, ",", fixed = TRUE)[[1]]
  rows <- lapply(parts, function(part) {
    part <- gsub("[<>]", "", part)
    if (!grepl("\\.\\.", part)) {
      pos <- suppressWarnings(as.integer(part))
      return(data.frame(start = pos, end = pos, strand = strand, stringsAsFactors = FALSE))
    }
    bounds <- strsplit(part, "\\.\\.", perl = TRUE)[[1]]
    start <- suppressWarnings(as.integer(bounds[1]))
    end <- suppressWarnings(as.integer(bounds[2]))
    if (is.na(start) || is.na(end)) {
      return(NULL)
    }
    data.frame(start = start, end = end, strand = strand, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(start = integer(), end = integer(), strand = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

parse_genbank_text <- function(text) {
  lines <- strsplit(as.character(text), "\n", fixed = TRUE)[[1]]
  if (length(lines) == 0) {
    stop("GenBank text is empty.")
  }
  trimmed_lines <- trimws(lines)

  locus <- "sequence"
  sequence_lines <- character()
  feature_rows <- list()
  current_type <- NULL
  current_location <- NULL
  current_qualifiers <- list()
  in_features <- FALSE
  in_origin <- FALSE
  origin_index <- which(startsWith(lines, "ORIGIN"))[1]

  flush_feature <- function() {
    if (is.null(current_type) || is.null(current_location)) {
      return(NULL)
    }
    rows <- parse_genbank_location(current_location)
    if (nrow(rows) == 0) {
      return(NULL)
    }
    gene_name <- NA_character_
    note_value <- NA_character_
    if (length(current_qualifiers) > 0) {
      for (qualifier in current_qualifiers) {
        if (startsWith(qualifier, "/gene=")) {
          gene_name <- sub("^/gene=", "", qualifier)
        } else if (startsWith(qualifier, "/locus_tag=")) {
          gene_name <- sub("^/locus_tag=", "", qualifier)
        } else if (startsWith(qualifier, "/label=")) {
          gene_name <- sub("^/label=", "", qualifier)
        } else if (startsWith(qualifier, "/note=")) {
          note_value <- sub("^/note=", "", qualifier)
        }
      }
    }
    if (is.na(gene_name)) {
      gene_name <- note_value
    }
    for (i in seq_len(nrow(rows))) {
      feature_rows[[length(feature_rows) + 1]] <<- data.frame( # nolint: assignment_linter.
        type = current_type,
        start = as.integer(rows$start[i]),
        end = as.integer(rows$end[i]),
        strand = as.character(rows$strand[i]),
        gene = ifelse(is.na(gene_name), "", gene_name),
        stringsAsFactors = FALSE
      )
    }
    NULL
  }

  for (line in lines) {
    if (startsWith(line, "LOCUS")) {
      tokens <- strsplit(trimws(line), "\\s+", perl = TRUE)[[1]]
      if (length(tokens) >= 2) {
        locus <- tokens[2]
      }
    } else if (startsWith(line, "FEATURES")) {
      in_features <- TRUE
      in_origin <- FALSE
      flush_feature()
      current_type <- NULL
      current_location <- NULL
      current_qualifiers <- list()
    } else if (startsWith(line, "ORIGIN")) {
      in_features <- FALSE
      in_origin <- TRUE
      flush_feature()
      current_type <- NULL
      current_location <- NULL
      current_qualifiers <- list()
    } else if (startsWith(line, "//")) {
      in_origin <- FALSE
    } else if (in_origin) {
      sequence_lines <- c(sequence_lines, gsub("[^A-Za-z]", "", line))
    } else if (in_features) {
      if (grepl("^ {5}[A-Za-z_]+", line, perl = TRUE)) {
        flush_feature()
        current_type <- trimws(substr(line, 6, 20))
        current_location <- trimws(substr(line, 21, nchar(line)))
        current_qualifiers <- list()
      } else if (grepl("^ {21}/", line, perl = TRUE)) {
        current_qualifiers <- c(current_qualifiers, trimws(substr(line, 22, nchar(line))))
      }
    }
  }
  flush_feature()

  origin_line_index <- which(trimmed_lines == "ORIGIN")[1]
  terminator_index <- which(trimmed_lines == "//")[1]
  if (!is.na(origin_line_index)) {
    seq_start <- origin_line_index + 1
    seq_end <- if (!is.na(terminator_index) && terminator_index > origin_line_index) terminator_index - 1 else length(lines)
    if (seq_start <= seq_end) {
      sequence_lines <- gsub("[^A-Za-z]", "", lines[seq_start:seq_end])
    }
  }

  if (length(sequence_lines) == 0 && !is.na(origin_index)) {
    tail_lines <- lines[(origin_index + 1):length(lines)]
    tail_lines <- tail_lines[!startsWith(tail_lines, "//")]
    sequence_lines <- gsub("[^A-Za-z]", "", tail_lines)
  }

  if (length(sequence_lines) == 0) {
    origin_match <- regexpr("ORIGIN[\\s\\S]*?//", text, perl = TRUE)
    if (origin_match[1] != -1) {
      origin_block <- regmatches(text, origin_match)
      origin_block <- sub("^.*?ORIGIN\\s*", "", origin_block, perl = TRUE)
      origin_block <- sub("\\/\\/.*$", "", origin_block, perl = TRUE)
      tail_lines <- strsplit(origin_block, "\n", fixed = TRUE)[[1]]
      sequence_lines <- gsub("[^A-Za-z]", "", tail_lines)
    }
  }

  sequence <- normalize_dna_sequence(sequence_lines)
  features <- if (length(feature_rows) > 0) do.call(rbind, feature_rows) else data.frame(type = character(), start = integer(), end = integer(), strand = character(), gene = character(), stringsAsFactors = FALSE)
  list(
    sequence = sequence,
    id = locus,
    name = locus,
    description = locus,
    features = features,
    annotations = list(molecule_type = "DNA"),
    topology = "linear",
    format = "genbank"
  )
}

read_genbank_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  parse_genbank_text(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

write_genbank_text <- function(record, width = 60) {
  sequence <- normalize_dna_sequence(record$sequence)
  locus <- sio_or(record$id, "sequence")
  topology <- tolower(as.character(sio_or(record$topology, sio_or(record$annotations$topology, "linear"))))
  molecule_type <- sio_or(record$annotations$molecule_type, "DNA")
  length_value <- nchar(sequence)

  lines <- c(
    sprintf("LOCUS       %-16s %11d bp    %s    %s", locus, length_value, molecule_type, topology),
    sprintf("DEFINITION  %s.", sio_or(record$description, locus)),
    sprintf("ACCESSION   %s", locus),
    "VERSION     ",
    "KEYWORDS    .",
    "SOURCE      synthetic DNA construct",
    "  ORGANISM  synthetic construct",
    "FEATURES             Location/Qualifiers"
  )

  features <- record$features
  if (is.null(features) || !is.data.frame(features)) {
    features <- data.frame()
  }
  if (nrow(features) > 0) {
    for (i in seq_len(nrow(features))) {
      feature <- features[i, , drop = FALSE]
      strand <- ifelse("strand" %in% names(features), as.character(feature$strand), "+")
      location <- sprintf("%s..%s", as.integer(feature$start), as.integer(feature$end))
      if (strand == "-") {
        location <- paste0("complement(", location, ")")
      }
      lines <- c(lines, sprintf("     %-16s %s", as.character(feature$type), location))
      if ("gene" %in% names(features) && nzchar(trimws(as.character(feature$gene)))) {
        lines <- c(lines, sprintf("                     /gene=\"%s\"", as.character(feature$gene)))
      }
      if ("attributes" %in% names(features) && nzchar(trimws(as.character(feature$attributes)))) {
        attrs <- parse_gff3_attributes(feature$attributes)
        for (attr_name in names(attrs)) {
          if (identical(attr_name, "gene")) {
            next
          }
          lines <- c(lines, sprintf("                     /%s=\"%s\"", attr_name, as.character(attrs[[attr_name]])))
        }
      }
    }
  }

  lines <- c(lines, "ORIGIN")
  seq_chunks <- wrap_sequence(tolower(sequence), width = 60)
  if (length(seq_chunks) > 0) {
    positions <- seq(1, length.out = length(seq_chunks), by = 60)
    for (i in seq_along(seq_chunks)) {
      block <- seq_chunks[i]
      line_number <- sprintf("%9d", positions[i])
      chunk_starts <- seq(1, nchar(block), by = 10)
      grouped <- paste(vapply(chunk_starts, function(start) {
        substr(block, start, min(start + 9, nchar(block)))
      }, character(1)), collapse = " ")
      lines <- c(lines, paste(line_number, grouped))
    }
  }
  lines <- c(lines, "//")
  paste(lines, collapse = "\n")
}

write_genbank_file <- function(record, path, width = 60) {
  writeLines(write_genbank_text(record, width = width), path)
  invisible(path)
}

sequence_file_to_record <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("fa", "fasta", "fna")) {
    return(read_fasta_file(path))
  }
  if (ext %in% c("gb", "gbk", "genbank")) {
    return(read_genbank_file(path))
  }
  if (ext %in% c("gff", "gff3")) {
    return(list(sequence = "", id = basename(path), features = read_gff3_file(path), annotations = list(molecule_type = "DNA"), topology = "linear", format = "gff3"))
  }
  stop(sprintf("Unsupported file format: %s", path))
}

read_fasta_gff3_to_genbank <- function(fasta_path, gff3_path, output_path = NULL, topology = "linear") {
  fasta_record <- read_fasta_file(fasta_path)
  gff3_features <- normalize_gff3_features(read_gff3_file(gff3_path))
  record <- sequence_record(
    sequence = fasta_record$sequence,
    id = fasta_record$id,
    name = fasta_record$name,
    description = fasta_record$description,
    features = gff3_features,
    annotations = list(molecule_type = "DNA"),
    topology = topology
  )
  if (!is.null(output_path)) {
    write_genbank_file(record, output_path)
  }
  record
}

record_to_fasta <- function(record, output_path = NULL) {
  fasta_text <- write_fasta_text(record)
  if (!is.null(output_path)) {
    writeLines(fasta_text, output_path)
    return(invisible(output_path))
  }
  fasta_text
}

record_to_gff3 <- function(record, output_path = NULL, include_fasta = FALSE) {
  gff3_text <- write_gff3_text(record, include_fasta = include_fasta)
  if (!is.null(output_path)) {
    writeLines(gff3_text, output_path)
    return(invisible(output_path))
  }
  gff3_text
}

record_to_genbank <- function(record, output_path = NULL) {
  genbank_text <- write_genbank_text(record)
  if (!is.null(output_path)) {
    writeLines(genbank_text, output_path)
    return(invisible(output_path))
  }
  genbank_text
}

convert_sequence_file <- function(input_path, output_format = c("genbank", "fasta", "gff3"), output_path = NULL, gff3_path = NULL, topology = "linear", include_fasta = FALSE) {
  output_format <- match.arg(output_format)
  input_ext <- tolower(tools::file_ext(input_path))

  if (input_ext %in% c("faa", "pep", "protein")) {
    if (output_format != "genbank") {
      stop("Protein FASTA input currently supports conversion to GenBank only.")
    }
    if (is.null(gff3_path)) {
      stop("A nucleotide FASTA file is required when converting a protein FASTA to GenBank.")
    }
    return(protein_fasta_to_genbank(input_path, gff3_path, output_path = output_path, topology = topology))
  }

  record <- sequence_file_to_record(input_path)
  if (record$format == "gff3") {
    if (is.null(gff3_path)) {
      stop("A FASTA file is required to convert GFF3 input into GenBank/FASTA.")
    }
    fasta_record <- read_fasta_file(gff3_path)
    record <- sequence_record(fasta_record$sequence, id = fasta_record$id, name = fasta_record$name, description = fasta_record$description, features = normalize_gff3_features(record$features), annotations = list(molecule_type = "DNA"), topology = topology)
  }

  if (output_format == "fasta") {
    return(record_to_fasta(record, output_path = output_path))
  }
  if (output_format == "gff3") {
    return(record_to_gff3(record, output_path = output_path, include_fasta = include_fasta))
  }
  record_to_genbank(record, output_path = output_path)
}

read_faa_fna_to_genbank <- function(faa_path, fna_path, output_path = NULL, topology = "linear") {
  protein_fasta_to_genbank(faa_path, fna_path, output_path = output_path, topology = topology)
}
