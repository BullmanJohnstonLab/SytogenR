`%||%` <- function(x, y) {
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

empty_mymotif_records <- function() {
  data.frame(
    rec_seq = character(),
    enz_type = character(),
    meth_base = character(),
    meth_type = character(),
    comp_meth_base = character(),
    comp_meth_type = character(),
    stringsAsFactors = FALSE
  )
}

parse_delimited_motif_table <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2) {
    return(NULL)
  }

  separators <- c("\t", ",", ";", "|")
  for (sep in separators) {
    if (!grepl(sep, text, fixed = TRUE)) {
      next
    }
    parsed <- try(
      utils::read.table(
        text = paste(lines, collapse = "\n"),
        sep = sep,
        header = TRUE,
        stringsAsFactors = FALSE,
        comment.char = "",
        quote = '"',
        fill = TRUE
      ),
      silent = TRUE
    )
    if (!inherits(parsed, "try-error") && ncol(parsed) >= 1) {
      return(parsed)
    }
  }

  NULL
}

normalize_motif_columns <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }

  names(df) <- trimws(names(df))
  normalized <- tolower(names(df))

  motif_aliases <- c("motif", "rec_seq", "recognition_motif", "recognition_sequence", "sequence", "seq")
  motif_idx <- which(normalized %in% motif_aliases)
  if (length(motif_idx) == 0) {
    stop("Could not find a recognition-sequence field in the motif table.")
  }
  names(df)[motif_idx[1]] <- "motif"

  defaults <- c("enz_type", "meth_base", "meth_type", "comp_meth_base", "comp_meth_type")
  current <- tolower(names(df))
  for (column in defaults) {
    if (!column %in% current) {
      df[[column]] <- NA_character_
      current <- c(current, column)
    }
  }

  df
}

parse_rebase_record <- function(chunk) {
  matches <- gregexpr("<([A-Za-z0-9_\\-]+)>([^<>]*)", chunk, perl = TRUE)
  tokens <- regmatches(chunk, matches)[[1]]
  if (length(tokens) == 0) {
    return(NULL)
  }

  fields <- list()
  for (token in tokens) {
    parsed <- regmatches(token, regexec("<([A-Za-z0-9_\\-]+)>([^<>]*)", token, perl = TRUE))[[1]]
    if (length(parsed) < 3) {
      next
    }
    fields[[parsed[2]]] <- parsed[3]
  }

  motif <- fields$rec_seq %||% fields$recognition_sequence %||% fields$sequence %||% fields$seq
  if (!nzchar(trimws(as.character(motif)))) {
    return(NULL)
  }

  list(
    motif = toupper(trimws(as.character(motif))),
    enz_type = trimws(as.character(fields$enz_type %||% NA_character_)),
    meth_base = trimws(as.character(fields$meth_base %||% NA_character_)),
    meth_type = trimws(as.character(fields$meth_type %||% NA_character_)),
    comp_meth_base = trimws(as.character(fields$comp_meth_base %||% NA_character_)),
    comp_meth_type = trimws(as.character(fields$comp_meth_type %||% NA_character_))
  )
}

parse_rebase_motifs <- function(text) {
  chunks <- strsplit(text, "<>", fixed = TRUE)[[1]]
  chunks <- trimws(chunks)
  chunks <- chunks[nzchar(chunks)]
  if (length(chunks) == 0) {
    return(list())
  }

  parsed <- lapply(chunks, parse_rebase_record)
  Filter(Negate(is.null), parsed)
}

parse_motif_text <- function(text) {
  if (is.null(text)) {
    stop("Motif text is empty.")
  }

  text <- paste(as.character(text), collapse = "\n")
  text <- enc2utf8(text)
  text <- trimws(text)
  if (!nzchar(text)) {
    stop("Motif text is empty.")
  }

  guessed <- try(normalize_motif_columns(parse_delimited_motif_table(text)), silent = TRUE)
  if (!inherits(guessed, "try-error") && !is.null(guessed) && nrow(guessed) > 0) {
    guessed$motif <- toupper(trimws(as.character(guessed$motif)))
    guessed <- guessed[nzchar(guessed$motif), , drop = FALSE]
    if (nrow(guessed) > 0) {
      return(guessed)
    }
  }

  rebase_motifs <- parse_rebase_motifs(text)
  if (length(rebase_motifs) == 0) {
    stop("Could not parse the restriction motif table.")
  }

  do.call(
    rbind,
    lapply(rebase_motifs, function(entry) {
      data.frame(
        motif = entry$motif,
        enz_type = entry$enz_type,
        meth_base = entry$meth_base,
        meth_type = entry$meth_type,
        comp_meth_base = entry$comp_meth_base,
        comp_meth_type = entry$comp_meth_type,
        stringsAsFactors = FALSE
      )
    })
  )
}

parse_motif_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }

  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  parse_motif_text(text)
}

mymotif_records <- function(motif_df) {
  if (is.null(motif_df) || nrow(motif_df) == 0) {
    return(empty_mymotif_records())
  }

  aliases <- list(
    rec_seq = c("rec_seq", "motif", "recognition_motif", "recognition_sequence", "sequence", "seq"),
    enz_type = c("enz_type", "type"),
    meth_base = c("meth_base", "methylated_base_plus"),
    meth_type = c("meth_type", "methylated_base_plus_type"),
    comp_meth_base = c("comp_meth_base", "methylated_base_minus"),
    comp_meth_type = c("comp_meth_type", "methylated_base_minus_type")
  )

  columns <- tolower(names(motif_df))

  get_value <- function(field, i) {
    source <- aliases[[field]][aliases[[field]] %in% columns]
    if (length(source) == 0) {
      return("")
    }
    source_name <- names(motif_df)[match(source[1], columns)]
    value <- motif_df[[source_name]][i]
    if (is.na(value)) {
      return("")
    }
    trimws(as.character(value))
  }

  rows <- lapply(seq_len(nrow(motif_df)), function(i) {
    rec_seq <- toupper(get_value("rec_seq", i))
    if (!nzchar(rec_seq)) {
      return(NULL)
    }

    data.frame(
      rec_seq = rec_seq,
      enz_type = get_value("enz_type", i) %||% "-1",
      meth_base = get_value("meth_base", i) %||% "-",
      meth_type = get_value("meth_type", i) %||% "-",
      comp_meth_base = get_value("comp_meth_base", i) %||% "-",
      comp_meth_type = get_value("comp_meth_type", i) %||% "-",
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(empty_mymotif_records())
  }

  do.call(rbind, rows)
}
