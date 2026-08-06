# SyToGen
#
# This is the main program for these tools. This will pull in motifs from parse_motif,
# and DNA annotated with motif_finder. 
#
# Sytogen will find the optimal mutation of a motif so that the motif is destroyed.
# Preferring wobble bases of codons (using the codon_bias of the DNA contruct as a 
# guide). In regions that don't fall inside a protein (as guided by genbank) the base
# change should prefer a GC/AT-preserving switch. All potential changes should be scored
# and the best scoring mutation should be chosen.
# 
# Output includes a decision matrix and the genbank, gff3/fasta
#

sytogen_parse_motif_text <- function(...) {
	get("parse_motif_text", mode = "function", inherits = TRUE)(...)
}

sytogen_calculate_codon_bias <- function(...) {
	get("calculate_codon_bias", mode = "function", inherits = TRUE)(...)
}

sytogen_find_motifs <- function(...) {
	get("find_motifs", mode = "function", inherits = TRUE)(...)
}

sytogen_extract_sequence_interval <- function(...) {
	get("extract_sequence_interval", mode = "function", inherits = TRUE)(...)
}

sytogen_reverse_complement_iupac <- function(...) {
	get("reverse_complement_iupac", mode = "function", inherits = TRUE)(...)
}

sytogen_iupac_to_regex <- function(...) {
	get("iupac_to_regex", mode = "function", inherits = TRUE)(...)
}

sytogen_standard_genetic_code <- function() {
	get("standard_genetic_code", mode = "any", inherits = TRUE)
}

sytogen_clean_sequence <- function(sequence) {
	if (is.null(sequence) || is.na(sequence) || !nzchar(trimws(as.character(sequence)))) {
		stop("Sequence input is empty.")
	}

	cleaned <- toupper(gsub("[^A-Z]", "", as.character(sequence)))
	if (!nzchar(cleaned)) {
		stop("Sequence input is empty after normalization.")
	}

	cleaned
}

sytogen_extract_motifs <- function(motif_text = NULL, motifs = NULL) {
	if (!is.null(motifs) && length(motifs) > 0) {
		parsed <- unique(toupper(trimws(as.character(motifs))))
		parsed <- parsed[nzchar(parsed)]
		if (length(parsed) > 0) {
			return(parsed)
		}
	}

	if (is.null(motif_text) || !nzchar(trimws(as.character(motif_text)))) {
		stop("Provide either motif_text or a motifs vector.")
	}

	parsed_table <- try(sytogen_parse_motif_text(motif_text), silent = TRUE)
	if (!inherits(parsed_table, "try-error") && nrow(parsed_table) > 0) {
		return(unique(toupper(trimws(as.character(parsed_table$motif)))))
	}

	lines <- strsplit(as.character(motif_text), "\n", fixed = FALSE)[[1]]
	lines <- toupper(trimws(lines))
	lines <- lines[nzchar(lines)]
	lines <- lines[lines != "MOTIF"]
	lines <- lines[!startsWith(lines, "#")]
	lines <- unique(lines)
	if (length(lines) == 0) {
		stop("Could not parse motifs from motif_text.")
	}

	lines
}

mutate_base_at <- function(sequence, position, new_base) {
	paste0(
		if (position > 1) substr(sequence, 1, position - 1) else "",
		new_base,
		if (position < nchar(sequence)) substr(sequence, position + 1, nchar(sequence)) else ""
	)
}

base_class <- function(base) {
	if (base %in% c("G", "C")) {
		return("GC")
	}
	if (base %in% c("A", "T")) {
		return("AT")
	}
	"OTHER"
}

sytogen_is_empty <- function(value) {
	is.null(value) || length(value) == 0 || all(is.na(value)) || !nzchar(trimws(as.character(value[1])))
}

sytogen_normalize_ranges <- function(ranges, sequence_length) {
	if (is.null(ranges) || (is.character(ranges) && !nzchar(trimws(ranges)))) {
		return(data.frame(start = integer(), end = integer(), stringsAsFactors = FALSE))
	}

	if (is.data.frame(ranges)) {
		column_names <- tolower(names(ranges))
		if (!all(c("start", "end") %in% column_names)) {
			stop("Range data frame must have start and end columns.")
		}
		start_col <- names(ranges)[match("start", column_names)]
		end_col <- names(ranges)[match("end", column_names)]
		result <- data.frame(
			start = as.integer(ranges[[start_col]]),
			end = as.integer(ranges[[end_col]]),
			stringsAsFactors = FALSE
		)
		return(result[!is.na(result$start) & !is.na(result$end), , drop = FALSE])
	}

	if (is.list(ranges) && length(ranges) > 0 && is.numeric(ranges[[1]]) && length(ranges[[1]]) >= 2) {
		result <- do.call(rbind, lapply(ranges, function(range_item) {
			data.frame(start = as.integer(range_item[[1]]), end = as.integer(range_item[[2]]), stringsAsFactors = FALSE)
		}))
		return(result)
	}

	range_text <- paste(as.character(ranges), collapse = ",")
	parts <- strsplit(range_text, ",", fixed = TRUE)[[1]]
	rows <- lapply(parts, function(part) {
		part <- trimws(part)
		if (!nzchar(part)) {
			return(NULL)
		}
		if (!grepl("-", part, fixed = TRUE)) {
			stop(sprintf("Could not parse range '%s' — expected a format like '100-200'.", part))
		}
		bounds <- strsplit(part, "-", fixed = TRUE)[[1]]
		if (length(bounds) < 2) {
			stop(sprintf("Could not parse range '%s' — expected a format like '100-200'.", part))
		}
		start_value <- suppressWarnings(as.integer(trimws(bounds[1])))
		end_value <- suppressWarnings(as.integer(trimws(bounds[length(bounds)])))
		if (is.na(start_value) || is.na(end_value)) {
			stop(sprintf("Could not parse range '%s' — start and end must be whole numbers.", part))
		}
		if (start_value < 1 || end_value < 1) {
			stop(sprintf("Range '%s' must use positive, 1-based positions.", part))
		}
		if (start_value > end_value) {
			stop(sprintf("Range '%s' has a start position after its end position.", part))
		}
		if (end_value > sequence_length) {
			stop(sprintf("Range '%s' extends past the end of the sequence (length %s).", part, sequence_length))
		}
		data.frame(start = start_value, end = end_value, stringsAsFactors = FALSE)
	})
	rows <- Filter(Negate(is.null), rows)
	if (length(rows) == 0) {
		return(data.frame(start = integer(), end = integer(), stringsAsFactors = FALSE))
	}
	do.call(rbind, rows)
}

sytogen_position_in_ranges <- function(position, ranges) {
	if (is.null(ranges) || nrow(ranges) == 0) {
		return(FALSE)
	}
	any(position >= ranges$start & position <= ranges$end)
}

sytogen_parse_codon_usage <- function(codon_df, fallback_sequence = NULL) {
	if (is.null(codon_df) || nrow(codon_df) == 0) {
		if (is.null(fallback_sequence)) {
			return(numeric())
		}
		bias <- sytogen_calculate_codon_bias(fallback_sequence)
		return(stats::setNames(as.numeric(bias$codon_bias$relative_synonymous_usage), bias$codon_bias$codon))
	}

	column_names <- tolower(names(codon_df))
	codon_col <- names(codon_df)[match("codon", column_names)]
	if (is.na(codon_col)) {
		return(numeric())
	}
	score_candidates <- c("fraction", "frequency", "value", "usage", "proportion", "ranking_ratio", "ranking", "count")
	score_col <- names(codon_df)[match(score_candidates, column_names)]
	score_col <- score_col[!is.na(score_col)][1]
	if (is.na(score_col)) {
		return(numeric())
	}
	invert <- tolower(score_col) %in% c("ranking", "ranking_ratio")
	usage <- numeric()
	for (i in seq_len(nrow(codon_df))) {
		codon <- toupper(trimws(as.character(codon_df[[codon_col]][i])))
		value <- suppressWarnings(as.numeric(codon_df[[score_col]][i]))
		if (!nzchar(codon) || is.na(value)) {
			next
		}
		usage[[codon]] <- if (invert) -value else value
	}
	usage
}

sytogen_parse_motif_table <- function(motif_df, sequence, topology = "circular") {
	if (is.null(motif_df) || nrow(motif_df) == 0) {
		return(data.frame(motif = character(), start = integer(), end = integer(), strand = character(), enz_type = character(), stringsAsFactors = FALSE))
	}

	column_names <- tolower(names(motif_df))
	if (!"motif" %in% column_names) {
		stop("Motif table must include a motif column.")
	}
	motif_col <- names(motif_df)[match("motif", column_names)]
	strand_col <- names(motif_df)[match("strand", column_names)]
	start_col <- names(motif_df)[match("start", column_names)]
	end_col <- names(motif_df)[match("end", column_names)]
	enz_col <- names(motif_df)[match("enz_type", column_names)]
	if (is.na(enz_col)) {
		enz_col <- names(motif_df)[match("type", column_names)]
	}

	rows <- list()
	seen <- new.env(parent = emptyenv())
	for (i in seq_len(nrow(motif_df))) {
		motif <- toupper(trimws(as.character(motif_df[[motif_col]][i])))
		if (!nzchar(motif)) {
			next
		}
		strand <- if (!is.na(strand_col)) as.character(motif_df[[strand_col]][i]) else "+"
		strand <- ifelse(nzchar(trimws(strand)), strand, "+")
		enz_type <- if (!is.na(enz_col)) as.character(motif_df[[enz_col]][i]) else ""

		if (!is.na(start_col) && !is.na(end_col) && !is.na(motif_df[[start_col]][i]) && !is.na(motif_df[[end_col]][i])) {
			start <- as.integer(motif_df[[start_col]][i])
			end <- as.integer(motif_df[[end_col]][i])
			key <- paste(motif, start, end, strand, sep = "|")
			if (is.null(seen[[key]])) {
				seen[[key]] <- TRUE
				rows[[length(rows) + 1]] <- data.frame(motif = motif, start = start, end = end, strand = strand, enz_type = enz_type, stringsAsFactors = FALSE)
			}
		} else {
			hits <- sytogen_find_motifs(sequence, c(motif), include_reverse = TRUE, topology = topology)
			if (nrow(hits) == 0) {
				next
			}
			for (j in seq_len(nrow(hits))) {
				key <- paste(hits$motif[j], hits$start[j], hits$end[j], hits$strand[j], sep = "|")
				if (is.null(seen[[key]])) {
					seen[[key]] <- TRUE
					rows[[length(rows) + 1]] <- data.frame(motif = hits$motif[j], start = hits$start[j], end = hits$end[j], strand = hits$strand[j], enz_type = enz_type, stringsAsFactors = FALSE)
				}
			}
		}
	}

	if (length(rows) == 0) {
		return(data.frame(motif = character(), start = integer(), end = integer(), strand = character(), enz_type = character(), stringsAsFactors = FALSE))
	}
	do.call(rbind, rows)
}

sytogen_is_type_iv_motif <- function(motif_row) {
	enz_type <- tolower(as.character(motif_row$enz_type %||% motif_row$type %||% ""))
	grepl("(^|[^0-9])4([^0-9]|$)", enz_type) || grepl("type\\s*iv|type\\s*4", enz_type)
}

sytogen_exact_motif_destroyed <- function(original_sequence, mutated_sequence, hit) {
	original_site <- sytogen_extract_sequence_interval(original_sequence, hit$start, hit$end, circular = FALSE)
	mutated_site <- sytogen_extract_sequence_interval(mutated_sequence, hit$start, hit$end, circular = FALSE)
	motif <- toupper(as.character(hit$motif))
	motif_rc <- sytogen_reverse_complement_iupac(motif)
	original_matches <- grepl(sytogen_iupac_to_regex(motif), original_site, perl = TRUE) || grepl(sytogen_iupac_to_regex(motif_rc), original_site, perl = TRUE)
	mutated_matches <- grepl(sytogen_iupac_to_regex(motif), mutated_site, perl = TRUE) || grepl(sytogen_iupac_to_regex(motif_rc), mutated_site, perl = TRUE)
	original_matches && !mutated_matches
}

sytogen_candidate_hits_created <- function(original_sequence, mutated_sequence, motifs, topology = "linear") {
	created <- character()
	for (motif in unique(toupper(motifs))) {
		original_hits <- sytogen_find_motifs(original_sequence, c(motif), include_reverse = TRUE, topology = topology)
		mutated_hits <- sytogen_find_motifs(mutated_sequence, c(motif), include_reverse = TRUE, topology = topology)
		if (nrow(mutated_hits) > nrow(original_hits)) {
			created <- c(created, motif)
		}
	}
	unique(created)
}

sytogen_collect_decisions <- function(sequence, motif_hits, motif_df, cds_ranges = NULL, protected_ranges = NULL, codon_usage = numeric(), preserve_gc = FALSE, topology = "linear") {
	if (is.null(motif_hits) || nrow(motif_hits) == 0) {
		return(list(decision_matrix = data.frame(), applied_mutations = list(), mutated_sequence = sequence, new_motifs = character(), resolved_motif_keys = character()))
	}

	if (is.null(cds_ranges)) {
		cds_ranges <- data.frame(start = integer(), end = integer(), stringsAsFactors = FALSE)
	}
	if (is.null(protected_ranges)) {
		protected_ranges <- data.frame(start = integer(), end = integer(), stringsAsFactors = FALSE)
	}

	decision_rows <- list()
	applied_mutations <- list()
	working_sequence <- sequence
	resolved_motif_keys <- character()
	motif_list <- unique(toupper(as.character(motif_df$motif)))
	row_index <- 1

	for (hit_index in seq_len(nrow(motif_hits))) {
		hit <- motif_hits[hit_index, , drop = FALSE]
		if (sytogen_is_type_iv_motif(hit)) {
			decision_rows[[row_index]] <- data.frame(
				motif = hit$motif,
				motif_start = hit$start,
				motif_end = hit$end,
				motif_strand = hit$strand,
				edit_position = NA_integer_,
				before = "",
				after = "",
				original_codon = "",
				replacement_codon = "",
				AA_LetterCode = "",
				synonymous = "",
				motifs_destroyed = 0,
				reasoning = "Type IV motif left unchanged.",
				motifs_created = 0,
				usage_score = 0,
				gc_preserving = FALSE,
				total_score = 0,
				chosen = FALSE,
				skip_reason = "type_iv",
				attempted_count = 0,
				rejected_count = 0,
				top_rejection_reason = "",
				top_rejection_count = 0,
				stringsAsFactors = FALSE
			)
			row_index <- row_index + 1
			next
		}

		candidate_rows <- list()
		candidate_scores <- numeric()
		candidate_counter <- 1
		for (position in seq(hit$start, hit$end)) {
			if (sytogen_position_in_ranges(position, protected_ranges)) {
				next
			}
			if (sytogen_position_in_ranges(position, cds_ranges)) {
				codon_window <- locate_codon_window(position, cds_ranges)
				if (is.null(codon_window)) {
					next
				}
				original_codon <- substr(working_sequence, codon_window$start, codon_window$end)
						if (!nzchar(original_codon) || nchar(original_codon) != 3 || !original_codon %in% names(sytogen_standard_genetic_code())) {
					next
				}
							aa <- sytogen_standard_genetic_code()[[original_codon]]
							synonymous_codons <- names(sytogen_standard_genetic_code())[sytogen_standard_genetic_code() == aa & names(sytogen_standard_genetic_code()) != original_codon]
				for (replacement_codon in synonymous_codons) {
					diff_positions <- which(strsplit(original_codon, "", fixed = TRUE)[[1]] != strsplit(replacement_codon, "", fixed = TRUE)[[1]])
					if (length(diff_positions) != 1) {
						next
					}
					mutation_position <- codon_window$start + diff_positions[1] - 1
					old_base <- substr(working_sequence, mutation_position, mutation_position)
					new_base <- substr(replacement_codon, diff_positions[1], diff_positions[1])
					mutated_sequence <- mutate_base_at(working_sequence, mutation_position, new_base)
					if (!sytogen_exact_motif_destroyed(working_sequence, mutated_sequence, hit)) {
						next
					}
					created_patterns <- sytogen_candidate_hits_created(working_sequence, mutated_sequence, motif_list, topology = topology)
					if (length(created_patterns) > 0) {
						next
					}
					usage_score <- if (length(codon_usage) > 0 && replacement_codon %in% names(codon_usage)) as.numeric(codon_usage[[replacement_codon]]) else 0
					gc_preserving <- base_class(old_base) == base_class(new_base)
					total_score <- 1000 + usage_score * 100 - 10 + if (preserve_gc && gc_preserving) 5 else 0
					candidate_rows[[candidate_counter]] <- data.frame(
						motif = hit$motif,
						motif_start = hit$start,
						motif_end = hit$end,
						motif_strand = hit$strand,
						edit_position = mutation_position,
						before = old_base,
						after = new_base,
						original_codon = original_codon,
						replacement_codon = replacement_codon,
						AA_LetterCode = aa,
						synonymous = TRUE,
						motifs_destroyed = 1,
						reasoning = "Synonymous coding edit",
						motifs_created = 0,
						usage_score = usage_score,
						gc_preserving = gc_preserving,
						total_score = total_score,
						chosen = FALSE,
						skip_reason = "",
						attempted_count = 1,
						rejected_count = 0,
						top_rejection_reason = "",
						top_rejection_count = 0,
						stringsAsFactors = FALSE
					)
					candidate_scores[candidate_counter] <- total_score
					candidate_counter <- candidate_counter + 1
				}
			} else {
				old_base <- substr(working_sequence, position, position)
				for (new_base in setdiff(c("A", "C", "G", "T"), old_base)) {
					mutated_sequence <- mutate_base_at(working_sequence, position, new_base)
					if (!sytogen_exact_motif_destroyed(working_sequence, mutated_sequence, hit)) {
						next
					}
					created_patterns <- sytogen_candidate_hits_created(working_sequence, mutated_sequence, motif_list, topology = topology)
					if (length(created_patterns) > 0) {
						next
					}
					gc_preserving <- base_class(old_base) == base_class(new_base)
					total_score <- 100 + if (preserve_gc && gc_preserving) 5 else 0
					candidate_rows[[candidate_counter]] <- data.frame(
						motif = hit$motif,
						motif_start = hit$start,
						motif_end = hit$end,
						motif_strand = hit$strand,
						edit_position = position,
						before = old_base,
						after = new_base,
						original_codon = "",
						replacement_codon = "",
						AA_LetterCode = "",
						synonymous = FALSE,
						motifs_destroyed = 1,
						reasoning = "Non-coding edit",
						motifs_created = 0,
						usage_score = 0,
						gc_preserving = gc_preserving,
						total_score = total_score,
						chosen = FALSE,
						skip_reason = "",
						attempted_count = 1,
						rejected_count = 0,
						top_rejection_reason = "",
						top_rejection_count = 0,
						stringsAsFactors = FALSE
					)
					candidate_scores[candidate_counter] <- total_score
					candidate_counter <- candidate_counter + 1
				}
			}
		}

		if (length(candidate_rows) == 0) {
			decision_rows[[row_index]] <- data.frame(
				motif = hit$motif,
				motif_start = hit$start,
				motif_end = hit$end,
				motif_strand = hit$strand,
				edit_position = NA_integer_,
				before = "",
				after = "",
				original_codon = "",
				replacement_codon = "",
				AA_LetterCode = "",
				synonymous = "",
				motifs_destroyed = 0,
				reasoning = "No valid candidate could be constructed for this motif.",
				motifs_created = 0,
				usage_score = 0,
				gc_preserving = FALSE,
				total_score = 0,
				chosen = FALSE,
				skip_reason = "no_valid_candidate",
				attempted_count = 0,
				rejected_count = 0,
				top_rejection_reason = "",
				top_rejection_count = 0,
				stringsAsFactors = FALSE
			)
			row_index <- row_index + 1
			next
		}

		candidate_df <- do.call(rbind, candidate_rows)
		best_idx <- order(-candidate_df$total_score, candidate_df$edit_position)[1]
		candidate_df$chosen <- FALSE
		candidate_df$chosen[best_idx] <- TRUE
		decision_rows <- c(decision_rows, split(candidate_df, seq_len(nrow(candidate_df))))
		row_index <- length(decision_rows) + 1
		best_row <- candidate_df[best_idx, , drop = FALSE]
		working_sequence <- mutate_base_at(working_sequence, as.integer(best_row$edit_position), as.character(best_row$after))
		applied_mutations[[length(applied_mutations) + 1]] <- list(position = as.integer(best_row$edit_position), old = as.character(best_row$before), new = as.character(best_row$after))
		resolved_motif_keys <- c(resolved_motif_keys, paste(hit$motif, hit$start, hit$end, hit$strand, sep = "|"))
	}

	decision_matrix <- if (length(decision_rows) > 0) do.call(rbind, decision_rows) else data.frame()
	new_motifs <- sytogen_candidate_hits_created(sequence, working_sequence, motif_list, topology = topology)
	list(
		decision_matrix = decision_matrix,
		applied_mutations = applied_mutations,
		mutated_sequence = working_sequence,
		new_motifs = new_motifs,
		resolved_motif_keys = unique(resolved_motif_keys)
	)
}

locate_codon_window <- function(position, cds_ranges = NULL) {
	if (is.null(cds_ranges) || nrow(cds_ranges) == 0) {
		return(NULL)
	}

	for (i in seq_len(nrow(cds_ranges))) {
		start <- as.integer(cds_ranges$start[i])
		end <- as.integer(cds_ranges$end[i])
		if (is.na(start) || is.na(end) || position < start || position > end) {
			next
		}

		codon_start <- start + ((position - start) %/% 3) * 3
		codon_end <- codon_start + 2
		if (codon_end <= end) {
			return(list(start = codon_start, end = codon_end, offset = position - codon_start + 1))
		}
	}

	NULL
}

score_candidate_mutation <- function(sequence, position, new_base, cds_ranges = NULL, preferred_codons = NULL) {
	old_base <- substr(sequence, position, position)
	if (old_base == new_base) {
		return(list(score = -Inf, synonymous = FALSE, codon = NA_character_, aa = NA_character_))
	}

	score <- 0
	synonymous <- FALSE
	mutated_codon <- NA_character_
	amino_acid <- NA_character_

	codon_window <- locate_codon_window(position, cds_ranges)
	if (!is.null(codon_window)) {
		original_codon <- substr(sequence, codon_window$start, codon_window$end)
		candidate_codon <- original_codon
		substr(candidate_codon, codon_window$offset, codon_window$offset) <- new_base

		if (original_codon %in% names(sytogen_standard_genetic_code()) && candidate_codon %in% names(sytogen_standard_genetic_code())) {
			aa_before <- sytogen_standard_genetic_code()[[original_codon]]
			aa_after <- sytogen_standard_genetic_code()[[candidate_codon]]
			amino_acid <- aa_before
			mutated_codon <- candidate_codon
			if (aa_before == aa_after) {
				synonymous <- TRUE
				score <- score + 3
				if (!is.null(preferred_codons) && !is.na(preferred_codons[[aa_before]]) && candidate_codon == preferred_codons[[aa_before]]) {
					score <- score + 2
				}
			} else {
				score <- score - 4
			}
		}
	} else {
		# Outside coding regions, prefer GC/AT-preserving substitutions.
		if (base_class(old_base) == base_class(new_base)) {
			score <- score + 2
		} else {
			score <- score + 1
		}
	}

	list(score = score, synonymous = synonymous, codon = mutated_codon, aa = amino_acid)
}

pick_best_mutation_for_hit <- function(sequence, hit, cds_ranges = NULL, preferred_codons = NULL) {
	motif_start <- as.integer(hit$start)
	motif_end <- as.integer(hit$end)
	candidate_bases <- c("A", "C", "G", "T")
	motif <- as.character(hit$motif)

	candidates <- lapply(seq(motif_start, motif_end), function(position) {
		current_base <- substr(sequence, position, position)
		lapply(setdiff(candidate_bases, current_base), function(new_base) {
			mutated_sequence <- mutate_base_at(sequence, position, new_base)
			mutated_fragment <- substr(mutated_sequence, motif_start, motif_end)
			if (mutated_fragment == motif) {
				return(NULL)
			}

			scored <- score_candidate_mutation(
				sequence = sequence,
				position = position,
				new_base = new_base,
				cds_ranges = cds_ranges,
				preferred_codons = preferred_codons
			)

			data.frame(
				motif = motif,
				hit_start = motif_start,
				hit_end = motif_end,
				strand = as.character(hit$strand),
				mutation_position = position,
				from_base = current_base,
				to_base = new_base,
				score = scored$score,
				synonymous = scored$synonymous,
				codon = ifelse(is.na(scored$codon), "", scored$codon),
				amino_acid = ifelse(is.na(scored$aa), "", scored$aa),
				stringsAsFactors = FALSE
			)
		})
	})

	candidates <- unlist(candidates, recursive = FALSE)
	candidates <- Filter(Negate(is.null), candidates)
	if (length(candidates) == 0) {
		return(NULL)
	}

	candidate_df <- do.call(rbind, candidates)
	best_idx <- order(-candidate_df$score, candidate_df$mutation_position)[1]
	candidate_df[best_idx, , drop = FALSE]
}

apply_non_overlapping_mutations <- function(sequence, decision_matrix) {
	if (is.null(decision_matrix) || nrow(decision_matrix) == 0) {
		return(sequence)
	}

	ordered <- decision_matrix[order(-decision_matrix$score, decision_matrix$mutation_position), , drop = FALSE]
	used_positions <- integer(0)
	mutated <- sequence

	for (i in seq_len(nrow(ordered))) {
		pos <- as.integer(ordered$mutation_position[i])
		if (pos %in% used_positions) {
			next
		}
		mutated <- mutate_base_at(mutated, pos, as.character(ordered$to_base[i]))
		used_positions <- c(used_positions, pos)
	}

	mutated
}

 sytogen_build_motif_summary <- function(motif_hits, resolved_motif_keys) {
	if (is.null(motif_hits) || nrow(motif_hits) == 0) {
		return(data.frame(motif = character(), type = character(), total_hits = integer(), resolved = integer(), unresolved = integer(), stringsAsFactors = FALSE))
	}

	rows <- list()
	group_keys <- unique(paste(motif_hits$motif, ifelse(nzchar(as.character(motif_hits$enz_type)), as.character(motif_hits$enz_type), "unknown"), sep = "|"))
	for (group_key in group_keys) {
		parts <- strsplit(group_key, "|", fixed = TRUE)[[1]]
		motif_name <- parts[1]
		type_name <- parts[2]
		group_hits <- motif_hits[toupper(motif_hits$motif) == toupper(motif_name) & ifelse(nzchar(as.character(motif_hits$enz_type)), as.character(motif_hits$enz_type), "unknown") == type_name, , drop = FALSE]
		resolved_count <- sum(paste(group_hits$motif, group_hits$start, group_hits$end, group_hits$strand, sep = "|") %in% resolved_motif_keys)
		rows[[length(rows) + 1]] <- data.frame(
			motif = motif_name,
			type = type_name,
			total_hits = nrow(group_hits),
			resolved = resolved_count,
			unresolved = nrow(group_hits) - resolved_count,
			stringsAsFactors = FALSE
		)
	}
	if (length(rows) == 0) {
		return(data.frame(motif = character(), type = character(), total_hits = integer(), resolved = integer(), unresolved = integer(), stringsAsFactors = FALSE))
	}
	do.call(rbind, rows)
}

run_sytogen_pipeline <- function(sequence,
										 codon_df = NULL,
										 motif_df = NULL,
										 params = list()) {
	cleaned_sequence <- sytogen_clean_sequence(sequence)
	topology <- tolower(as.character(params$topology %||% "circular"))
	if (!topology %in% c("linear", "circular")) {
		stop("topology must be 'linear' or 'circular'.")
	}
	preserve_gc <- isTRUE(params$preserve_gc)
	sequence_length <- nchar(cleaned_sequence)

	if (is.null(motif_df)) {
		stop("motif_df is required.")
	}

	motif_hits <- sytogen_parse_motif_table(motif_df, cleaned_sequence, topology = topology)
	if (!is.null(params$include_reverse) && !isTRUE(params$include_reverse)) {
		motif_hits <- motif_hits[motif_hits$strand != "-", , drop = FALSE]
	}
	codon_usage <- sytogen_parse_codon_usage(codon_df, cleaned_sequence)
	cds_ranges <- sytogen_normalize_ranges(params$cds_ranges, sequence_length)
	mask_ranges <- sytogen_normalize_ranges(params$mask_ranges, sequence_length)
	protected_ranges <- sytogen_normalize_ranges(params$protected_ranges, sequence_length)
	protected_override_ranges <- sytogen_normalize_ranges(params$protected_override_ranges, sequence_length)

	if (nrow(protected_override_ranges) > 0 && nrow(protected_ranges) > 0) {
		keep <- logical(nrow(protected_ranges))
		for (i in seq_len(nrow(protected_ranges))) {
			keep[i] <- !any(
				!(protected_ranges$end[i] < protected_override_ranges$start | protected_ranges$start[i] > protected_override_ranges$end)
			)
		}
		protected_ranges <- protected_ranges[keep, , drop = FALSE]
	}

	all_protected_ranges <- if (nrow(protected_ranges) == 0) mask_ranges else rbind(protected_ranges, mask_ranges)
	if (nrow(all_protected_ranges) > 0) {
		motif_hits <- motif_hits[!vapply(seq_len(nrow(motif_hits)), function(i) {
			sytogen_position_in_ranges(motif_hits$start[i], all_protected_ranges) || sytogen_position_in_ranges(motif_hits$end[i], all_protected_ranges)
		}, logical(1)), , drop = FALSE]
	}

	parsed_motif_df <- motif_df
	if (is.null(parsed_motif_df$enz_type)) {
		parsed_motif_df$enz_type <- ""
	}
	parsed_motif_df$motif <- toupper(trimws(as.character(parsed_motif_df$motif)))

	decisions <- sytogen_collect_decisions(
		sequence = cleaned_sequence,
		motif_hits = motif_hits,
		motif_df = parsed_motif_df,
		cds_ranges = cds_ranges,
		protected_ranges = all_protected_ranges,
		codon_usage = codon_usage,
		preserve_gc = preserve_gc,
		topology = topology
	)

	mutated_sequence <- decisions$mutated_sequence
	new_motifs <- decisions$new_motifs
	motif_summary <- sytogen_build_motif_summary(motif_hits, decisions$resolved_motif_keys)
	summary <- list(
		sequence_id = "sequence",
		topology = topology,
		original_length = sequence_length,
		altered_length = nchar(mutated_sequence),
		motifs_input = nrow(motif_hits),
		motifs_resolved = length(unique(decisions$resolved_motif_keys)),
		motifs_unresolved = max(0, nrow(motif_hits) - length(unique(decisions$resolved_motif_keys))),
		edits_applied = length(decisions$applied_mutations),
		candidates_total = if (is.null(decisions$decision_matrix)) 0 else nrow(decisions$decision_matrix),
		new_motifs_introduced = length(new_motifs),
		mask_regions_applied = nrow(mask_ranges),
		protected_override_ranges_applied = nrow(protected_override_ranges)
	)

	list(
		sequence = cleaned_sequence,
		motifs = parsed_motif_df$motif,
		motif_hits = motif_hits,
		codon_bias = codon_usage,
		decision_matrix = decisions$decision_matrix,
		mutated_sequence = mutated_sequence,
		applied_mutations = decisions$applied_mutations,
		new_motifs = new_motifs,
		summary = summary,
		motif_summary = motif_summary,
		assembly_plan = NULL,
		mask_regions = mask_ranges,
		protected_override_ranges = protected_override_ranges
	)
}

run_sytogen <- function(sequence,
						 motif_text = NULL,
						 motifs = NULL,
						 cds_ranges = NULL,
						 include_reverse = TRUE) {
	motif_df <- NULL
	if (!is.null(motif_text)) {
				parsed <- try(sytogen_parse_motif_text(motif_text), silent = TRUE)
		if (!inherits(parsed, "try-error") && nrow(parsed) > 0) {
			motif_df <- parsed
		}
	}
	if (is.null(motif_df)) {
		motif_df <- data.frame(motif = sytogen_extract_motifs(motif_text = motif_text, motifs = motifs), stringsAsFactors = FALSE)
	}

	params <- list(
		topology = "circular",
		preserve_gc = FALSE,
		include_reverse = include_reverse,
		cds_ranges = cds_ranges
	)
	run_sytogen_pipeline(sequence = sequence, codon_df = NULL, motif_df = motif_df, params = params)
}

sytogen_runner <- function(sequence,
						 motif_file = NULL,
						 motifs = NULL,
						 cds_ranges = NULL,
						 include_reverse = TRUE) {
	motif_text <- NULL
	if (!is.null(motif_file)) {
		if (!file.exists(motif_file)) {
			stop(sprintf("File not found: %s", motif_file))
		}
		motif_text <- paste(readLines(motif_file, warn = FALSE), collapse = "\n")
	}

	run_sytogen(
		sequence = sequence,
		motif_text = motif_text,
		motifs = motifs,
		cds_ranges = cds_ranges,
		include_reverse = include_reverse
	)
}