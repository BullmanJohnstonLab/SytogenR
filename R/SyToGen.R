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

	parsed_table <- try(parse_motif_text(motif_text), silent = TRUE)
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

		if (original_codon %in% names(standard_genetic_code) && candidate_codon %in% names(standard_genetic_code)) {
			aa_before <- standard_genetic_code[[original_codon]]
			aa_after <- standard_genetic_code[[candidate_codon]]
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

run_sytogen <- function(sequence,
												motif_text = NULL,
												motifs = NULL,
												cds_ranges = NULL,
												include_reverse = TRUE) {
	cleaned_sequence <- sytogen_clean_sequence(sequence)
	motif_list <- sytogen_extract_motifs(motif_text = motif_text, motifs = motifs)

	motif_hits <- find_motifs(cleaned_sequence, motif_list, include_reverse = include_reverse)
	codon_bias <- calculate_codon_bias(cleaned_sequence)

	if (is.null(cds_ranges)) {
		cds_ranges <- data.frame(start = integer(), end = integer(), stringsAsFactors = FALSE)
	}

	if (nrow(motif_hits) == 0) {
		return(list(
			sequence = cleaned_sequence,
			motifs = motif_list,
			motif_hits = motif_hits,
			codon_bias = codon_bias,
			decision_matrix = data.frame(),
			mutated_sequence = cleaned_sequence
		))
	}

	decisions <- lapply(seq_len(nrow(motif_hits)), function(i) {
		pick_best_mutation_for_hit(
			sequence = cleaned_sequence,
			hit = motif_hits[i, , drop = FALSE],
			cds_ranges = cds_ranges,
			preferred_codons = codon_bias$preferred_codons
		)
	})
	decisions <- Filter(Negate(is.null), decisions)

	decision_matrix <- if (length(decisions) > 0) do.call(rbind, decisions) else data.frame()
	mutated_sequence <- apply_non_overlapping_mutations(cleaned_sequence, decision_matrix)

	list(
		sequence = cleaned_sequence,
		motifs = motif_list,
		motif_hits = motif_hits,
		codon_bias = codon_bias,
		decision_matrix = decision_matrix,
		mutated_sequence = mutated_sequence
	)
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