#!/usr/bin/env Rscript
# ============================================================
# Step 2: Extract 40-nt junction fragments (Ex_In, In_Ex, Ex_Ex)
# Input : intron_info.tsv + 12 chromosome FASTA files
# Output: intron_junctions.tsv
# ============================================================

suppressPackageStartupMessages(library(Biostrings))

INPUT_DIR  <- "../input"
OUTPUT_DIR <- "../results"

intron_file <- file.path(OUTPUT_DIR, "intron_info.tsv")
out_file    <- file.path(OUTPUT_DIR, "intron_junctions.tsv")

# 12 chromosome FASTA files (chr01.fasta ... chr12.fasta)
chr_files <- file.path(INPUT_DIR, sprintf("chr%02d.fasta", 1:12))
names(chr_files) <- sprintf("chr%02d", 1:12)
missing <- chr_files[!file.exists(chr_files)]
if (length(missing) > 0) stop("Missing FASTA files:\n", paste(missing, collapse = "\n"))

introns <- read.delim(intron_file, stringsAsFactors = FALSE, check.names = FALSE)
introns$chr_num <- as.integer(sub(".*\\(([0-9]+):.*", "\\1", introns[["True intron ID"]]))

istart  <- as.integer(introns[["Start of intron_n"]])
iend    <- as.integer(introns[["End of intron_n"]])
a_start <- as.integer(introns[["Start of 20bp upstream the end of Exon_n"]])
b_end   <- as.integer(introns[["20bp downstream the start of intron_n"]])
c_start <- as.integer(introns[["Start of 20bp upstream the end of Intron_n"]])
d_end   <- as.integer(introns[["20bp downstream the start of Exon_n+1"]])

exon_intron <- rep(NA_character_, nrow(introns))
intron_exon <- rep(NA_character_, nrow(introns))
exon_exon   <- rep(NA_character_, nrow(introns))

for (n in sort(unique(introns$chr_num))) {
  key <- sprintf("chr%02d", n)
  seq_chr <- readDNAStringSet(chr_files[[key]])[[1]]
  seq_len <- length(seq_chr)
  idx <- which(introns$chr_num == n)
  cat("  ", key, "(", length(idx), "introns)... ")
  n_bad <- 0L
  for (j in idx) {
    a_s <- a_start[j]; a_e <- istart[j] - 1L
    b_s <- istart[j]; b_e <- b_end[j]
    c_s <- c_start[j]; c_e <- iend[j]
    d_s <- iend[j] + 1L; d_e <- d_end[j]
    if (a_s < 1 || a_e > seq_len || b_s < 1 || b_e > seq_len ||
        c_s < 1 || c_e > seq_len || d_s < 1 || d_e > seq_len) {
      n_bad <- n_bad + 1L; next
    }
    a <- as.character(subseq(seq_chr, a_s, a_e))
    b <- as.character(subseq(seq_chr, b_s, b_e))
    c <- as.character(subseq(seq_chr, c_s, c_e))
    d <- as.character(subseq(seq_chr, d_s, d_e))
    exon_intron[j] <- paste0(a, b)
    intron_exon[j] <- paste0(c, d)
    exon_exon[j]   <- paste0(a, d)
  }
  cat("done (", n_bad, " out-of-range)\n", sep = "")
  rm(seq_chr); gc()
}

result <- data.frame(
  "True intron ID"   = introns[["True intron ID"]],
  "Exon_intron_40bp" = exon_intron,
  "Intron_exon_40bp" = intron_exon,
  "Exon_exon_40bp"   = exon_exon,
  stringsAsFactors = FALSE, check.names = FALSE
)
len_ok <- nchar(result$Exon_intron_40bp) == 40 &
  nchar(result$Intron_exon_40bp) == 40 &
  nchar(result$Exon_exon_40bp) == 40
len_ok[is.na(len_ok)] <- FALSE
cat("\nTotal:", nrow(result), " valid 40-nt:", sum(len_ok), "\n")
write.table(result, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("Output:", out_file, "\n")