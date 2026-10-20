#!/usr/bin/env Rscript
# ============================================================
# Step 1: Extract intron information from a GTF annotation
# Input : annotation GTF (exon features)
# Output: intron_info.tsv (one row per transcript-level intron)
# ============================================================

INPUT_DIR  <- "../input"          # change to your local path
OUTPUT_DIR <- "../results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

gtf_file <- file.path(INPUT_DIR, "IRGSP-1.0_representative_transcript_exon_2022-03-11.gtf")
out_file <- file.path(OUTPUT_DIR, "intron_info.tsv")

# ---- Read exon lines ----
lines  <- readLines(gtf_file)
lines  <- lines[!grepl("^#", lines)]
fields <- strsplit(lines, "\t", fixed = TRUE)
fields <- fields[sapply(fields, function(x) length(x) >= 9 && x[3] == "exon")]

chrom  <- sapply(fields, `[`, 1)
start  <- as.integer(sapply(fields, `[`, 4))
end    <- as.integer(sapply(fields, `[`, 5))
strand <- sapply(fields, `[`, 7)
attr   <- sapply(fields, `[`, 9)
tx_id  <- sub('.*transcript_id "([^"]+)".*', '\\1', attr)

exons <- data.frame(chrom, start, end, strand, transcript_id = tx_id,
                    stringsAsFactors = FALSE)

# ---- Build introns per transcript ----
tx_list <- split(exons, exons$transcript_id)
result_list <- list()
k <- 1L

for (tx in names(tx_list)) {
  ex <- tx_list[[tx]]
  ex <- ex[order(ex$start), ]
  if (nrow(ex) < 2) next
  
  chrom_tx  <- ex$chrom[1]
  strand_tx <- ex$strand[1]
  chrom_disp <- if (grepl("^chr[0-9]+$", chrom_tx)) sub("^chr0*", "", chrom_tx) else chrom_tx
  
  for (i in seq_len(nrow(ex) - 1)) {
    exon_n  <- ex[i, ]
    exon_n1 <- ex[i + 1, ]
    intron_start <- exon_n$end + 1
    intron_end   <- exon_n1$start - 1
    if (intron_end < intron_start) next
    intron_len <- intron_end - intron_start + 1
    
    result_list[[k]] <- data.frame(
      True_intron_ID                      = paste0(tx, " intron", i, "  (", chrom_disp,
                                                   ":", intron_start, "-", intron_end, ")", strand_tx),
      Start_of_intron_n                   = intron_start,
      End_of_intron_n                     = intron_end,
      Length_of_intron_n                  = intron_len,
      Start_of_20bp_upstream_end_Exon_n   = exon_n$end - 19,
      X20bp_downstream_start_intron_n     = intron_start + 19,
      Start_of_20bp_upstream_end_Intron_n = intron_end - 19,
      X20bp_downstream_start_Exon_n1      = exon_n1$start + 19,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

result <- do.call(rbind, result_list)
result$tx_tmp    <- sub(" intron.*", "", result$True_intron_ID)
result$intron_num <- as.integer(sub(".* intron([0-9]+).*", "\\1", result$True_intron_ID))
result <- result[order(result$tx_tmp, result$intron_num), ]
result$tx_tmp <- NULL; result$intron_num <- NULL

colnames(result) <- c(
  "True intron ID",
  "Start of intron_n", "End of intron_n", "Length of intron_n",
  "Start of 20bp upstream the end of Exon_n",
  "20bp downstream the start of intron_n",
  "Start of 20bp upstream the end of Intron_n",
  "20bp downstream the start of Exon_n+1"
)

write.table(result, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Done. Introns extracted:", nrow(result), "\nOutput:", out_file, "\n")