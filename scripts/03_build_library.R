#!/usr/bin/env Rscript
# ============================================================
# Step 3: Build JustAS junction library
#   - Collapse introns by (gene_id, coordinates, strand)
#   - Negative strand: reverse-complement + swap Ex_In/In_Ex labels
#   - Intron numbering kept in genomic order
# Input : GTF + intron_junctions.tsv
# Output: justAS_junctions.fa + justAS_intron_map.tsv
# ============================================================

suppressPackageStartupMessages(library(Biostrings))

INPUT_DIR  <- "../input"
OUTPUT_DIR <- "../results"

gtf_file <- file.path(INPUT_DIR, "IRGSP-1.0_representative_transcript_exon_2022-03-11.gtf")
in_file  <- file.path(OUTPUT_DIR, "intron_junctions.tsv")
out_fa   <- file.path(OUTPUT_DIR, "justAS_junctions.fa")
out_map  <- file.path(OUTPUT_DIR, "justAS_intron_map.tsv")

# ---- transcript_id -> gene_id ----
lines <- readLines(gtf_file); lines <- lines[!grepl("^#", lines)]
f <- strsplit(lines, "\t", fixed = TRUE)
f <- f[sapply(f, length) >= 9 & sapply(f, `[`, 3) == "exon"]
attr <- sapply(f, `[`, 9)
tx2gene <- unique(data.frame(
  transcript_id = sub('.*transcript_id "([^"]+)".*', '\\1', attr),
  gene_id       = sub('.*gene_id "([^"]+)".*', '\\1', attr),
  stringsAsFactors = FALSE
))

# ---- parse intron_junctions.tsv ----
junc <- read.delim(in_file, stringsAsFactors = FALSE, check.names = FALSE)
id <- junc[["True intron ID"]]
junc$transcript_id <- sub(" .*$", "", id)
junc$coord         <- sub(".*\\(([0-9]+:[0-9]+-[0-9]+)\\).*", "\\1", id)
junc$strand        <- sub(".*\\)([+-])$", "\\1", id)
junc$coord_strand  <- paste0(junc$coord, junc$strand)
junc$gene_id       <- tx2gene$gene_id[match(junc$transcript_id, tx2gene$transcript_id)]
junc$gene_id[is.na(junc$gene_id)] <- junc$transcript_id[is.na(junc$gene_id)]

# ---- collapse by (gene_id, coord_strand) ----
junc <- junc[order(junc$gene_id, junc$coord_strand, junc$transcript_id), ]
keep <- !duplicated(paste(junc$gene_id, junc$coord_strand, sep = "|"))
junc_unique <- junc[keep, ]

# ---- expand to 3 fragments ----
long <- rbind(
  data.frame(intron_id = junc_unique[["True intron ID"]],
             sequence = junc_unique[["Exon_intron_40bp"]],
             junction_type = "Ex_In", stringsAsFactors = FALSE),
  data.frame(intron_id = junc_unique[["True intron ID"]],
             sequence = junc_unique[["Intron_exon_40bp"]],
             junction_type = "In_Ex", stringsAsFactors = FALSE),
  data.frame(intron_id = junc_unique[["True intron ID"]],
             sequence = junc_unique[["Exon_exon_40bp"]],
             junction_type = "Ex_Ex", stringsAsFactors = FALSE)
)
long$strand <- rep(junc_unique$strand, 3)
long <- long[!is.na(long$sequence) & long$sequence != "" & nchar(long$sequence) == 40, ]

# ---- negative strand: revcomp + swap labels ----
neg <- long$strand == "-"
if (any(neg)) {
  long$sequence[neg] <- as.character(reverseComplement(DNAStringSet(long$sequence[neg])))
  tmp_in  <- long$junction_type[neg] == "Ex_In"
  tmp_iex <- long$junction_type[neg] == "In_Ex"
  long$junction_type[neg][tmp_in]  <- "In_Ex"
  long$junction_type[neg][tmp_iex] <- "Ex_In"
}

# ---- write outputs ----
long$fragment_id <- sprintf("JX%07d", seq_len(nrow(long)))
con <- file(out_fa, "w")
for (i in seq_len(nrow(long))) {
  writeLines(paste0(">", long$fragment_id[i]), con)
  writeLines(long$sequence[i], con)
}
close(con)

write.table(long[, c("fragment_id", "sequence", "intron_id", "junction_type")],
            out_map, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Introns:", nrow(junc_unique), "\n")
cat("Fragments:", nrow(long), "\n")
cat("Output:", out_fa, "and", out_map, "\n")