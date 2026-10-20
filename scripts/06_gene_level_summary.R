#!/usr/bin/env Rscript
# ============================================================
# Step 6: Map intron-level IRI to gene level
# Input : GTF + results/JustAS_IRI_all_samples.tsv
# Output: results/JustAS_gene_level_summary.tsv
# ============================================================

INPUT_DIR  <- "../input"
OUTPUT_DIR <- "../results"

gtf_file <- file.path(INPUT_DIR, "IRGSP-1.0_representative_transcript_exon_2022-03-11.gtf")
iri_file <- file.path(OUTPUT_DIR, "JustAS_IRI_all_samples.tsv")
out_file <- file.path(OUTPUT_DIR, "JustAS_gene_level_summary.tsv")

lines <- readLines(gtf_file); lines <- lines[!grepl("^#", lines)]
f <- strsplit(lines, "\t", fixed = TRUE)
f <- f[sapply(f, length) >= 9 & sapply(f, `[`, 3) == "exon"]
attr <- sapply(f, `[`, 9)
tx2gene <- unique(data.frame(
  transcript_id = sub('.*transcript_id "([^"]+)".*', '\\1', attr),
  gene_id       = sub('.*gene_id "([^"]+)".*', '\\1', attr),
  stringsAsFactors = FALSE
))

iri <- read.delim(iri_file, stringsAsFactors = FALSE, check.names = FALSE)
iri$transcript_id <- sub(" .*$", "", iri$intron_id)
iri$gene_id <- tx2gene$gene_id[match(iri$transcript_id, tx2gene$transcript_id)]

samples <- c("LD1", "LD2", "D1", "D2")
gene_summary <- data.frame(gene_id = unique(iri$gene_id[!is.na(iri$gene_id)]),
                           stringsAsFactors = FALSE)
for (s in samples) {
  col <- paste0(s, "_IRI")
  max_iri <- tapply(iri[[col]], iri$gene_id, max, na.rm = TRUE)
  n_pos   <- tapply(iri[[col]] > 0, iri$gene_id, sum, na.rm = TRUE)
  gene_summary[[paste0(s, "_maxIRI")]]          <- max_iri[gene_summary$gene_id]
  gene_summary[[paste0(s, "_n_intron_IRI_gt0")]] <- n_pos[gene_summary$gene_id]
  cat(s, ": genes with IRI>0 =", length(unique(iri$gene_id[!is.na(iri$gene_id) & iri[[col]] > 0])), "\n")
}

write.table(gene_summary, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Output:", out_file, "\n")