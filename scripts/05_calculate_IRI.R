#!/usr/bin/env Rscript
# ============================================================
# Step 5: Calculate IRI for each intron in each sample
# Input : results/quant/{LD1,LD2,D1,D2}/abundance.tsv + justAS_intron_map.tsv
# Output: results/JustAS_IRI_all_samples.tsv
# ============================================================

OUTPUT_DIR <- "../results"
samples <- c("LD1", "LD2", "D1", "D2")

map_file <- file.path(OUTPUT_DIR, "justAS_intron_map.tsv")
map <- read.delim(map_file, stringsAsFactors = FALSE, check.names = FALSE)

iri_list <- list()
for (s in samples) {
  ab <- read.delim(file.path(OUTPUT_DIR, "quant", s, "abundance.tsv"),
                   stringsAsFactors = FALSE)
  ab <- ab[, c("target_id", "est_counts")]
  colnames(ab)[2] <- s
  
  m <- merge(map, ab, by.x = "fragment_id", by.y = "target_id", all.x = TRUE)
  m[[s]][is.na(m[[s]])] <- 0
  
  agg <- aggregate(m[[s]], by = list(intron_id = m$intron_id, junction_type = m$junction_type),
                   FUN = sum)
  colnames(agg)[3] <- "count"
  wide <- reshape(agg, idvar = "intron_id", timevar = "junction_type", direction = "wide")
  colnames(wide) <- sub("^count\\.", "", colnames(wide))
  for (jt in c("Ex_In", "In_Ex", "Ex_Ex")) {
    if (!jt %in% colnames(wide)) wide[[jt]] <- 0
    wide[[jt]][is.na(wide[[jt]])] <- 0
  }
  
  ex_in <- wide$Ex_In; in_ex <- wide$In_Ex; ex_ex <- wide$Ex_Ex
  denom <- ex_ex + 0.5 * ex_in + 0.5 * in_ex
  iri   <- ifelse(denom > 0, 0.5 * (ex_in + in_ex) / denom, 0)
  
  res <- data.frame(intron_id = wide$intron_id,
                    Ex_In = ex_in, In_Ex = in_ex, Ex_Ex = ex_ex, IRI = iri,
                    stringsAsFactors = FALSE)
  colnames(res)[2:5] <- paste0(s, "_", colnames(res)[2:5])
  iri_list[[s]] <- res
  cat(s, ": introns =", nrow(res), " IRI>0 =", sum(res[[paste0(s, "_IRI")]] > 0), "\n")
}

iri_all <- iri_list[[1]]
for (s in samples[-1]) iri_all <- merge(iri_all, iri_list[[s]], by = "intron_id", all = TRUE)
iri_all <- iri_all[order(iri_all$intron_id), ]
out <- file.path(OUTPUT_DIR, "JustAS_IRI_all_samples.tsv")
write.table(iri_all, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Output:", out, "\n")