# JustAS: Intron Retention Detection from RNA-seq

JustAS is a lightweight pipeline for detecting intron retention (IR) events from RNA-seq data. It builds a custom junction library from annotated introns, quantifies exon–intron, intron–exon, and exon–exon junction fragments using kallisto, and calculates Intron Retention Intensity (IRI) for each intron.

## Requirements

- R >= 4.2.1
- Biostrings (Bioconductor)
- kallisto v0.46.1
- Rice reference: IRGSP-1.0 genome + MSU 7.0 / IRGSP-1.0 representative transcript annotation

## Preparing input files

Place the following files in a directory named `input/` next to `scripts/` and `results/`:

- `IRGSP-1.0_representative_transcript_exon_2022-03-11.gtf` (annotation)
- `chr01.fasta` ... `chr12.fasta` (genome, one FASTA per chromosome)
- Paired-end FASTQ files for each sample (referenced in `scripts/04_kallisto_commands.txt`)

A small example (`example_data/chr01_first100introns.tsv`) is provided to illustrate the `intron_junctions.tsv` format; it can be used to test Steps 3, 5, and 7 without the full genome.

## Workflow

### Step 1: Extract introns from GTF
```bash
Rscript scripts/01_extract_introns.R