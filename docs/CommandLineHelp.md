# Command-Line Help for `snpdep`

This document contains the help content for the `snpdep` command-line program.

**Command Overview:**

* [`snpdep`↴](#snpdep)

## `snpdep`

Program to calculate haplotype-specific depth from BAM/CRAM files at positions provided in a haplotype-resolved VCF file.

**Usage:** `snpdep [OPTIONS] --output <OUTPUT_VCP> <INPUT_VCF> <READS>`

###### **Arguments:**

* `<INPUT_VCF>` — Path to VCF file containing haplotype-resolved variants
* `<READS>` — Path to BAM/CRAM reads file

###### **Options:**

* `-o`, `--output <OUTPUT_VCP>` — Path to Output annotated VCF file
* `-f`, `--reads-format <READS_FORMAT>` — Reads file format

  Default value: `bam`

  Possible values: `bam`, `cram`

* `-r`, `--reference <FASTA>` — Path to the reference FASTA file (required if `--reads-format cram`)
* `--min-mapq <MIN_MAPQ>` — Minimum read mapping quality (MAPQ) filtration threshold

  Default value: `1`
* `--min-count <MIN_COUNT>` — Minimum read coverage required at a position to be annotated

  Default value: `1`
* `-c`, `--chunksize <CHUNK_SIZE>` — Number of SNPs processed by a thread each iteration

  Default value: `1000000`
* `--format-field-id <FORMAT_ID>` — Unique ID for new VCF format field (e.g RNA, K4me3), max 5 characters

  Default value: `RC`
* `--format-field-name <FORMAT_NAME>` — Short description for new VCF format field (e.g RNAseq, H3K4me3)

  Default value: `BAM/CRAM`
* `-t`, `--threads <NUM_THREADS>` — Number of threads to use

  Default value: `1`
* `--markdown-help`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>