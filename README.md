# snpdep
Rust Program to calculate haplotype-specific depth from BAM/CRAM files at positions provided in a haplotype-resolved VCF file.

## Usage 

```text
Program to calculate haplotype-specific depth from BAM/CRAM files at positions provided in a haplotype-resolved VCF file.

Usage: snpdep [OPTIONS] --output <OUTPUT_VCF> <INPUT_VCF> <READS>

Arguments:
  <INPUT_VCF>  Path to VCF file containing haplotype-resolved variants. Both uncompressed `.vcf` and compressed `.vcf.gz` files are supported
  <READS>      Path to BAM/CRAM reads file

Options:
  -o, --output <OUTPUT_VCF>
          Path to the output annotated VCF file. Compression will be inferred from the file extension: `.vcf` for uncompressed and `.vcf.gz` for compressed output
  -f, --reads-format <READS_FORMAT>
          Specifies the reads file format [default: bam] [possible values: bam, cram]
  -r, --reference <FASTA>
          Path to the reference FASTA file (required if `--reads-format cram`)
      --min-mapq <MIN_MAPQ>
          Minimum read mapping quality (MAPQ) filtration threshold [default: 1]
      --min-count <MIN_COUNT>
          Minimum read coverage required at a position to be annotated [default: 1]
  -c, --chunksize <CHUNK_SIZE>
          Number of SNPs processed by a thread in each iteration [default: 1000]
      --format-field-id <FORMAT_ID>
          Unique ID for the new VCF format field (e.g RNA, K4me3), max 5 characters [default: RC]
      --format-field-name <FORMAT_NAME>
          Short description for new VCF format field (e.g RNAseq, H3K4me3) [default: BAM/CRAM]
  -t, --threads <NUM_THREADS>
          Number of threads to use [default: 1]
  -h, --help
          Print help
  -V, --version
          Print version
```