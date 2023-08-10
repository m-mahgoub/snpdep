# snpdep

snpdep is a high-speed bioinformatics tool to calculate haplotype-specific read depth from BAM/CRAM files using positions provided in a VCF file. Developed using the Rust programming language, snpdep is a standalone binary tool with no external dependencies and with multithreading support. This significantly accelerates runtime, allowing it to scale in tandem with the availability of computing resources.

## Installation

Pre-compiled binaries are provided for Linux and MacOS from the [release page](https://github.com/m-mahgoub/snpdep/releases)

### Building from source

To build `snpdep` from source [cargo](https://www.rust-lang.org/learn/get-started) should be used.

```bash
git clone https://github.com/m-mahgoub/snpdep.git
cd snpdep
cargo install --path .
# or
cargo install --git https://github.com/m-mahgoub/snpdep.git
```

## Usage

To run, snpdep requires (1) a genotype VCF file and (2) a BAM/CRAM reads file. It then calculates the number of reads supporting each allele. Typically, the genotype VCF file is produced by calling germline variants from Whole Genome Sequencing reads on the sample to be examined with the BAM/CRAM read count. The BAM/CRAM can be of any "functional" data type, such as RNAseq, ChIP-seq, CUT&TAG, ATACseq, etc., for which haplotype-specific reads are required. This proves useful when investigating haplotype bias in these datasets.

While not mandatory, the usage of a haplotype-resolved VCF file as input is highly recommended. This is because it provides additional insight for studying haplotype bias across multiple functional datasets, such as DNA methylation vs. gene expression or chromatin accessibility.

snpdep operates with a single-sample VCF file (either compressed or uncompressed) and one BAM/CRAM file, both of which should have the same reference genome. It runs through all records in the VCF file and counts the reads supporting the reference and alternative allele. If there are multiple alternative alleles, only the first is considered. Only single nucleotide substitution variants are taken into account.

The output VCF file will include a new format field tag for each record (user-defined tag), along with haplotype counts in the sample's field. Also, the VCF file's header will contain a new line to define and describe the new format field. The default behavior is to annotate records with at least one supporting read (irrespective of whether it's from a reference or alternative allele). However, users can modify this by choosing to annotate all records or restrict annotation to variants with a specific count number using the `--min-count <MIN_COUNT>` option.

If the reads input file is in the CRAM format, the user must provide a path for the reference fasta `--reference <FASTA>` and specify the CRAM format option `--reads-format CRAM`. At present, CRAM format support is only available for the Linux OS.

```text
Usage: snpdep [OPTIONS] --output <OUTPUT_VCF> <INPUT_VCF> <READS>

Arguments:
  <INPUT_VCF>  Path to VCF file containing sample genotype. Both uncompressed `.vcf` and compressed `.vcf.gz` files are supported
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
          Minimum read coverage required at a position to be annotated. if set to 0, all positions will be annotated [default: 1]
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

## Example: Haplotype-Specific RNA-Seq Read Counting
This is a mini-test example for counting haplotype-specific RNA-seq reads at each biallelic position within a sample. Follow the steps below to download the required software and test data, then execute the analysis using different options.
#### Step 1: Download and Install snpdep
First, download the `snpdep` executable and add it to your system `$PATH`.
```bash
mkdir test
cd test
wget https://github.com/m-mahgoub/snpdep/releases/download/v0.1.0-alpha/snpdep_v0.1.0-alpha_linux.tar.gz
tar -xvf snpdep_v0.1.0-alpha_linux.tar.gz
chmod +x linux/snpdep
export PATH=$PATH:$PWD/linux
```
#### Step 2: Download Test Data
You will need the AWS Command Line Interface (CLI) installed to download the test data from AWS S3
```bash
aws s3 cp s3://davidspencerlab/test_data/snpdep/ . --recursive
```
#### Step 3: Execute Analysis
You can run snpdep with different options as demonstrated below.

3.1 Run with Default Parameters
```bash
snpdep --output output.vcf test.vcf test.bam
```
3.2 Run with CRAM Reads as Input
```bash
snpdep --reads-format cram --reference hg38.fa --output output.vcf test.vcf test.cram
```
3.3 Run with Multiple Options
```bash
snpdep  \
--threads 4                            `# Use 4 threads` \
--chunksize 10                         `# Chunk size of 10 (thread will process 10 VCF records per iteration)` \
--format-field-id K4me3                `# Format-ID that will be added in the record's format field` \
--format-field-name "H3K4me3 ChIPseq"  `# Format Description that will be added in the header` \
--reads-format cram                    `# Input reads format` \
--reference hg38.fa                    `# Reference FASTA file` \
--output output.vcf.gz                 `# Output file (using ".gz" will automatically compress it)` \
test.vcf  \                            `# input vcf (positional)` \
test.cram                              `# input reads (positional)`
```
