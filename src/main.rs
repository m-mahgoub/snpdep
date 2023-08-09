#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use clap::Parser;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::{prelude::*, ThreadPoolBuilder};
use rust_htslib::bam::{
    ext::BamRecordExtensions, IndexedReader as bamIndexedReader, Read as bamRead,
};
use rust_htslib::bcf::{
    header::Header as bcfHeader, Format::Vcf, Read as bcfRead, Reader as bcfReader,
    Writer as bcfWriter,
};

use separator::Separatable;
use std::{collections::BTreeMap, process, str, sync::Arc, sync::Mutex};

// Define CLI validation Functions:
// Function to validate bam/cram files from the cli arguments. CRAM files in MacOS will not be accepted
fn validate_bam_file(val: &str) -> Result<String, String> {
    #[cfg(target_os = "linux")]
    if val.to_lowercase().ends_with(".bam") | val.to_lowercase().ends_with(".cram") {
        Ok(val.to_string())
    } else {
        Err("Reads File must end with '.bam' or '.cram'".to_string())
    }
    #[cfg(not(target_os = "linux"))]
    if val.to_lowercase().ends_with(".bam") {
        Ok(val.to_string())
    } else if val.to_lowercase().ends_with(".cram") {
        Err("CRAM files are supported only in Linux OS".to_string())
    } else {
        Err("Reads File must end with '.bam'".to_string())
    }
}

// Function to validate VCF file from the cli arguments.
fn validate_vcf_file(val: &str) -> Result<String, String> {
    if val.to_lowercase().ends_with(".vcf") | val.to_lowercase().ends_with("vcf.gz") {
        Ok(val.to_string())
    } else {
        Err("VCF File must end with '.vcf' or '.vcf.gz'".to_string())
    }
}

// Function to validate format id from the cli arguments to ensure it is not more than 5 characters. Long characters can take a lot of space in the new annotated vcf output.
fn validate_format_id(val: &str) -> Result<String, String> {
    if val.len() <= 5 {
        Ok(val.to_string())
    } else {
        Err("Format ID length must not exceed 5 characters.".to_string())
    }
}
#[derive(Parser)]
#[command(
    author,
    version,
    about = "Program to calculate haplotype-specific depth from BAM/CRAM files at positions provided in a genotype VCF file.",
    arg_required_else_help(true)
)]
// Define CLI arguments
struct Cli {
    /// Path to VCF file containing sample genotype. Both uncompressed `.vcf` and compressed `.vcf.gz` files are supported
    #[arg(value_name = "INPUT_VCF", required = true, index = 1, value_parser=clap::builder::ValueParser::new(validate_vcf_file))]
    input_vcf: Option<String>,
    /// Path to BAM/CRAM reads file.
    #[arg(value_name = "READS",required=true, index=2,value_parser=clap::builder::ValueParser::new(validate_bam_file))]
    reads: Option<String>,
    /// Path to the output annotated VCF file. Compression will be inferred from the file extension: `.vcf` for uncompressed and `.vcf.gz` for compressed output.
    #[arg(
        short = 'o',
        long = "output",
        value_name = "OUTPUT_VCF",
        required = true,
        value_parser=clap::builder::ValueParser::new(validate_vcf_file)
    )]
    output_vcf: Option<String>,
    /// Specifies the reads file format.
    #[arg(short='f', long="reads-format", value_name = "READS_FORMAT", value_parser=["bam","cram"], default_value="bam")]
    reads_format: Option<String>,
    /// Path to the reference FASTA file (required if `--reads-format cram`).
    #[arg(
        short,
        long,
        value_name = "FASTA",
        required_if_eq("reads_format", "cram")
    )]
    reference: Option<String>,
    /// Minimum read mapping quality (MAPQ) filtration threshold.
    #[arg(long = "min-mapq", value_name = "MIN_MAPQ", default_value = "1")]
    min_mapq: Option<u8>,
    /// Minimum read coverage required at a position to be annotated. if set to 0, all positions will be annotated.
    #[arg(long = "min-count", value_name = "MIN_COUNT", default_value = "1")]
    min_count: Option<i32>,
    /// Number of SNPs processed by a thread in each iteration.
    #[arg(short, long, value_name = "CHUNK_SIZE", default_value = "1000")]
    chunksize: Option<usize>,
    /// Unique ID for the new VCF format field (e.g RNA, K4me3), max 5 characters.
    #[arg(long="format-field-id", value_name = "FORMAT_ID", default_value = "RC",value_parser=clap::builder::ValueParser::new(validate_format_id))]
    format_id: Option<String>,
    /// Short description for new VCF format field (e.g RNAseq, H3K4me3).
    #[arg(
        long = "format-field-name",
        value_name = "FORMAT_NAME",
        default_value = "BAM/CRAM"
    )]
    format_name: Option<String>,

    /// Number of threads to use.
    #[arg(short, long, value_name = "NUM_THREADS", default_value = "1")]
    threads: Option<usize>,
}

// Define Functions provoked in main()
// Function to validate that user provided format id is not already existing in the input vcf file.
// If provided id already exits the program will exit.
fn validate_unique_format_id(
    bcf_header_view: &rust_htslib::bcf::header::HeaderView,
    format_id: &str,
) {
    for header_record in bcf_header_view.header_records() {
        match header_record {
            rust_htslib::bcf::header::HeaderRecord::Format { key, values } => {
                for (tag, value) in values.iter() {
                    if tag == "ID" && value == format_id {
                        eprintln!("{} is already used as a unique Format ID in the provided VCF. Provide a unique Format ID using `--format-field-id <FORMAT_ID>`", format_id);
                        std::process::exit(0);
                    } else {
                        continue;
                    }
                }
            }
            _ => {
                // Continue to the next record if it's not of type HeaderRecord::Format
                continue;
            }
        }
    }
}

// Define enum RecordData: this is just re-formatting a snp record.
// RecordData the output format from process_record() function, and will be passed in this new format to process_bam_data() function
enum RecordData {
    Success {
        chr: String,
        pos: i64,
        ref_allele_str: String,
        alt_allele_str: String,
    },
}

// function that reformat snp record from `rust_htslib::bcf::record::Record` to `RecordData`
fn process_record(record: &rust_htslib::bcf::record::Record) -> Result<RecordData, &'static str> {
    // Get record chromosome
    let chr: String = match record.rid() {
        Some(rid) => {
            let record_header = record.header();
            let seq = record_header.rid2name(rid).unwrap();
            str::from_utf8(seq).unwrap().to_owned()
        }
        None => "UNKNOWN".to_owned(),
    };

    // Get record position
    let pos = record.pos();
    // Get record ref and alt as strings
    let alleles_binding = record.alleles();
    let ref_allele = alleles_binding.get(0).unwrap();
    let alt_allele = alleles_binding.get(1).unwrap();
    let ref_allele_str = str::from_utf8(ref_allele)
        .map_err(|_| "Invalid UTF-8 in reference allele")?
        .to_owned();
    let alt_allele_str = str::from_utf8(alt_allele)
        .map_err(|_| "Invalid UTF-8 in alternate allele")?
        .to_owned();
    // Return result as RecordData enum
    Ok(RecordData::Success {
        chr,
        pos,
        ref_allele_str,
        alt_allele_str,
    })
}

// This is the main logic of the whole analysis.
// It will take (1) a snp record and (2) bam file reader. It will then count the reads supporting the ref and alt alleles.
// The function will return modified record (rust_htslib::bcf::record::Record) by adding new foramt id in the format field and alleles count in the sample field.
fn process_bam_data(
    // indexed bam reader
    bam: &mut rust_htslib::bam::IndexedReader,
    // contig
    chr: &str,
    // snp position
    pos: i64,
    // ref allele
    ref_allele_str: &str,
    // alt allele
    alt_allele_str: &str,
    // record as `rust_htslib::bcf::record::Record` Struct
    record: &mut rust_htslib::bcf::record::Record,
    // out vcf writer defined within the scope of a single thread.
    // It is similar to input vcf writer, only with a modified header (+ new format field header appended)
    chunk_vcf: &mut rust_htslib::bcf::Writer,
    // new format id provided by the cli
    format_id: &str,
    // minimum mapq provided by cli
    min_mapq: u8,
    // minimum reads count to annotate the record, provided by cli
    min_count: i32,
) -> Result<rust_htslib::bcf::record::Record, &'static str> {
    // fetch bam region overlapping the snp
    let _ = bam.fetch((chr, pos, pos + 1));
    // set ref/alt counts to 0
    let mut ref_count = 0;
    let mut alt_count = 0;
    // iterate over pileups in the fetched bam region
    for p in bam.pileup() {
        let pileup = p.unwrap();
        // make sure the pileup overlaps snp position
        if i64::from(pileup.pos()) == pos {
            // iterate over pileup alignments
            for alignment in pileup.alignments() {
                // filter the alignment
                if !alignment.is_del()
                    && !alignment.is_refskip()
                    && !alignment.record().is_duplicate()
                    && alignment.record().mapq() >= min_mapq
                {
                    let bam_record: rust_htslib::bam::Record = alignment.record();
                    // extract the nucleotide in the query position in the current bam record, then add it to the ref/alt counts if matches either.
                    for ref_pos in bam_record.reference_positions() {
                        if ref_pos == pos {
                            let read_pos = alignment.qpos().unwrap();
                            let pos_seq_byte = bam_record.seq()[read_pos];
                            let pos_seq_str =
                                std::str::from_utf8(&[pos_seq_byte]).unwrap().to_string();
                            if pos_seq_str == ref_allele_str {
                                ref_count = ref_count + 1;
                            } else if pos_seq_str == alt_allele_str {
                                alt_count = alt_count + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    // modify the record only if the total counts are more then min_count threshold.
    if ref_count + alt_count >= min_count {
        // reformat the counts to Vec<i32> as expected by push_format_integer() function
        // N.B: Because there is one sample, the flattened_array has one element (length of 1)
        let flattened_array: Vec<i32> =
            vec![[ref_count, alt_count]].into_iter().flatten().collect();
        // Translate record to header of the new output writer, i.e. it will inherit the format fields from the new header
        chunk_vcf.translate(record);
        let _ = record.push_format_integer(format_id.as_bytes(), &flattened_array);
    }
    // return new record
    let new_record = record.clone();
    Ok(new_record)
}

fn main() {
    // retrieve arguments and define conditional exits based on arguments values
    let cli = Cli::parse();
    let bcf_path = &cli.input_vcf.unwrap();
    let bam_path = &cli.reads.unwrap();
    let out_vcf_path = cli.output_vcf.unwrap();
    // make sure to not override the input file by mistake!
    if bcf_path.as_str() == out_vcf_path.as_str() {
        eprintln!("Input and Output VCF paths are the same! Use different paths");
        std::process::exit(0)
    };
    // set whether output should be compressed or not
    let is_out_vcf_uncompressed: bool = out_vcf_path.to_lowercase().ends_with(".vcf");
    let reads_format = cli.reads_format.unwrap();
    // if input in cram, force the usage of --reads-format & --reference
    if bam_path.ends_with("cram") && reads_format.as_str() == "bam" {
        eprintln!("CRAM file is provided as input! Set foramat input to cram `--reads-format cram` and provide referecne fasta using `--reference <FASTA>`");
        std::process::exit(0)
    }
    // if input is cram in MaOS exit. There is a bug that leads to a segmentation error when trying to use cram in MacOS.
    #[cfg(not(target_os = "linux"))]
    if reads_format.as_str() == "cram" {
        eprintln!("CRAM format is supported only in Linux OS ");
        std::process::exit(0)
    }

    let ref_path = match reads_format.as_str() {
        "cram" => cli.reference.unwrap(),
        _ => String::new(),
    };
    let chunk_size = cli.chunksize.unwrap();
    let num_threads = cli.threads.unwrap();
    let min_count = cli.min_count.unwrap();
    let min_mapq = cli.min_mapq.unwrap();
    let format_id = cli.format_id.unwrap();
    let format_name = cli.format_name.unwrap();
    let new_format_line = format!(
        r#"##FORMAT=<ID={},Number=.,Type=Integer,Description="Number of {} Reads Originating from (REF,ALT) Alleles">"#,
        format_id, format_name
    );
    // Define files and readers
    let mut bcf_reader = bcfReader::from_path(bcf_path).unwrap_or_else(|err| {
        eprintln!("{}", err);
        process::exit(1);
    });
    let bcf_header_view = bcf_reader.header();
    // Validate uniqueness of Format ID. It will exit if validation failed
    validate_unique_format_id(&bcf_header_view, &format_id);
    let mut out_vcf_header = bcfHeader::from_template(bcf_header_view);
    // confirm the new header ID is not already present
    out_vcf_header.push_record(new_format_line.as_bytes());
    let mut out_vcf =
        bcfWriter::from_path(out_vcf_path, &out_vcf_header, is_out_vcf_uncompressed, Vcf).unwrap();

    let mut records: Vec<_> = bcf_reader.records().collect::<Vec<_>>();
    // Print a start message.
    println!("");
    println!("Using {} Threads", num_threads);
    println!(
        "Processing {} Records: split into {} Record chunks...",
        records.len(),
        chunk_size,
    );
    println!("");
    // Initiate multi-progress bar. This progress bar will be suppressed by default if stdout is redirected to a file (e.g log file)
    let mpb = MultiProgress::new();
    mpb.set_draw_target(ProgressDrawTarget::stdout());
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] [{bar:60.cyan/blue}] {pos:>6}/{len:6} {msg}",
    )
    .unwrap()
    .progress_chars("##-");
    // initiate progress bar (1) for processing the reads
    let pb1 = ProgressBar::new(records.len().try_into().unwrap());
    pb1.set_style(sty.clone());
    pb1.set_message("Annotating Records");

    // Define all_chunks_map which is `std::collections::BTreeMap` data type. It is BTreeMap of nested BTreeMap(s)
    // all_chunks_map is Arc<Mutex<_>> type. This is because it will be shared by multiple threads and should be locked/unlocked to allow exactly one thread to gain access to it at a time
    // In all_chunks_map the keys are chunk numbers and values are BTreeMap of records (keys: record number, values: records)
    // One nested BTreeMap is the result of one thread processing a chunk of records. This nested BTreeMap (records BTreeMap) is defined in the scope of the thread
    // A thread will iterate over records in a chunk and insert each processed record in the records BTreeMap.
    // After a thread consumes all records in a chunk it will:
    // (1) lock all_chunks_map
    // (2) insert the the records BTreeMap in all_chunks_map
    // (3) Unlock all_chunks_map
    // (4) destroy current chunk records BTreeMap
    // (5) move on to process anther chunk
    // BTreeMap was chosen over HashMap because BTreeMap is sorted by its keys. This ensures the the records are sorted in the output vcf file.
    let all_chunks_map: Arc<
        Mutex<BTreeMap<usize, BTreeMap<usize, rust_htslib::bcf::record::Record>>>,
    > = Arc::new(Mutex::new(BTreeMap::new()));
    // split the records into chunks. chunk size is determined by `--chunksize` option
    let chunks = records.chunks_mut(chunk_size);
    // Set records counter to be used in the logging. It is also Arc<Mutex<_>> type because it will be shared between threads
    let records_count = Arc::new(Mutex::new(0));
    // Configure the rayon thread pool with the desired number of threads
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();
    // Iterate over vcf records in chunks. Multithreading is provoked in this loop by par_iter_mut(). Each chunk will be passed to a thread in the above defined ThreadPool
    chunks
        .collect::<Vec<_>>()
        .par_iter_mut()
        .enumerate()
        .for_each(|(chunk_num, chunk)| {
            // define vcf and bam readers within the scope of the thread
            let chunk_header = chunk.first().as_ref().unwrap().as_ref().unwrap().header(); // get header from first item in the chunk
            let mut new_bcf_header = bcfHeader::from_template(chunk_header);
            new_bcf_header.push_record(new_format_line.as_bytes());
            let mut chunk_vcf =
                bcfWriter::from_path("/dev/null", &new_bcf_header, true, Vcf).unwrap();
            let mut bam_reader =
                bamIndexedReader::from_path(bam_path).expect("Failed to read bam file");
            if reads_format == "cram" {
                bam_reader.set_reference(ref_path.to_string());
            }
            // count the records within a chunk, for logging purpose
            let mut chunk_records = 0;
            // This is the inner BTreeMap to be added in all_chunks_map. It is not accessible by multiple threads and therefore not necessarily Arc<Mutex<_>> type
            let mut chunk_map: BTreeMap<usize, rust_htslib::bcf::record::Record> = BTreeMap::new();
            // iterate over each record in the chunk
            for (i, record_result) in chunk.iter_mut().enumerate() {
                chunk_records += 1;
                let mut record = record_result.as_mut().expect("Failed to get Ref allele");
                // Process Record using the process_record function:
                match process_record(&mut record) {
                    Ok(RecordData::Success {
                        chr,
                        pos,
                        ref_allele_str,
                        alt_allele_str,
                    }) => {
                        // This closure for successful process_record() call
                        match process_bam_data(
                            &mut bam_reader,
                            &chr,
                            pos,
                            &ref_allele_str,
                            &alt_allele_str,
                            &mut record,
                            &mut chunk_vcf,
                            format_id.as_str(),
                            min_mapq,
                            min_count,
                        ) {
                            Ok(new_record) => {
                                chunk_map.insert(i, new_record);
                            }
                            Err(err) => {
                                // Handle the error case
                                eprintln!("Error processing BAM data for record {}: {}", i, err);
                            }
                        }
                    }
                    Err(e) => {
                        // This closure for failed process_record() call
                        eprintln!("Error processing record: {}", e);
                    }
                }
                // This closure is to add count to the total number of processed records, and log after finishing  a chunk
                if chunk_records == chunk_size {
                    let mut records_count_lock = records_count.lock().unwrap();
                    *records_count_lock += chunk_size;
                    // print to stout if progress bar is suppressed
                    if mpb.is_hidden() {
                        println!(
                            "{} Records Processed",
                            records_count_lock.separated_string()
                        );
                    }

                    chunk_records = 0;
                }
            }
            // add current chunk BTreeMap to global all_chunks_map
            if !chunk_map.is_empty() {
                let mut all_chunks_map_lock = all_chunks_map.lock().unwrap();
                all_chunks_map_lock.insert(chunk_num, chunk_map);
            }
            // increment records-processing progress bar
            pb1.inc(chunk.len().try_into().unwrap());
        });
    // This is the end of all records processing!
    pb1.finish_with_message("Annotation Done!");
    // now lock the all_chunks_map, and iterate over it to write new records into the new vcf file
    let all_chunks_map_lock = &*all_chunks_map.lock().unwrap();
    // add new progress bar for records' writing
    let pb2 = ProgressBar::new(records.len().try_into().unwrap());
    pb2.set_style(sty.clone());
    pb2.set_message("Writing Records");
    for (chunk_num, chunk_map) in all_chunks_map_lock {
        for (i, record) in chunk_map {
            out_vcf.write(record);
        }
        pb2.inc(chunk_map.len().try_into().unwrap());
    }
    pb2.finish_with_message("Writing Done!");
}
