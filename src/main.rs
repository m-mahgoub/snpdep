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

fn validate_vcf_file(val: &str) -> Result<String, String> {
    if val.to_lowercase().ends_with(".vcf") | val.to_lowercase().ends_with("vcf.gz") {
        Ok(val.to_string())
    } else {
        Err("VCF File must end with '.vcf' or '.vcf.gz'".to_string())
    }
}

fn validate_format_id(val: &str) -> Result<String, String> {
    if val.len() <= 5 {
        Ok(val.to_string())
    } else {
        Err("Format ID length must not exceed 5 characters.".to_string())
    }
}

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

#[derive(Parser)]
#[command(
    author,
    version,
    about = "Program to calculate haplotype-specific depth from BAM/CRAM files at positions provided in a genotype VCF file.",
    arg_required_else_help(true)
)]
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
enum RecordData {
    Success {
        chr: String,
        pos: i64,
        ref_allele_str: String,
        alt_allele_str: String,
    },
}

fn process_record(record: &rust_htslib::bcf::record::Record) -> Result<RecordData, &'static str> {
    // Get record chromosome
    let chr = match record.rid() {
        Some(rid) => {
            let record_header = record.header();
            let seq = record_header.rid2name(rid).unwrap();
            str::from_utf8(seq).unwrap().to_owned()
        }
        None => "UNKNOWN".to_owned(),
    };

    // Get record position
    let pos = record.pos();
    // binding suggested by compiler to avoid short living variables
    let alleles_binding = record.alleles();
    let ref_allele = alleles_binding.get(0).unwrap();
    let alt_allele = alleles_binding.get(1).unwrap();
    let ref_allele_str = str::from_utf8(ref_allele)
        .map_err(|_| "Invalid UTF-8 in reference allele")?
        .to_owned();
    let alt_allele_str = str::from_utf8(alt_allele)
        .map_err(|_| "Invalid UTF-8 in alternate allele")?
        .to_owned();

    Ok(RecordData::Success {
        chr,
        pos,
        ref_allele_str,
        alt_allele_str,
    })
}

fn process_bam_data(
    bam: &mut rust_htslib::bam::IndexedReader,
    chr: &str,
    pos: i64,
    ref_allele_str: &str,
    alt_allele_str: &str,
    record: &mut rust_htslib::bcf::record::Record,
    chunk_vcf: &mut rust_htslib::bcf::Writer,
    format_id: &str,
    min_mapq: u8,
    min_count: i32,
) -> Result<rust_htslib::bcf::record::Record, &'static str> {
    let _ = bam.fetch((chr, pos, pos + 1));
    let mut ref_count = 0;
    let mut alt_count = 0;
    for p in bam.pileup() {
        let pileup = p.unwrap();
        if i64::from(pileup.pos()) == pos {
            for alignment in pileup.alignments() {
                if !alignment.is_del()
                    && !alignment.is_refskip()
                    && !alignment.record().is_duplicate()
                    && alignment.record().mapq() >= min_mapq
                {
                    let bam_record = alignment.record();
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
    if ref_count + alt_count >= min_count {
        let flattened_array: Vec<i32> =
            vec![[ref_count, alt_count]].into_iter().flatten().collect();
        chunk_vcf.translate(record);
        let _ = record.push_format_integer(format_id.as_bytes(), &flattened_array);
        // chunk_vcf.write(record).unwrap();
    }
    let new_record = record.clone();

    Ok(new_record)
}

fn main() {
    let cli = Cli::parse();
    let bcf_path = &cli.input_vcf.unwrap();
    let bam_path = &cli.reads.unwrap();
    let out_vcf_path = cli.output_vcf.unwrap();
    // make sure to not override the input file by mistake!
    if bcf_path.as_str() == out_vcf_path.as_str() {
        eprintln!("Input and Output VCF paths are the same! Use different paths");
        std::process::exit(0)
    };
    // set output compression
    let is_out_vcf_uncompressed: bool = out_vcf_path.to_lowercase().ends_with(".vcf");
    let reads_format = cli.reads_format.unwrap();
    // validate cram input
    if bam_path.ends_with("cram") && reads_format.as_str() == "bam" {
        eprintln!("CRAM file is provided as input! Set foramat input to cram `--reads-format cram` and provide referecne fasta using `--reference <FASTA>`");
        std::process::exit(0)
    }

    #[cfg(not(target_os = "linux"))]
    if reads_format.as_str() == "cram" {
        eprintln!("CRAM format is supported only in Linux OS ");
        std::process::exit(0)
    }

    let mut ref_path = match reads_format.as_str() {
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
    // Configure the rayon thread pool with the desired number of threads
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();
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
    println!("");
    println!("Using {} Threads", num_threads);
    println!(
        "Processing {} Records: split into {} Record chunks...",
        records.len(),
        chunk_size,
    );
    println!("");
    let mpb = MultiProgress::new();
    mpb.set_draw_target(ProgressDrawTarget::stdout());
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] [{bar:60.cyan/blue}] {pos:>6}/{len:6} {msg}",
    )
    .unwrap()
    .progress_chars("##-");
    let pb1 = ProgressBar::new(records.len().try_into().unwrap());
    pb1.set_style(sty.clone());
    pb1.set_message("Annotating Records");
    let all_chunks_map: Arc<
        Mutex<BTreeMap<usize, BTreeMap<usize, rust_htslib::bcf::record::Record>>>,
    > = Arc::new(Mutex::new(BTreeMap::new()));
    let chunks = records.chunks_mut(chunk_size);
    let mut records_count = Arc::new(Mutex::new(0));

    // Iteration over vcf records
    chunks
        .collect::<Vec<_>>()
        .par_iter_mut()
        .enumerate()
        .for_each(|(chunk_num, chunk)| {
            let chunk_header = chunk.first().as_ref().unwrap().as_ref().unwrap().header(); // get header from first item in the chunk
            let mut new_bcf_header = bcfHeader::from_template(chunk_header);
            // let format_rc_line = new_format_line;
            new_bcf_header.push_record(new_format_line.as_bytes());
            let mut chunk_vcf =
                bcfWriter::from_path("/dev/null", &new_bcf_header, true, Vcf).unwrap();
            let mut bam_reader =
                bamIndexedReader::from_path(bam_path).expect("Failed to read bam file");
            if reads_format == "cram" {
                bam_reader.set_reference(ref_path.to_string());
            }
            let mut chunk_records = 0;
            let mut chunk_map: BTreeMap<usize, rust_htslib::bcf::record::Record> = BTreeMap::new();
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
                        // Here, you have the processed data for each record.
                        // You can perform any further processing or actions you want with this data.
                        // For example, printing the processed data:
                        // println!(
                        //     "{:?}:{:?} {:?}>{:?} , Thread: {:?}",
                        //     chr,
                        //     pos,
                        //     ref_allele_str,
                        //     alt_allele_str,
                        //     thread::current().id()
                        // );
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
                                // chunk_vcf.write(record).unwrap();
                                // ... (rest of the code for processing the record data)
                                // Here you can use chr, pos, ref_allele_str, and alt_allele_str as needed
                            }
                            Err(err) => {
                                // Handle the error case
                                eprintln!("Error processing BAM data for record {}: {}", i, err);
                                // continue;
                            }
                        }
                    }
                    Err(e) => {
                        // Handle the error if the processing fails (process_record).
                        // You can choose to print an error message or perform other actions.
                        eprintln!("Error processing record: {}", e);
                    }
                }
                if chunk_records == chunk_size {
                    let mut records_count_lock = records_count.lock().unwrap();
                    *records_count_lock += chunk_size;
                    if mpb.is_hidden() {
                        println!(
                            "{} Records Processed",
                            records_count_lock.separated_string()
                        );
                    }

                    chunk_records = 0;
                }
            }
            // let mut all_chunks_map_lock;
            if !chunk_map.is_empty() {
                let mut all_chunks_map_lock = all_chunks_map.lock().unwrap();
                all_chunks_map_lock.insert(chunk_num, chunk_map);
            }
            pb1.inc(chunk.len().try_into().unwrap());
        });
    pb1.finish_with_message("Annotation Done!");
    let all_chunks_map_lock = &*all_chunks_map.lock().unwrap();
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
