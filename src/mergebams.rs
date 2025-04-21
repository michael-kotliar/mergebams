extern crate simple_log;
extern crate clap;
extern crate bam;
extern crate csv;
extern crate flate2;
extern crate differ;
extern crate itertools;

use itertools::Itertools;
use differ::{Differ, Tag};
use bam::RecordWriter;
use bam::record::tags::TagValue;
use clap::{App, load_yaml};
use std::str;
use std::process::Command;


struct Params {
    inputs: String,
    output: String,
    threads: usize,
}

fn main() {
    let params = load_params();
    if let Ok((header, params)) = checkheaders(params){
        process_bams(params, header);
        return;
    }else{
        eprintln!("BAM header sequences do not match. Exiting.");
    }
}

fn sanitize_and_replace_cb(cb: &str, suffix: &str) -> String {
    let re = regex::Regex::new(r"-\d+$").unwrap();
    if re.is_match(cb) {
        re.replace(cb, suffix).to_string()
    } else {
        format!("{}{}", cb, suffix)
    }
}

fn load_params() -> Params {
    let yaml = load_yaml!("params_mergebams.yml");
    let params = App::from_yaml(yaml).get_matches();
    let inputs = params.value_of("inputs").unwrap();
    let threads: usize = params.value_of("threads").unwrap_or("1").parse().unwrap();
    let output = params.value_of("output").unwrap();
    Params{
        inputs: inputs.to_string(),
        output: output.to_string(),
        threads: threads
    }
}

fn checkheaders(params: Params) -> Result<(bam::Header, Params), &'static str>{
    let inputs = params.inputs.to_string();
    let bam_vec = inputs.split(",").collect::<Vec<&str>>();
    let bam_vec_to_header = inputs.split(",").collect::<Vec<&str>>();
    let bam_it = bam_vec.into_iter().combinations(2);
    let mut grand_count = 0;
    for bam_pairs in bam_it {
        let mut header_vec = Vec::new();
        for inbam in bam_pairs.iter() {
            let hreader = bam::BamReader::from_path(inbam, 0).unwrap();
            let header = hreader.header().clone();
            header_vec.push(header);
        }
        let mut count = 0;
        let a = header_vec[0].reference_names();
        let b = header_vec[1].reference_names();
        let differ = Differ::new(&a, &b);
        for span in differ.spans() {
            match span.tag {
                Tag::Equal => (),
                _ => count+=1,
            }
        }
        eprintln!("Found {} discrepencies between sequence names in the header of:\n\n{}\nand:\n\n{}\n", count, bam_pairs[0], bam_pairs[1]);
        grand_count+=count;
    }
    if grand_count == 0{
        let new_header = make_new_header(bam_vec_to_header);
        return Ok((new_header, params));
    } else {
        return Err("Discrepent headers")
    }
}

fn make_new_header(bam_vec: Vec<&str>) -> bam::Header {
    let hreader = bam::BamReader::from_path(bam_vec[0], 0).unwrap();
    let mut out_header = hreader.header().clone();
    let mergebam_line = bam_vec.join(", ");
    let _msg = out_header.push_line(&("@CO\tmergebams has included the BAM records from the following files (using the header from the first): ".to_owned()+&mergebam_line));
    return out_header;
}

fn process_bams(params: Params, header: bam::Header) -> Params{
    let out_bam = "_temp.bam";
    let inputs = params.inputs.to_string();
    let bam_vec = inputs.split(",").collect::<Vec<&str>>();
    let (read_threads, write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let mut fail_count = 0;
    let mut pass_count = 0;
    let mut other_count = 0;
    let mut pass_writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(write_threads)
        .from_path(out_bam, header.clone()).unwrap();

    for (pos, inbam) in bam_vec.iter().enumerate() {
        let reader = bam::BamReader::from_path(inbam.to_string(), read_threads).unwrap();
        let suffix = format!("-{}", pos + 1);
        for record in reader {
            let mut newrecord = record.as_ref().unwrap().clone();
            match record.unwrap().tags().get(b"CB") {
                Some(TagValue::String(array_view, _)) => {
                    let cb_str = String::from_utf8_lossy(array_view.as_ref()).to_string();
                    let new_cb = sanitize_and_replace_cb(&cb_str, &suffix);
                    newrecord.tags_mut().remove(b"CB");
                    newrecord.tags_mut().push_string(b"CB", new_cb.as_bytes());
                    pass_writer.write(&newrecord).unwrap();
                    pass_count+=1;
                },
                Some(TagValue::Char(value)) => {
                    let cb_str = value.to_string();
                    let new_cb = sanitize_and_replace_cb(&cb_str, &suffix);
                    newrecord.tags_mut().remove(b"CB");
                    newrecord.tags_mut().push_string(b"CB", new_cb.as_bytes());
                    pass_writer.write(&newrecord).unwrap();
                    other_count+=1;
                },
                _ => {
                    fail_count+=1;
                }
            }
        }
    }
    pass_writer.finish()?;
    eprintln!("Processed all reads!!\nFound:\n{} - reads PASSING\n{} - reads PASSING but with issues\n{} - reads FAILING", pass_count, other_count, fail_count);

    let status = Command::new("samtools")
        .args(&["sort", "-@", &params.threads.to_string(), "-o", &params.output.to_string(), "_temp.bam"])
        .status()
        .expect("samtools not found");
    if !status.success() {
        panic!("samtools sort failed");
    }

    let status = Command::new("samtools")
        .args(&["index", "-@", &params.threads.to_string(), &params.output.to_string()])
        .status()
        .expect("samtools not found");
    if !status.success() {
        panic!("samtools index failed");
    }

    return params;
}