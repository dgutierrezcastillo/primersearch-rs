use anyhow::{Context, Result};
use clap::Parser;
use needletail::parse_fastx_file;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use rayon::prelude::*;

#[derive(Parser, Debug)]
#[command(author, version, about = "A Rust implementation of EMBOSS primersearch")]
struct Args {
    /// Input FASTA file
    #[arg(short, long)]
    seqall: PathBuf,

    /// Input primer file (Name\tForward\tReverse)
    #[arg(short, long)]
    infile: PathBuf,

    /// Output file
    #[arg(short, long)]
    outfile: PathBuf,

    /// Maximum mismatch percentage (0-100)
    #[arg(short, long, default_value_t = 0)]
    mismatchpercent: u8,
}

#[derive(Debug, Clone)]
struct PrimerPair {
    name: String,
    forward: Vec<u8>,
    reverse: Vec<u8>,
    reverse_rc: Vec<u8>,
}

struct SequenceRecord {
    id: String,
    seq: Vec<u8>,
}

struct Amplimer {
    seq_id: String,
    start: usize,
    end: usize,
    length: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // 1. Read primers
    let primers = read_primers(&args.infile)?;
    println!("Read {} primer pairs", primers.len());

    // 2. Read all sequences into memory
    println!("Loading sequences into memory...");
    let mut sequences = Vec::new();
    let mut reader = parse_fastx_file(&args.seqall).context("Failed to parse FASTA")?;
    while let Some(record) = reader.next() {
        let seq_record = record.context("Invalid FASTA record")?;
        sequences.push(SequenceRecord {
            id: String::from_utf8_lossy(seq_record.id()).to_string(),
            seq: seq_record.seq().to_ascii_uppercase(),
        });
    }

    // 3. Process primers in parallel
    println!("Searching for primers in parallel...");
    let results: Vec<(&PrimerPair, Vec<Amplimer>)> = primers.par_iter().map(|primer| {
        let mut amplimers = Vec::new();
        for seq_rec in &sequences {
            let matches = find_amplimers(&seq_rec.seq, primer, args.mismatchpercent, &seq_rec.id);
            amplimers.extend(matches);
        }
        (primer, amplimers)
    }).collect();

    // 4. Write output sequentially
    let mut out_file = File::create(&args.outfile)?;
    for (primer, amplimers) in results {
        write_emboss_output(&mut out_file, primer, &amplimers)?;
    }

    Ok(())
}

fn read_primers(path: &PathBuf) -> Result<Vec<PrimerPair>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut primers = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            let fwd = parts[1].to_ascii_uppercase().as_bytes().to_vec();
            let rev = parts[2].to_ascii_uppercase().as_bytes().to_vec();
            let reverse_rc = reverse_complement(&rev);
            primers.push(PrimerPair {
                name: parts[0].to_string(),
                forward: fwd,
                reverse: rev,
                reverse_rc,
            });
        }
    }
    Ok(primers)
}

fn find_amplimers(sequence: &[u8], primer: &PrimerPair, mismatch_pct: u8, seq_id: &str) -> Vec<Amplimer> {
    let mut results = Vec::new();
    
    let fwd_matches = find_matches_hamming(sequence, &primer.forward, mismatch_pct);
    let rev_rc_matches = find_matches_hamming(sequence, &primer.reverse_rc, mismatch_pct);

    for &f_start in &fwd_matches {
        for &r_start in &rev_rc_matches {
            if r_start > f_start {
                let end = r_start + primer.reverse.len();
                results.push(Amplimer {
                    seq_id: seq_id.to_string(),
                    start: f_start + 1, // 1-based
                    end: end,
                    length: end - f_start,
                });
            }
        }
    }

    results
}

fn find_matches_hamming(seq: &[u8], pattern: &[u8], mismatch_pct: u8) -> Vec<usize> {
    let mut matches = Vec::new();
    if pattern.is_empty() || seq.len() < pattern.len() {
        return matches;
    }

    let max_dist = (pattern.len() * mismatch_pct as usize) / 100;

    for i in 0..=(seq.len() - pattern.len()) {
        let mut dist = 0;
        let chunk = &seq[i..i + pattern.len()];
        let mut match_found = true;
        for j in 0..pattern.len() {
            if chunk[j] != pattern[j] {
                dist += 1;
                if dist > max_dist {
                    match_found = false;
                    break;
                }
            }
        }
        if match_found {
            matches.push(i);
        }
    }
    matches
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b,
    }).collect()
}

fn write_emboss_output(file: &mut File, primer: &PrimerPair, amplimers: &[Amplimer]) -> Result<()> {
    writeln!(file, "Primer name {}", primer.name)?;
    let primer_name_str = &primer.name;
    for (i, amp) in amplimers.iter().enumerate() {
        writeln!(file, "Amplimer {}", i + 1)?;
        writeln!(file, "\tSequence: {}", amp.seq_id)?;
        writeln!(file, "\t{} hits forward strand at {} with 0 mismatches", primer_name_str, amp.start)?;
        writeln!(file, "\t{} hits reverse strand at {} with 0 mismatches", primer_name_str, amp.end)?;
        writeln!(file, "\tAmplimer length: {} bp", amp.length)?;
        writeln!(file, "")?;
    }
    Ok(())
}
