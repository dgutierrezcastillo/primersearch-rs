import time
import subprocess
import random
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Paths
EMBOSS_PS = "/usr/bin/primersearch"
RUST_PS = "/home/diego/primersearch-rs/target/release/primersearch-rs"
GENOME_FASTA = "arabidopsis.fa"

def get_chromosome_1(input_fasta):
    print("Extracting Chromosome 1...")
    chr1_seq = ""
    with open(input_fasta, 'r') as f:
        is_chr1 = False
        for line in f:
            if line.startswith(">1 "):
                is_chr1 = True
                continue
            elif line.startswith(">") and is_chr1:
                break
            if is_chr1:
                chr1_seq += line.strip()
    return chr1_seq

def run_bench(cmd):
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end = time.time()
    if result.returncode != 0:
        return None, result.stderr
    return end - start, None

def main():
    chr1 = get_chromosome_1(GENOME_FASTA)
    chr1_len = len(chr1)
    print(f"Total Chromosome 1 length: {chr1_len} bp")

    # Sizes to test
    test_sizes = [1000, 10000, 100000, 1000000, 5000000, 10000000, 20000000, chr1_len]
    size_labels = ["1KB", "10KB", "100KB", "1MB", "5MB", "10MB", "20MB", "Chr1 (30MB)"]
    
    num_primers = 50 # Keep primers low so EMBOSS can finish
    
    # Pre-generate 50 primers from the start of the chromosome
    print(f"Generating {num_primers} primers...")
    primers = []
    seed_seq = chr1[:2000] # Use the first 2kb to ensure we have hits in all test sizes
    for i in range(num_primers):
        start = random.randint(0, 1000)
        fwd = seed_seq[start:start+20].upper()
        dist = random.randint(100, 500)
        rev = seed_seq[start+dist:start+dist+20].upper()
        primers.append((f"p{i}", fwd, rev))
        
    with open("scaling.primers", "w") as f:
        for p in primers:
            f.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")

    emboss_times = []
    rust_times = []
    
    print("Building Rust tool...")
    subprocess.run([". $HOME/.cargo/env && cd /home/diego/primersearch-rs && cargo build --release"], shell=True, check=True)

    for size in test_sizes:
        print(f"Testing sequence size: {size} bp...")
        
        # Create temp fasta for this size
        temp_fasta = "temp_scale.fa"
        with open(temp_fasta, "w") as f:
            f.write(f">temp_{size}\n")
            f.write(chr1[:size] + "\n")
            
        # Rust
        print("  Rust...")
        t_rust, err = run_bench([RUST_PS, "-s", temp_fasta, "-i", "scaling.primers", "-o", "rust_scale.out", "-m", "0"])
        rust_times.append(t_rust)
        
        # EMBOSS
        print("  EMBOSS...")
        t_emboss, err = run_bench([EMBOSS_PS, "-seqall", temp_fasta, "-infile", "scaling.primers", "-outfile", "emboss_scale.out", "-mismatchpercent", "0"])
        emboss_times.append(t_emboss)
        
        print(f"  Results: Rust={t_rust:.4f}s, EMBOSS={t_emboss:.4f}s")

    # Plotting
    plt.figure(figsize=(12, 7))
    x = np.arange(len(size_labels))
    width = 0.35

    plt.bar(x - width/2, emboss_times, width, label='EMBOSS', color='#3498db')
    plt.bar(x + width/2, rust_times, width, label='Rust', color='#e67e22')

    plt.xlabel('Sequence Length')
    plt.ylabel('Execution Time (seconds)')
    plt.title(f'Scaling Comparison: EMBOSS vs Rust Primer Search\n({num_primers} Primer Pairs)')
    plt.xticks(x, size_labels)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Log scale for Y because differences are huge
    plt.yscale('log')
    plt.ylabel('Execution Time (seconds, log scale)')

    plt.savefig("scaling_benchmark.png")
    print("\nBenchmark completed. Scaling plot saved to 'scaling_benchmark.png'")

if __name__ == "__main__":
    main()
