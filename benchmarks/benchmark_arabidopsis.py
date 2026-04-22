import time
import subprocess
import random
import os
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend
import matplotlib.pyplot as plt

# Paths
EMBOSS_PS = "/usr/bin/primersearch"
RUST_PS = "/home/diego/primersearch-rs/target/release/primersearch-rs"
GENOME_FASTA = "arabidopsis_subset.fa"

def prepare_subset(input_fasta, output_fasta, size_mb=10):
    print(f"Creating a {size_mb}MB subset of the genome for benchmarking...")
    with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        total_size = 0
        for line in fin:
            fout.write(line)
            if not line.startswith(">"):
                total_size += len(line.strip())
            if total_size > size_mb * 1024 * 1024:
                break

def generate_primers_from_genome(fasta_path, num_primers=100):
    print(f"Extracting {num_primers} realistic primers from {fasta_path}...")
    seq = ""
    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    
    primers = []
    attempts = 0
    while len(primers) < num_primers and attempts < num_primers * 10:
        attempts += 1
        start = random.randint(0, len(seq) - 1000)
        fwd = seq[start:start+20].upper()
        dist = random.randint(100, 500)
        rev = seq[start+dist:start+dist+20].upper()
        if 'N' not in fwd and 'N' not in rev and len(fwd)==20 and len(rev)==20:
            primers.append((f"p{len(primers)}", fwd, rev))

    with open("arabidopsis.primers", "w") as f:
        for p in primers:
            f.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")

def run_bench(cmd):
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end = time.time()
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    return end - start

def main():
    if not os.path.exists("arabidopsis_subset.fa"):
        prepare_subset("arabidopsis.fa", "arabidopsis_subset.fa", 10)
    
    num_p = 100
    generate_primers_from_genome(GENOME_FASTA, num_primers=num_p)
    
    # Ensure Rust is built
    print("Building Rust tool...")
    subprocess.run([". $HOME/.cargo/env && cd /home/diego/primersearch-rs && cargo build --release"], shell=True, check=True)

    results = {"EMBOSS": [], "Rust": []}
    iterations = 1
    
    print(f"Starting Benchmark on Arabidopsis Subset (10MB) with {num_p} primers...")
    
    for i in range(iterations):
        print(f"Iteration {i+1}...")
        
        # Rust first
        print("  Running Rust...")
        rust_cmd = [RUST_PS, "-s", GENOME_FASTA, "-i", "arabidopsis.primers", "-o", "rust_ara.out", "-m", "0"]
        t_rust = run_bench(rust_cmd)
        results["Rust"].append(t_rust)
        print(f"  Rust took: {t_rust:.2f}s")
        
        # EMBOSS second
        print("  Running EMBOSS...")
        emboss_cmd = [EMBOSS_PS, "-seqall", GENOME_FASTA, "-infile", "arabidopsis.primers", "-outfile", "emboss_ara.out", "-mismatchpercent", "0"]
        t_emboss = run_bench(emboss_cmd)
        results["EMBOSS"].append(t_emboss)
        print(f"  EMBOSS took: {t_emboss:.2f}s")
        
    avg_emboss = sum(results["EMBOSS"]) / iterations
    avg_rust = sum(results["Rust"]) / iterations
    
    print(f"\nFinal Results:")
    print(f"EMBOSS Average: {avg_emboss:.2f}s")
    print(f"Rust Average:   {avg_rust:.2f}s")
    print(f"Speedup: {avg_emboss/avg_rust:.1f}x")
    
    # Create Plot
    plt.figure(figsize=(10, 6))
    bars = plt.bar(["EMBOSS", "Rust"], [avg_emboss, avg_rust], color=['#3498db', '#e67e22'])
    plt.ylabel("Time (seconds)")
    plt.title("Primer Search Performance: Arabidopsis Subset (10MB)\n100 Primer Pairs")
    
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.1, f"{yval:.2f}s", ha='center', va='bottom', fontsize=12, fontweight='bold')
        
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig("arabidopsis_benchmark.png")
    print("\nBenchmark plot saved to 'arabidopsis_benchmark.png'")

if __name__ == "__main__":
    main()
