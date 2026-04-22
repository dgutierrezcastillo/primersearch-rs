import time
import subprocess
import random
import os
import matplotlib.pyplot as plt

# Paths
EMBOSS_PS = "/usr/bin/primersearch"
RUST_PS = "/home/diego/primersearch-rs/target/release/primersearch-rs"

# Ensure Rust is built in release mode
print("Building Rust primersearch in release mode...")
subprocess.run([". $HOME/.cargo/env && cd /home/diego/primersearch-rs && cargo build --release"], shell=True, check=True)

def generate_data(seq_len=1000000, num_primers=100):
    print(f"Generating {seq_len}bp sequence and {num_primers} primers...")
    bases = ['A', 'C', 'G', 'T']
    seq = "".join(random.choice(bases) for _ in range(seq_len))
    
    with open("bench.fasta", "w") as f:
        f.write(">bench_seq\n")
        f.write(seq + "\n")
        
    primers = []
    for i in range(num_primers):
        # Pick a random 20bp window from the sequence to ensure some hits
        start = random.randint(0, seq_len - 500)
        fwd = seq[start:start+20]
        # Pick a reverse primer some distance away
        dist = random.randint(100, 500)
        rev_raw = seq[start+dist:start+dist+20]
        primers.append((f"p{i}", fwd, rev_raw))
        
    with open("bench.primers", "w") as f:
        for p in primers:
            f.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")

def run_bench(cmd, name):
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end = time.time()
    if result.returncode != 0:
        print(f"Error running {name}: {result.stderr}")
    return end - start

def main():
    # Larger benchmark to see real differences
    seq_size = 10000000 # 10MB
    num_p = 200
    generate_data(seq_size, num_p) 
    
    results = {"EMBOSS": [], "Rust": []}
    iterations = 3
    
    print(f"Running benchmark on {seq_size/1e6}MB sequence with {num_p} primers ({iterations} iterations)...")
    
    for i in range(iterations):
        print(f"Iteration {i+1}...")
        
        # EMBOSS command
        emboss_cmd = [EMBOSS_PS, "-seqall", "bench.fasta", "-infile", "bench.primers", "-outfile", "emboss.out", "-mismatchpercent", "0"]
        t_emboss = run_bench(emboss_cmd, "EMBOSS")
        results["EMBOSS"].append(t_emboss)
        
        # Rust command
        rust_cmd = [RUST_PS, "-s", "bench.fasta", "-i", "bench.primers", "-o", "rust.out", "-m", "0"]
        t_rust = run_bench(rust_cmd, "Rust")
        results["Rust"].append(t_rust)
        
    avg_emboss = sum(results["EMBOSS"]) / iterations
    avg_rust = sum(results["Rust"]) / iterations
    
    print(f"\nResults (Average of {iterations} runs):")
    print(f"EMBOSS: {avg_emboss:.4f}s")
    print(f"Rust:   {avg_rust:.4f}s")
    if avg_rust > 0:
        print(f"Speedup: {avg_emboss/avg_rust:.2f}x")
    
    print("\nBenchmark completed successfully.")

if __name__ == "__main__":
    main()
