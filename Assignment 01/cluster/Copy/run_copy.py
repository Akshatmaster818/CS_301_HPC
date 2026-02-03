import os
import csv
import sys
import subprocess
from collections import defaultdict

# User-configurable parameters
EXECUTABLE = "./output"        # compiled binary file
NUM_RUNS = 10                         # number of repetitions
OUTPUT_CSV = "results/avg_algo_times.csv"

# Ensure results directory exists
os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

# Storage for results, Keyed by problem size
results = defaultdict(lambda: {
    "runs": None,
    "total_particles": None,
    "e2e_times": [],
    "algo_times": []
})

print(f"Running benchmark {NUM_RUNS} times...\n")

# Run executable multiple times
for run_id in range(NUM_RUNS):
    print(f"=== Run {run_id + 1}/{NUM_RUNS} ===")

    proc = subprocess.run(
        [EXECUTABLE],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )

    lines = proc.stdout.strip().split("\n")

    # Print raw output for this run
    for line in lines:
        print(line)

    # Skip header, parse data
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(",")]
        if len(parts) != 5:
            continue

        problem_size = int(parts[0])
        runs = int(parts[1])
        total_particles = int(parts[2])
        e2e_time = float(parts[3])
        algo_time = float(parts[4])

        results[problem_size]["runs"] = runs
        results[problem_size]["total_particles"] = total_particles
        results[problem_size]["e2e_times"].append(e2e_time)
        results[problem_size]["algo_times"].append(algo_time)

    print()

# Calculate averages and write to CSV
print(f"Writing results to {OUTPUT_CSV}\n")

with open(OUTPUT_CSV, mode='w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Np-ProblemSize", "RUNS", "TotalParticles", "AvgE2ETime", "AvgAlgoTime"])

    for problem_size in sorted(results.keys()):
        data = results[problem_size]
        avg_e2e_time = sum(data["e2e_times"]) / len(data["e2e_times"])
        avg_algo_time = sum(data["algo_times"]) / len(data["algo_times"])

        writer.writerow([
            problem_size,
            data["runs"],
            data["total_particles"],
            avg_e2e_time,
            avg_algo_time
        ])

print("Done!")
