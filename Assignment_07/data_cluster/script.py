import re
import csv
import os

def parse_benchmark_data(input_file, output_file):
    # Regular expressions for the new data format
    # Matches "250 100 900000:" -> Grid X, Grid Y, Particles
    config_pattern = re.compile(r"^(\d+)\s+(\d+)\s+(\d+):")
    
    # Matches either a standalone number "16" OR a command line prompt "... ./output16 ..."
    # Extracts the thread count
    thread_pattern_simple = re.compile(r"^(\d+)$")
    thread_pattern_cmd = re.compile(r"\./output(\d+)\s")
    
    # Matches the timing and metric lines
    interp_pattern = re.compile(r"Total Interpolation Time = ([0-9.]+)")
    norm_pattern = re.compile(r"Total Normalization Time = ([0-9.]+)")
    mover_pattern = re.compile(r"Total Mover Time = ([0-9.]+)")
    denorm_pattern = re.compile(r"Total Denormalization Time = ([0-9.]+)")
    algo_pattern = re.compile(r"Total Algorithm Time = ([0-9.]+)")
    voids_pattern = re.compile(r"Total Number of Voids = (\d+)")

    parsed_data = []

    # State variables
    current_grid = None
    current_particles = None
    current_threads = None
    
    # Dictionary to temporarily hold the current block's metrics
    current_metrics = {}

    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    with open(input_file, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
                
            # 1. Check for Configuration (Grid & Particles)
            config_match = config_pattern.search(line)
            if config_match:
                # Combine X and Y into "X x Y"
                current_grid = f"{config_match.group(1)} x {config_match.group(2)}"
                current_particles = config_match.group(3)
                continue
            
            # 2. Check for Thread counts
            thread_simple_match = thread_pattern_simple.search(line)
            thread_cmd_match = thread_pattern_cmd.search(line)
            
            if thread_simple_match:
                current_threads = thread_simple_match.group(1)
                current_metrics = {} # Reset metrics for new thread run
                continue
            elif thread_cmd_match:
                current_threads = thread_cmd_match.group(1)
                current_metrics = {} # Reset metrics for new thread run
                continue
                
            # 3. Check for Timings
            if m := interp_pattern.search(line):
                current_metrics['Interpolation Time'] = m.group(1)
            elif m := norm_pattern.search(line):
                current_metrics['Normalization Time'] = m.group(1)
            elif m := mover_pattern.search(line):
                current_metrics['Mover Time'] = m.group(1)
            elif m := denorm_pattern.search(line):
                current_metrics['Denormalization Time'] = m.group(1)
            elif m := algo_pattern.search(line):
                current_metrics['Algorithm Time'] = m.group(1)
            elif m := voids_pattern.search(line):
                current_metrics['Voids'] = m.group(1)
                
                # "Total Number of Voids" is the last line of a block.
                # Once we hit this, we save the accumulated row to our list.
                parsed_data.append([
                    current_grid,
                    current_threads,
                    current_particles,
                    current_metrics.get('Interpolation Time', ''),
                    current_metrics.get('Normalization Time', ''),
                    current_metrics.get('Mover Time', ''),
                    current_metrics.get('Denormalization Time', ''),
                    current_metrics.get('Algorithm Time', ''),
                    current_metrics.get('Voids', '')
                ])

    # Write the extracted data to a CSV file
    headers = [
        'Grid Dimensions', 'Threads', 'Total Particles', 
        'Interpolation Time', 'Normalization Time', 'Mover Time', 
        'Denormalization Time', 'Algorithm Time', 'Voids'
    ]
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(parsed_data)

    print(f"Successfully extracted {len(parsed_data)} rows to '{output_file}'.")

if __name__ == "__main__":
    # Define your input and output filenames here
    INPUT_FILENAME = 'times_parallel_cluster.txt'
    OUTPUT_FILENAME = 'data_cluster.csv'
    
    parse_benchmark_data(INPUT_FILENAME, OUTPUT_FILENAME)