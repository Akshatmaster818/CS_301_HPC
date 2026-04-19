import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def generate_hpc_plots(csv_filename):
    if not os.path.exists(csv_filename):
        print(f"Error: File '{csv_filename}' not found.")
        return

    # Load the data
    df = pd.read_csv(csv_filename)

    # Clean up column names just in case there are trailing spaces
    df.columns = df.columns.str.strip()

    # Ensure correct data types
    df['Threads'] = df['Threads'].astype(int)
    df['Total Particles'] = df['Total Particles'].astype(int)
    df['Algorithm Time'] = df['Algorithm Time'].astype(float)
    df['Interpolation Time'] = df['Interpolation Time'].astype(float)
    df['Mover Time'] = df['Mover Time'].astype(float)

    # Get unique configurations (Grid Dimensions + Total Particles)
    configs = df[['Grid Dimensions', 'Total Particles']].drop_duplicates()

    print(f"Found {len(configs)} configurations. Generating plots...")

    for index, config in configs.iterrows():
        grid = config['Grid Dimensions']
        particles = config['Total Particles']
        
        # Filter data for the current configuration and sort by threads
        config_data = df[(df['Grid Dimensions'] == grid) & 
                         (df['Total Particles'] == particles)].copy()
        config_data = config_data.sort_values(by='Threads')
        
        threads = config_data['Threads'].values
        algo_times = config_data['Algorithm Time'].values
        interp_times = config_data['Interpolation Time'].values
        mover_times = config_data['Mover Time'].values
        
        # Check if we have serial data (1 thread) to calculate Absolute Speedup
        serial_row = config_data[config_data['Threads'] == 1]
        if not serial_row.empty:
            t_serial = serial_row['Algorithm Time'].values[0]
            speedups = t_serial / algo_times
            baseline_threads = 1
        else:
            print(f"Warning: No 1-thread data for {grid} with {particles} particles. Using lowest thread count as baseline.")
            t_serial = algo_times[0]
            speedups = t_serial / algo_times
            baseline_threads = threads[0]

        # Create a figure with 3 subplots (1 row, 3 columns)
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
        
        # Format particles for the title (e.g., 900,000 or 5.0M)
        if particles >= 1000000:
            part_label = f"{particles/1000000:.1f}M"
        else:
            part_label = f"{particles:,}"
            
        fig.suptitle(f'Performance Analysis: Grid {grid} | Particles: {part_label}', 
                     fontsize=16, fontweight='bold', y=1.02)

        # --- Subplot 1: Execution Time vs Cores ---
        ax1.plot(threads, algo_times, marker='o', linestyle='-', color='#1f77b4', linewidth=2, markersize=8)
        ax1.set_title('Execution Time vs Cores', fontsize=13)
        ax1.set_xlabel('Number of Threads', fontsize=12)
        ax1.set_ylabel('Total Algorithm Time (seconds)', fontsize=12)
        ax1.set_xticks(threads)
        ax1.grid(True, linestyle='--', alpha=0.7)

        # --- Subplot 2: Speedup vs Cores ---
        ax2.plot(threads, speedups, marker='s', linestyle='-', color='#2ca02c', linewidth=2, markersize=8, label='Actual Speedup')
        # Ideal speedup line
        ideal_speedups = threads / baseline_threads
        ax2.plot(threads, ideal_speedups, marker='', linestyle='--', color='#d62728', linewidth=1.5, label='Ideal Speedup')
        
        ax2.set_title('Absolute Speedup vs Cores', fontsize=13)
        ax2.set_xlabel('Number of Threads', fontsize=12)
        ax2.set_ylabel('Speedup ($T_1 / T_p$)', fontsize=12)
        ax2.set_xticks(threads)
        ax2.legend()
        ax2.grid(True, linestyle='--', alpha=0.7)

        # --- Subplot 3: Phase Breakdown (Interpolation vs Mover) ---
        # Using a stacked bar chart to show the bottleneck
        width = 0.6 if len(threads) < 6 else 0.8
        
        # Convert threads array to string for categorical x-axis spacing in bar charts
        thread_labels = [str(t) for t in threads]
        x_indices = np.arange(len(thread_labels))
        
        ax3.bar(x_indices, interp_times, width, label='Interpolation', color='#ff7f0e', alpha=0.8)
        ax3.bar(x_indices, mover_times, width, bottom=interp_times, label='Mover', color='#9467bd', alpha=0.8)
        
        ax3.set_title('Phase Breakdown: Interpolation vs Mover', fontsize=13)
        ax3.set_xlabel('Number of Threads', fontsize=12)
        ax3.set_ylabel('Time (seconds)', fontsize=12)
        ax3.set_xticks(x_indices)
        ax3.set_xticklabels(thread_labels)
        ax3.legend()
        ax3.grid(True, axis='y', linestyle='--', alpha=0.7)

        # Finalize and Save
        plt.tight_layout()
        safe_grid = grid.replace(' ', '').replace('x', 'X')
        filename = f"plot_Assign07_{safe_grid}_{particles}_particles.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved: {filename}")
        
        plt.close()

if __name__ == "__main__":
    # Ensure your CSV is named 'results.csv' or change the variable below
    CSV_FILENAME = 'data_lab.csv' 
    generate_hpc_plots(CSV_FILENAME)