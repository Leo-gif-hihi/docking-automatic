import os
import csv
from pathlib import Path

def extract_free_energy(log_file):
    """Parses a Vina log file and returns the top free energy (affinity) score."""
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
            
        for i, line in enumerate(lines):
            if "mode |   affinity" in line:
                # Top score is typically 3 lines after the header
                # mode |   affinity | dist from best mode
                #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
                # -----+------------+----------+----------
                #    1        -12.9          0          0
                if i + 3 < len(lines):
                    data_line = lines[i + 3]
                    parts = data_line.split()
                    if parts and parts[0] == "1":
                        return float(parts[1])
    except Exception as e:
        print(f"Error reading {log_file}: {e}")
    return None

def rank_complexes(output_dir):
    """Finds all log files in output_dir, extracts energies, and returns a sorted list."""
    output_path = Path(output_dir)
    log_files = list(output_path.glob("*.log"))
    
    results = []
    for log_file in log_files:
        energy = extract_free_energy(log_file)
        if energy is not None:
            # Complex name can be inferred from log file name (e.g., 1IEP_imatinib_vina.log -> 1IEP_imatinib)
            complex_name = log_file.stem
            if complex_name.endswith("_vina"):
                complex_name = complex_name[:-5]
            
            # Split into protein and ligand (first name before first _, rest is ligand)
            parts = complex_name.split("_", 1)
            if len(parts) == 2:
                protein, ligand = parts
            else:
                protein, ligand = complex_name, "N/A"
                
            results.append((protein, ligand, energy))
            
    # Sort by free energy (lowest/most negative first)
    results.sort(key=lambda x: x[2])
    return results

def print_ranking(results, output_csv=None):
    if not results:
        print("No valid log files or energy scores found.")
        return
        
    print("\n--- Ranking of Complexes by Free Energy ---")
    print(f"{'Protein':<15} | {'Ligand':<25} | {'Affinity (kcal/mol)':<20}")
    print("-" * 66)
    for protein, ligand, energy in results:
        print(f"{protein:<15} | {ligand:<25} | {energy:<20.2f}")

    if output_csv:
        try:
            with open(output_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Protein', 'Ligand', 'Affinity (kcal/mol)'])
                for protein, ligand, energy in results:
                    writer.writerow([protein, ligand, energy])
            print(f"\nRanking saved to {output_csv}")
        except Exception as e:
            print(f"Error saving to CSV: {e}")

