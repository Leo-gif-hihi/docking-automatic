import os
import csv
import logging
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
        logging.error(f"Error reading {log_file}: {e}")
    return None

def rank_complexes(output_dir):
    """Finds all log files in output_dir, extracts energies, and returns a sorted list."""
    output_path = Path(output_dir)
    log_files = list(output_path.rglob("*.log"))
    
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
        logging.warning("No valid log files or energy scores found.")
        return
        
    logging.info("\n--- Ranking of Complexes by Free Energy ---")
    logging.info(f"{'Protein':<15} | {'Ligand':<25} | {'Affinity (kcal/mol)':<20}")
    logging.info("-" * 66)
    for protein, ligand, energy in results:
        logging.info(f"{protein:<15} | {ligand:<25} | {energy:<20.2f}")

    if output_csv:
        try:
            with open(output_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Protein', 'Ligand', 'Affinity (kcal/mol)'])
                for protein, ligand, energy in results:
                    writer.writerow([protein, ligand, energy])
            logging.info(f"\nRanking saved to {output_csv}")
            
            # Check for isomers and create a best isomers ranking
            has_isomers = any("_isomer_" in ligand for _, ligand, _ in results)
            if has_isomers:
                best_isomers = {}
                for protein, ligand, energy in results:
                    if "_isomer_" in ligand:
                        base_ligand = ligand.split("_isomer_")[0]
                    else:
                        base_ligand = ligand
                        
                    key = (protein, base_ligand)
                    # Keep the one with the lowest energy
                    if key not in best_isomers or energy < best_isomers[key][1]:
                        best_isomers[key] = (ligand, energy)
                        
                best_results = [(p, data[0], data[1]) for (p, _), data in best_isomers.items()]
                best_results.sort(key=lambda x: x[2])  # Sort by energy
                
                best_csv = Path(output_csv).with_name(f"best_isomers_{Path(output_csv).name}")
                with open(best_csv, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Protein', 'Best_Isomer_Ligand', 'Affinity (kcal/mol)'])
                    for protein, ligand, energy in best_results:
                        writer.writerow([protein, ligand, energy])
                logging.info(f"Best isomers ranking saved to {best_csv}")
                
        except Exception as e:
            logging.error(f"Error saving to CSV: {e}")

