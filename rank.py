import os
import csv
import logging
import re
from pathlib import Path
from logger_utils import log_step

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

def rank_complexes(output_dir, known_ligands=None):
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
            
            protein_pocket = complex_name
            ligand = "N/A"

            # If we know the ligand names, use them to split accurately from the end
            if known_ligands:
                # Sort ligands by length descending to match longest possible ligand name first
                for lig in sorted(known_ligands, key=len, reverse=True):
                    if complex_name.endswith("_" + lig):
                        ligand = lig
                        protein_pocket = complex_name[:-len(lig)-1]
                        break
            
            # Fallback if no known_ligands or no match
            if ligand == "N/A":
                # Try to extract pocket ID first to see where the protein ends
                match = re.match(r'^(.*?)_pocket_(\d+)_(.*)$', complex_name)
                if match:
                    protein = match.group(1)
                    pocket = match.group(2)
                    ligand = match.group(3)
                    results.append((protein, pocket, ligand, energy))
                    continue
                else:
                    # Fallback: split by first underscore
                    parts = complex_name.split("_", 1)
                    if len(parts) == 2:
                        protein, ligand = parts
                    else:
                        protein, ligand = complex_name, "N/A"
                    pocket = "N/A"
                    results.append((protein, pocket, ligand, energy))
                    continue

            # If we successfully extracted protein_pocket and ligand using known_ligands
            match = re.match(r'^(.*?)_pocket_(\d+)$', protein_pocket)
            if match:
                protein = match.group(1)
                pocket = match.group(2)
            else:
                protein = protein_pocket
                pocket = "N/A"
                
            results.append((protein, pocket, ligand, energy))
            
    # Sort by free energy (lowest/most negative first)
    results.sort(key=lambda x: x[3])
    return results

def print_ranking(results, output_csv=None):
    if not results:
        logging.warning("No valid log files or energy scores found.")
        return
        
    print()
    display_limit = 20
    log_step(None, f"--- Top {min(display_limit, len(results))} Ranking of Complexes by Free Energy (Total: {len(results)}) ---", color="magenta")
    log_step(None, f"{'Protein':<15} | {'Pocket ID':<10} | {'Ligand':<25} | {'Affinity (kcal/mol)':<20}", color="magenta")
    log_step(None, "-" * 79, color="magenta")
    for protein, pocket, ligand, energy in results[:display_limit]:
        log_step(None, f"{protein:<15} | {pocket:<10} | {ligand:<25} | {energy:<20.2f}", color="magenta")

    if output_csv:
        try:
            with open(output_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Protein', 'Pocket ID', 'Ligand', 'Affinity (kcal/mol)'])
                for protein, pocket, ligand, energy in results:
                    writer.writerow([protein, pocket, ligand, energy])
            print()
            log_step(None, f"Ranking saved to {output_csv}", color="magenta")
            
            # Check for isomers and create a best isomers ranking
            has_isomers = any("_isomer_" in ligand for _, _, ligand, _ in results)
            if has_isomers:
                best_isomers = {}
                for protein, pocket, ligand, energy in results:
                    if "_isomer_" in ligand:
                        base_ligand = ligand.split("_isomer_")[0]
                    else:
                        base_ligand = ligand
                        
                    key = (protein, pocket, base_ligand)
                    # Keep the one with the lowest energy
                    if key not in best_isomers or energy < best_isomers[key][2]:
                        best_isomers[key] = (ligand, energy)
                        
                best_results = [(p, pock, data[0], data[1]) for (p, pock, _), data in best_isomers.items()]
                best_results.sort(key=lambda x: x[3])  # Sort by energy
                
                best_csv = Path(output_csv).with_name(f"best_isomers_{Path(output_csv).name}")
                with open(best_csv, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Protein', 'Pocket ID', 'Best_Isomer_Ligand', 'Affinity (kcal/mol)'])
                    for protein, pocket, ligand, energy in best_results:
                        writer.writerow([protein, pocket, ligand, energy])
                log_step(None, f"Best isomers ranking saved to {best_csv}", color="magenta")
                
        except Exception as e:
            logging.error(f"Error saving to CSV: {e}")

def _convert_to_pdb(input_file, in_format, output_file, extra_args=None):
    import subprocess
    cmd = ["obabel", "-i", in_format, str(input_file)]
    if extra_args:
        cmd.extend(extra_args)
    cmd.extend(["-o", "pdb", "-O", str(output_file)])
    subprocess.run(cmd, check=True, capture_output=True)

def _find_available_chain(pdb_file):
    import string
    used_chains = set()
    if pdb_file.exists():
        with open(pdb_file, 'r') as p_f:
            for line in p_f:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                    used_chains.add(line[21])
    
    # Legacy PDB allows A-Z, a-z, and 0-9
    all_candidates = "ZYXWVUTSRQPONMLKJIHGFEDCBA" + string.ascii_lowercase[::-1] + string.digits[::-1]
    for candidate in all_candidates:
        if candidate not in used_chains:
            return candidate
    return None

def _combine_complex_pdbs(prot_pdb, lig_pdb, complex_pdb, lig_chain):
    max_atom_serial = 0
    with open(complex_pdb, 'w') as out_f:
        # Write protein, skipping END and CONECT lines
        with open(prot_pdb, 'r') as p_f:
            for line in p_f:
                if not line.startswith(("END", "CONECT", "MASTER")):
                    if line.startswith(("ATOM", "HETATM")) and len(line) >= 11:
                        try:
                            serial = int(line[6:11].strip())
                            if serial > max_atom_serial:
                                max_atom_serial = serial
                        except ValueError:
                            pass
                    out_f.write(line)
        # Write ligand, skipping CONECT lines and adding a chain ID
        with open(lig_pdb, 'r') as l_f:
            for line in l_f:
                if not line.startswith(("END", "CONECT", "MASTER", "COMPND", "AUTHOR")):
                    # Assign the unique chain to the ligand to differentiate it from the protein
                    if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                        line = line[:21] + lig_chain + line[22:]
                        
                        # Renumber atom serial
                        max_atom_serial += 1
                        serial_str = f"{max_atom_serial:>5}"
                        # If it exceeds 5 characters (e.g. 100,000), typical PDB format overflows.
                        if len(serial_str) > 5:
                            raise ValueError(f"Atom serial number {max_atom_serial} exceeds 99999 (legacy PDB limit).")
                        line = line[:6] + serial_str + line[11:]
                        
                    out_f.write(line)
        # Write END
        out_f.write("END\n")

def generate_complexes(results, output_dir, protein_clean_dir, display_limit=20):
    if not results:
        return
    
    import subprocess
    import tempfile
    import logging
    import os
    from pathlib import Path
    from logger_utils import log_step

    vis_dir = Path(output_dir) / "visualization"
    os.makedirs(vis_dir, exist_ok=True)
    
    log_step(None, f"Generating PDB complex files for top {min(display_limit, len(results))} results...", color="white")
    
    for protein, pocket, ligand, energy in results[:display_limit]:
        protein_pocket_base = f"{protein}_pocket_{pocket}" if pocket != "N/A" else protein
        
        # Paths
        protein_cif = Path(protein_clean_dir) / f"{protein}FH.cif"
        
        ligand_sdf = Path(output_dir) / "vina_output" / f"{protein_pocket_base}_{ligand}" / f"{protein_pocket_base}_{ligand}_vina_out.sdf"
        
        complex_pdb = vis_dir / f"{protein_pocket_base}_{ligand}_complex.pdb"
        
        if not protein_cif.exists():
            logging.warning(f"Protein file missing: {protein_cif}. Skipping complex generation.")
            continue
            
        if not ligand_sdf.exists():
            logging.warning(f"Ligand SDF file missing: {ligand_sdf}. Skipping complex generation.")
            continue
            
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_prot_pdb = Path(tmpdir) / "prot.pdb"
            tmp_lig_pdb = Path(tmpdir) / "lig.pdb"
            
            # Convert CIF to PDB
            try:
                _convert_to_pdb(protein_cif, "cif", tmp_prot_pdb)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to convert protein CIF to PDB: {e}")
                continue
                
            # Convert SDF to PDB
            try:
                _convert_to_pdb(ligand_sdf, "sdf", tmp_lig_pdb, extra_args=["-f", "1", "-l", "1"])
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to convert ligand SDF to PDB: {e}")
                continue
                
            # Determine available chain ID for ligand
            lig_chain = _find_available_chain(tmp_prot_pdb)
            if lig_chain is None:
                logging.warning(f"No available chain IDs left for {protein_pocket_base}_{ligand} (Legacy PDB limit reached). Skipping complex generation.")
                continue

            # Combine PDBs
            try:
                _combine_complex_pdbs(tmp_prot_pdb, tmp_lig_pdb, complex_pdb, lig_chain)
                logging.debug(f"Created complex: {complex_pdb}")
            except Exception as e:
                logging.warning(f"Skipping complex {protein_pocket_base}_{ligand}: {e}")
                if complex_pdb.exists():
                    try:
                        complex_pdb.unlink()
                    except OSError:
                        pass
                
    log_step(None, f"Complexes saved to {vis_dir}", color="cyan")
