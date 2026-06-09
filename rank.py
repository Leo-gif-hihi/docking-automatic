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
                
            run_index = "1"
            for part in log_file.parts:
                run_match = re.match(r'^run_(\d+)$', part)
                if run_match:
                    run_index = run_match.group(1)
                    break
            
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
                    results.append((protein, pocket, ligand, run_index, energy))
                    continue
                else:
                    # Fallback: split by first underscore
                    parts = complex_name.split("_", 1)
                    if len(parts) == 2:
                        protein, ligand = parts
                    else:
                        protein, ligand = complex_name, "N/A"
                    pocket = "N/A"
                    results.append((protein, pocket, ligand, run_index, energy))
                    continue

            # If we successfully extracted protein_pocket and ligand using known_ligands
            match = re.match(r'^(.*?)_pocket_(\d+)$', protein_pocket)
            if match:
                protein = match.group(1)
                pocket = match.group(2)
            else:
                protein = protein_pocket
                pocket = "N/A"
                
            results.append((protein, pocket, ligand, run_index, energy))
            
    # Sort by free energy (lowest/most negative first)
    results.sort(key=lambda x: x[4])
    return results

def print_ranking(results, output_csv=None):
    if not results:
        logging.warning("No valid log files or energy scores found.")
        return []
        
    best_dict = {}
    for protein, pocket, ligand, run, energy in results:
        if "_isomer_" in ligand:
            base_ligand = ligand.split("_isomer_")[0]
        else:
            base_ligand = ligand
            
        key = (protein, pocket, base_ligand)
        # Keep the one with the lowest energy
        if key not in best_dict or energy < best_dict[key][4]:
            best_dict[key] = (protein, pocket, ligand, run, energy)
            
    curated_results = list(best_dict.values())
    curated_results.sort(key=lambda x: x[4])  # Sort by energy

    print()
    display_limit = 10
    log_step(None, f"--- Top {min(display_limit, len(curated_results))} Curated Best Complexes by Free Energy (Total: {len(curated_results)}) ---", color="magenta")
    log_step(None, f"{'Protein':<15} | {'Pocket ID':<10} | {'Ligand':<25} | {'Run':<5} | {'Affinity (kcal/mol)':<20}", color="magenta")
    log_step(None, "-" * 87, color="magenta")
    for protein, pocket, ligand, run, energy in curated_results[:display_limit]:
        log_step(None, f"{protein:<15} | {pocket:<10} | {ligand:<25} | {run:<5} | {energy:<20.2f}", color="magenta")

    if output_csv:
        try:
            with open(output_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Protein', 'Pocket ID', 'Ligand', 'Run', 'Affinity (kcal/mol)'])
                for protein, pocket, ligand, run, energy in results:
                    writer.writerow([protein, pocket, ligand, run, energy])
            print()
            log_step(None, f"Raw ranking saved to {output_csv}", color="magenta")
            
            best_csv = Path(output_csv).with_name(f"curated_best_{Path(output_csv).name}")
            with open(best_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Protein', 'Pocket ID', 'Best_Ligand', 'Run', 'Affinity (kcal/mol)'])
                for protein, pocket, ligand, run, energy in curated_results:
                    writer.writerow([protein, pocket, ligand, run, energy])
            log_step(None, f"Curated best ranking saved to {best_csv}", color="magenta")
                
        except Exception as e:
            logging.error(f"Error saving to CSV: {e}")

    return curated_results

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

def _map_protein_chains(tmp_prot_pdb, mapped_prot_pdb, protein, asym_id_mapping):
    import logging
    import string
    
    protein_mapping = asym_id_mapping.get(protein) or asym_id_mapping.get(f"{protein}.cif") or {}
    prody_to_orig = {v: k for k, v in protein_mapping.items()}
    
    existing_chains = set()
    with open(tmp_prot_pdb, 'r') as in_f:
        for line in in_f:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                existing_chains.add(line[21])
                
    used_chains = set()
    for current_chain, orig_chain in prody_to_orig.items():
        if len(orig_chain) == 1:
            used_chains.add(orig_chain)
            
    for current_chain in existing_chains:
        if current_chain not in prody_to_orig:
            used_chains.add(current_chain)
            
    all_candidates = "ZYXWVUTSRQPONMLKJIHGFEDCBA" + string.ascii_lowercase[::-1] + string.digits[::-1]
    
    final_mapping = {}
    for current_chain, orig_chain in prody_to_orig.items():
        if len(orig_chain) == 1:
            final_mapping[current_chain] = orig_chain
        else:
            assigned = False
            for candidate in all_candidates:
                if candidate not in used_chains:
                    final_mapping[current_chain] = candidate
                    used_chains.add(candidate)
                    logging.warning(f"Original chain ID '{orig_chain}' for {protein} is > 1 char. Reassigned to '{candidate}'.")
                    assigned = True
                    break
            if not assigned:
                logging.warning(f"Cannot reassign chain for '{orig_chain}' in {protein}. Legacy PDB limit reached.")
                return False, None
    
    with open(tmp_prot_pdb, 'r') as in_f, open(mapped_prot_pdb, 'w') as out_f:
        for line in in_f:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                current_chain = line[21]
                if current_chain in final_mapping:
                    line = line[:21] + final_mapping[current_chain] + line[22:]
            out_f.write(line)
            
    orig_to_final = {orig: final_mapping[curr] for curr, orig in prody_to_orig.items() if curr in final_mapping}
    return True, orig_to_final

def generate_complexes(results, output_dir, protein_clean_dir, display_limit=20):
    if not results:
        return
    
    import subprocess
    import tempfile
    import logging
    import os
    import json
    from pathlib import Path
    from logger_utils import log_step

    vis_dir = Path(output_dir) / "visualization"
    os.makedirs(vis_dir, exist_ok=True)
    
    mapping_file = Path(protein_clean_dir) / "asym_id_mapping.json"
    asym_id_mapping = {}
    if mapping_file.exists():
        try:
            with open(mapping_file, "r") as f:
                asym_id_mapping = json.load(f)
        except json.JSONDecodeError:
            logging.warning(f"Failed to parse {mapping_file}")

    log_step(None, f"Generating PDB complex files for top {min(display_limit, len(results))} results...", color="white")
    
    complex_chain_mappings = {}
    for protein, pocket, ligand, run, energy in results[:display_limit]:
        protein_pocket_base = f"{protein}_pocket_{pocket}" if pocket != "N/A" else protein
        
        # Paths
        protein_cif = Path(protein_clean_dir) / f"{protein}FH.cif"
        
        ligand_sdf = Path(output_dir) / "vina_output" / f"run_{run}" / f"{protein_pocket_base}_{ligand}" / f"{protein_pocket_base}_{ligand}_vina_out.sdf"
        
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
                
            # Apply original chain ID mapping
            mapped_prot_pdb = Path(tmpdir) / "mapped_prot.pdb"
            success, orig_to_final = _map_protein_chains(tmp_prot_pdb, mapped_prot_pdb, protein, asym_id_mapping)
            if not success:
                logging.warning(f"Skipping complex for {protein} due to legacy PDB chain limitation.")
                continue

            # Determine available chain ID for ligand
            lig_chain = _find_available_chain(mapped_prot_pdb)
            if lig_chain is None:
                logging.warning(f"No available chain IDs left for {protein_pocket_base}_{ligand} (Legacy PDB limit reached). Skipping complex generation.")
                continue

            # Combine PDBs
            try:
                _combine_complex_pdbs(mapped_prot_pdb, tmp_lig_pdb, complex_pdb, lig_chain)
                logging.debug(f"Created complex: {complex_pdb}")
                orig_to_final["LIGAND"] = lig_chain
                complex_chain_mappings[f"{protein_pocket_base}_{ligand}"] = orig_to_final
            except Exception as e:
                logging.warning(f"Skipping complex {protein_pocket_base}_{ligand}: {e}")
                if complex_pdb.exists():
                    try:
                        complex_pdb.unlink()
                    except OSError:
                        pass
                
    if complex_chain_mappings:
        mapping_out_file = vis_dir / "final_chain_mapping.json"
        try:
            with open(mapping_out_file, "w") as f:
                json.dump(complex_chain_mappings, f, indent=4)
            log_step(None, f"Final chain mappings saved to {mapping_out_file}", color="white")
        except Exception as e:
            logging.error(f"Failed to save final chain mappings: {e}")
            
    log_step(None, f"Complexes saved to {vis_dir}", color="white")

def visualize_prolif_complexes(results, output_dir, protein_clean_dir, display_limit=20):
    if not results:
        return
        
    import os
    import logging
    from pathlib import Path
    import tempfile
    import time
    from logger_utils import log_step
    import warnings
    import sys
    
    # Aggressively ignore all warnings during the visualization
    warnings.simplefilter("ignore")
    os.environ["PYTHONWARNINGS"] = "ignore"
    
    # Ultimate silence: redirect stderr to /dev/null during imports
    devnull = open(os.devnull, 'w')
    old_stderr = sys.stderr
    sys.stderr = devnull
    try:
        import MDAnalysis as mda
        import prolif as plf
        from selenium import webdriver
        from selenium.webdriver.chrome.options import Options
    except ImportError as e:
        sys.stderr = old_stderr
        devnull.close()
        logging.error(f"Missing dependency for ProLif visualization: {e}")
        log_step(None, f"Skipping ProLif visualization due to missing dependency: {e}", color="yellow")
        return
    finally:
        sys.stderr = old_stderr
        devnull.close()

    vis_dir = Path(output_dir) / "visualization" / "prolif"
    os.makedirs(vis_dir, exist_ok=True)
    
    log_step(None, f"Generating ProLif visualizations for top {min(display_limit, len(results))} results...", color="white")

    # Setup headless chrome
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--window-size=1200,800")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    
    try:
        driver = webdriver.Chrome(options=chrome_options)
    except Exception as e:
        logging.error(f"Failed to initialize Selenium Chrome driver: {e}")
        log_step(None, "Skipping ProLif visualization because Selenium Chrome driver could not be initialized.", color="yellow")
        return

    try:
        for protein, pocket, ligand, run, energy in results[:display_limit]:
            protein_pocket_base = f"{protein}_pocket_{pocket}" if pocket != "N/A" else protein
            
            # Load the FH.cif file
            protein_cif = Path(protein_clean_dir) / f"{protein}FH.cif"
            ligand_sdf = Path(output_dir) / "vina_output" / f"run_{run}" / f"{protein_pocket_base}_{ligand}" / f"{protein_pocket_base}_{ligand}_vina_out.sdf"
            
            out_png = vis_dir / f"{protein_pocket_base}_{ligand}_prolif.png"
            
            if not protein_cif.exists():
                logging.warning(f"Protein file missing: {protein_cif}. Skipping ProLif.")
                continue
                
            if not ligand_sdf.exists():
                logging.warning(f"Ligand SDF file missing: {ligand_sdf}. Skipping ProLif.")
                continue

            try:
                # 1. Protein Preparation (Convert CIF to PDB for MDAnalysis)
                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_pdb:
                    tmp_prot_pdb = Path(tmp_pdb.name)
                
                try:
                    # _convert_to_pdb is defined earlier in rank.py
                    _convert_to_pdb(protein_cif, "cif", tmp_prot_pdb)
                    u = mda.Universe(str(tmp_prot_pdb))
                    protein_mol = plf.Molecule.from_mda(u)
                finally:
                    try:
                        os.unlink(tmp_prot_pdb)
                    except OSError:
                        pass

                # 2. Docking Poses Preparation
                pose_iterable = plf.sdf_supplier(str(ligand_sdf))

                # 3. Fingerprint Generation
                fp = plf.Fingerprint()
                fp.run_from_iterable(pose_iterable, protein_mol, progress=False)

                # 4. Visualization (Target Pose 0)
                pose_index = 0
                view = fp.plot_lignetwork(pose_iterable[pose_index], kind="frame", frame=pose_index)
                
                with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as tmp_html:
                    tmp_html_path = tmp_html.name
                    
                # Save Pyvis/HTML output
                if hasattr(view, 'save'):
                    view.save(tmp_html_path)
                elif hasattr(view, 'write_html'):
                    view.write_html(tmp_html_path)
                elif hasattr(view, 'data'):
                    with open(tmp_html_path, 'w', encoding='utf-8') as f:
                        f.write(view.data)
                else:
                    with open(tmp_html_path, 'w', encoding='utf-8') as f:
                        f.write(str(view))
                
                # Use Selenium to take a screenshot
                driver.get(f"file://{tmp_html_path}")
                time.sleep(2)  # Wait for Pyvis network to stabilize and render
                driver.save_screenshot(str(out_png))
                
                # Cleanup temp HTML
                try:
                    os.unlink(tmp_html_path)
                except OSError:
                    pass
                    
                logging.debug(f"Saved ProLif visualization: {out_png}")
            except Exception as e:
                logging.error(f"Failed to generate ProLif visualization for {protein_pocket_base}_{ligand}: {e}")
                
        log_step(None, f"ProLif visualizations saved to {vis_dir}", color="white")
    finally:
        driver.quit()
