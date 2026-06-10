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

def _get_base_ligand(ligand):
    """Extracts the base ligand name, ignoring isomer suffixes."""
    return ligand.split("_isomer_")[0] if "_isomer_" in ligand else ligand

def _load_cid_to_name(ligand_path=None):
    """Loads a mapping of CID to compound names from PubChem_compound_summary_list.csv."""
    import csv
    from pathlib import Path
    import logging
    
    cid_to_name = {}
    summary_csv_path = Path(ligand_path if ligand_path else "ligand-test") / "PubChem_compound_summary_list.csv"
    if summary_csv_path.exists():
        try:
            with open(summary_csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if "Compound_CID" in row and "Name" in row:
                        cid_to_name[row["Compound_CID"]] = row["Name"]
        except Exception as e:
            logging.warning(f"Failed to read summary CSV: {e}")
    return cid_to_name

from contextlib import contextmanager

@contextmanager
def _prepare_prolif_protein(protein_cif):
    """Context manager to convert CIF to PDB and yield a ProLIF Molecule."""
    import tempfile
    import os
    from pathlib import Path
    import MDAnalysis as mda
    import prolif as plf
    
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_pdb:
        tmp_prot_pdb = Path(tmp_pdb.name)
        
    try:
        _convert_to_pdb(protein_cif, "cif", tmp_prot_pdb)
        u = mda.Universe(str(tmp_prot_pdb))
        protein_mol = plf.Molecule.from_mda(u)
        yield protein_mol
    finally:
        try:
            os.unlink(tmp_prot_pdb)
        except OSError:
            pass

def _generate_network_plot(fp, pose, driver, out_png):
    """Generates a 2D ProLIF interaction network and saves a screenshot using Selenium."""
    import tempfile
    import os
    import time
    import logging
    
    view = fp.plot_lignetwork(pose, kind="frame", frame=0)
    
    with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as tmp_html:
        tmp_html_path = tmp_html.name
        
    try:
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
        
        driver.get(f"file://{tmp_html_path}")
        time.sleep(2)
        driver.save_screenshot(str(out_png))
        logging.debug(f"Saved ProLif visualization: {out_png}")
    finally:
        try:
            os.unlink(tmp_html_path)
        except OSError:
            pass

def _generate_barcode_plot(df_combined, pdf_path, tiff_path, protein_pocket_base):
    """Generates a Barcode plot from a combined DataFrame."""
    import matplotlib.pyplot as plt
    from prolif.plotting.barcode import Barcode
    from logger_utils import log_step
    
    sorted_columns = df_combined.sum().sort_values(ascending=False).index
    df_sorted = df_combined[sorted_columns]

    ax = Barcode(df_sorted).display(
        figsize=(10, min(8, max(4, len(sorted_columns) * 0.3))),
        dpi=300,
        xlabel="Compound"
    )
    
    fig = ax.figure
    fig.savefig(str(pdf_path), format="pdf", bbox_inches="tight")
    fig.savefig(str(tiff_path), format="tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)
    
    log_step(None, f"Saved barcode for {protein_pocket_base}", color="white")

def print_ranking(results, output_csv=None, ligand_path=None):
    if not results:
        logging.warning("No valid log files or energy scores found.")
        return []
        
    best_dict = {}
    for protein, pocket, ligand, run, energy in results:
        base_ligand = _get_base_ligand(ligand)
            
        key = (protein, pocket, base_ligand)
        # Keep the one with the lowest energy
        if key not in best_dict or energy < best_dict[key][4]:
            best_dict[key] = (protein, pocket, ligand, run, energy)
            
    curated_results = list(best_dict.values())
    curated_results.sort(key=lambda x: x[4])  # Sort by energy

    # Load CID to Name mapping
    cid_to_name = _load_cid_to_name(ligand_path)

    import statistics

    energies_dict = {}
    for protein, pocket, ligand, run, energy in results:
        key = (protein, pocket, ligand)
        if key not in energies_dict:
            energies_dict[key] = []
        energies_dict[key].append(energy)

    def format_row(row, is_curated=False):
        protein, pocket, ligand, run, energy = row
        cid = _get_base_ligand(ligand)
        if "_isomer_" in ligand:
            parts = ligand.split("_isomer_")
            isomer = parts[1] if len(parts) > 1 else ""
        else:
            isomer = "N/A"
        ligand_name = cid_to_name.get(cid, "N/A")
        
        base_res = [protein, pocket, cid, ligand_name, isomer, run, energy]
        if is_curated:
            key = (protein, pocket, ligand)
            energies = energies_dict.get(key, [energy])
            mean_e = statistics.mean(energies)
            if len(energies) > 1:
                sd_e = statistics.stdev(energies)
                mean_sd_str = f"{mean_e:.2f} +- {sd_e:.2f}"
            else:
                mean_sd_str = f"{mean_e:.2f} +- N/A"
            base_res.append(mean_sd_str)
        return base_res

    header_raw = ['Protein', 'Pocket ID', 'CID', 'Ligand_name', 'Isomer', 'Run', 'Affinity (kcal/mol)']
    header_curated = header_raw + ['Mean+-SD (kcal/mol)']

    formatted_results = [format_row(row, is_curated=False) for row in results]
    formatted_curated = [format_row(row, is_curated=True) for row in curated_results]

    print()
    display_limit = 10
    log_step(None, f"--- Top {min(display_limit, len(formatted_curated))} Curated Best Complexes by Free Energy (Total: {len(formatted_curated)}) ---", color="magenta")
    log_step(None, f"{'Protein':<15} | {'Pocket ID':<10} | {'CID':<15} | {'Affinity (kcal/mol)':<20} | {'Mean+-SD (kcal/mol)':<20}", color="magenta")
    log_step(None, "-" * 92, color="magenta")
    for row in formatted_curated[:display_limit]:
        protein, pocket, cid, lname, iso, run, energy, mean_sd = row
        log_step(None, f"{protein:<15} | {pocket:<10} | {cid:<15} | {energy:<20.2f} | {mean_sd:<20}", color="magenta")

    if output_csv:
        try:
            with open(output_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header_raw)
                writer.writerows(formatted_results)
            print()
            log_step(None, f"Raw ranking saved to {output_csv}", color="magenta")
            
            best_csv = Path(output_csv).with_name(f"curated_best_{Path(output_csv).name}")
            with open(best_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header_curated)
                writer.writerows(formatted_curated)
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

def visualize_prolif_results(results, output_dir, protein_clean_dir, ligand_path=None, display_limit=20):
    if not results:
        return
        
    import os
    import logging
    import pandas as pd
    import matplotlib.pyplot as plt
    from pathlib import Path
    from logger_utils import log_step
    import warnings
    import sys
    
    warnings.simplefilter("ignore")
    os.environ["PYTHONWARNINGS"] = "ignore"
    
    devnull = open(os.devnull, 'w')
    old_stderr = sys.stderr
    sys.stderr = devnull
    try:
        import MDAnalysis as mda
        import prolif as plf
        from prolif.plotting.barcode import Barcode
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

    # Map CIDs to standard names
    cid_to_name = _load_cid_to_name(ligand_path)

    # Set publication typography globally
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 10,
        'legend.fontsize': 10
    })
    Barcode.COLORS["Hydrophobic"] = "#4477AA"
    Barcode.COLORS["HBAcceptor"] = "#EE6677"
    Barcode.COLORS["HBDonor"] = "#228833"
    Barcode.COLORS["PiStacking"] = "#CCBB44"

    # Group results by Protein and Pocket
    targets = {}
    for row in results[:display_limit]:
        protein, pocket, ligand, run, energy = row
        key = (protein, pocket)
        if key not in targets:
            targets[key] = []
        targets[key].append(row)

    try:
        all_interactions = plf.Fingerprint.list_available()
        for (protein, pocket), target_results in targets.items():
            protein_pocket_base = f"{protein}_pocket_{pocket}" if pocket != "N/A" else protein
            protein_cif = Path(protein_clean_dir) / f"{protein}FH.cif"
            
            if not protein_cif.exists():
                logging.warning(f"Missing {protein_cif}. Skipping visualization for {protein_pocket_base}.")
                continue

            fp_dataframes = []
            ligand_labels = []

            try:
                with _prepare_prolif_protein(protein_cif) as protein_mol:
                    for row in target_results:
                        _, _, ligand, run, _ = row
                        ligand_sdf = Path(output_dir) / "vina_output" / f"run_{run}" / f"{protein_pocket_base}_{ligand}" / f"{protein_pocket_base}_{ligand}_vina_out.sdf"
                        
                        if not ligand_sdf.exists():
                            continue

                        base_ligand = _get_base_ligand(ligand)
                        label = cid_to_name.get(base_ligand, base_ligand)
                        
                        poses = list(plf.sdf_supplier(str(ligand_sdf)))
                        if not poses:
                            continue
                            
                        fp = plf.Fingerprint(all_interactions)
                        fp.run_from_iterable([poses[0]], protein_mol, progress=False)
                        
                        # Generate 2D network screenshot
                        out_png = vis_dir / f"{protein_pocket_base}_{ligand}_prolif.png"
                        try:
                            _generate_network_plot(fp, poses[0], driver, out_png)
                        except Exception as e:
                            logging.error(f"Failed to generate network plot for {protein_pocket_base}_{ligand}: {e}")
                        
                        df_i = fp.to_dataframe()
                        
                        if 'ligand' in df_i.columns.names:
                            idx = df_i.columns.names.index('ligand')
                            new_tuples = [tuple('LIG' if i == idx else val for i, val in enumerate(tup)) for tup in df_i.columns]
                            df_i.columns = pd.MultiIndex.from_tuples(new_tuples, names=df_i.columns.names)
                            
                        fp_dataframes.append(df_i)
                        ligand_labels.append(label)
                        
            except Exception as e:
                logging.error(f"Failed to process protein {protein} or its ligands for ProLif visualization: {e}")
                continue

            if fp_dataframes:
                df_combined = pd.concat(fp_dataframes, ignore_index=True).fillna(False)
                df_combined.index = ligand_labels
                
                pdf_path = vis_dir / f"{protein_pocket_base}_barcode.pdf"
                tiff_path = vis_dir / f"{protein_pocket_base}_barcode.tiff"
                try:
                    _generate_barcode_plot(df_combined, pdf_path, tiff_path, protein_pocket_base)
                except Exception as e:
                    logging.error(f"Failed to generate barcode plot for {protein_pocket_base}: {e}")

        log_step(None, f"All ProLif visualizations saved to {vis_dir}", color="magenta")
    finally:
        driver.quit()