import os
import glob
import argparse
import sys
import logging
from prody import parsePDB, writePDB, confProDy

# Turn off ProDy's progress text to keep your terminal clean
confProDy(verbosity='none')

def get_hetatms(structure):
    """Extracts unique HETATM residue names from a structure, excluding water."""
    hetatms = structure.select('not protein and not water')
    if hetatms is not None:
        return set(hetatms.getResnames())
    return set()

def get_selection_string(to_eliminate, all_hetatms):
    """Generates the ProDy selection string based on residues to eliminate."""
    if set(to_eliminate) == all_hetatms and all_hetatms:
        return "protein"
    elif to_eliminate:
        exclude_str = " ".join(to_eliminate)
        return f"protein or (not water and not resname {exclude_str})"
    return "protein or (not water)"

def prompt_elimination(hetatms_set):
    """Prompts the user for which HETATMs to eliminate."""
    logging.info(f"Found the following HETATM residues: {', '.join(hetatms_set)}")
    logging.info("Which ones would you like to ELIMINATE?")
    try:
        user_input = input("Enter them separated by commas (type 'all' to delete all HETATM, press Enter to delete only water): ").strip()
    except KeyboardInterrupt:
        logging.error("\n\nProcess interrupted by user (Ctrl+C). Exiting...")
        sys.exit(0)
        
    if user_input == "":
        return []
    if user_input.lower() == "all":
        return list(hetatms_set)
    return [res.strip() for res in user_input.split(',')]

def clean_and_save_pdb(structure, original_filepath, output_dir, sel_str):
    """Applies the selection string to clean the structure and saves it to the output directory."""
    filename = os.path.basename(original_filepath)
    clean_selection = structure.select(sel_str)
    out_filepath = os.path.join(output_dir, filename)
    
    if clean_selection:
            writePDB(out_filepath, clean_selection)
            
            # --- Restore CRYST1 line for reduce2 ---
            cryst1_line = None
            # Read the original PDB file to find the CRYST1 line
            with open(original_filepath, 'r') as f:
                for line in f:
                    if line.startswith("CRYST1"):
                        cryst1_line = line
                        break
            
            # If there is CRYST1, 
            # insert it on the first line of the cleaned file
            if cryst1_line:
                with open(out_filepath, 'r') as f:
                    content = f.read()
                with open(out_filepath, 'w') as f:
                    f.write(cryst1_line)
                    f.write(content)

            logging.debug(f" -> Saved: {out_filepath}")
            return out_filepath
    else:
        logging.debug(f" -> Skipped {filename}: No atoms left after cleaning.")
        return None

def clean_global_mode(pdb_files, output_dir):
    """Runs the cleaning process in global mode (one prompt for all files)."""
    logging.info(f"Scanning {len(pdb_files)} files for HETATM residues...")
    all_hetatms = set()
    structures = []
    
    # Pre-parse and collect all HETATMs
    for filepath in pdb_files:
        structure = parsePDB(filepath)
        if structure is None:
            continue
        structures.append((filepath, structure))
        all_hetatms.update(get_hetatms(structure))
    
    if not all_hetatms:
        logging.debug("No HETATM residues found in any of the PDB files.")
        to_eliminate = []
    else:
        to_eliminate = prompt_elimination(all_hetatms)
        action_desc = "ALL HETATM residues" if set(to_eliminate) == all_hetatms else ", ".join(to_eliminate)
        logging.debug(f"Action: Eliminating {action_desc}.")
        
    sel_str = get_selection_string(to_eliminate, all_hetatms)
    
    cleaned_paths = []
    for filepath, structure in structures:
        cleaned_pdb_path = clean_and_save_pdb(structure, filepath, output_dir, sel_str)
        if cleaned_pdb_path:
            cleaned_paths.append(cleaned_pdb_path)
            
    return cleaned_paths

def clean_local_mode(pdb_files, output_dir):
    """Runs the cleaning process in local mode (prompt per file)."""
    logging.info(f"Running in LOCAL mode. Scanning {len(pdb_files)} files individually.")
    logging.info("Tip: You can press Ctrl+C at any prompt to escape and stop processing.\n")
    
    cleaned_paths = []
    for filepath in pdb_files:
        filename = os.path.basename(filepath)
        structure = parsePDB(filepath)
        if structure is None:
            continue
            
        file_hetatms = get_hetatms(structure)
        logging.debug(f"--- {filename} ---")
        
        if not file_hetatms:
            logging.debug("No HETATM residues found. Just cleaning water.")
            to_eliminate = []
        else:
            to_eliminate = prompt_elimination(file_hetatms)
            action_desc = "ALL HETATM residues in this file" if set(to_eliminate) == file_hetatms else ", ".join(to_eliminate)
            logging.debug(f"Action: Eliminating {action_desc}.")
            
        sel_str = get_selection_string(to_eliminate, file_hetatms)
        cleaned_pdb_path = clean_and_save_pdb(structure, filepath, output_dir, sel_str)
        if cleaned_pdb_path:
            cleaned_paths.append(cleaned_pdb_path)
        logging.debug("")  # Blank line for readability between files if in debug mode
        
    return cleaned_paths

def run_reduce2(protein_path, protein_protonated):
    """Runs mmtbx.reduce2 to add hydrogens and optimize the structure."""
    import subprocess
    from pathlib import Path

    if Path(protein_protonated).exists():
        logging.debug(f"Skipping REDUCE2: {protein_protonated} already exists.")
        return

    cmd_reduce2 = [
        "mmtbx.reduce2", str(protein_path), "approach=add", 
        "add_flip_movers=True", f"output.filename={protein_protonated}"
    ]
    subprocess.run(cmd_reduce2, check=True, capture_output=True, text=True)


def run_meeko_receptor(protein_protonated, protein_prep_out):
    """Runs mk_prepare_receptor.py (Meeko) to prepare the receptor."""
    import subprocess
    from pathlib import Path

    if Path(protein_prep_out).exists():
        logging.debug(f"Skipping Meeko: {protein_prep_out} already exists.")
        return

    cmd_meeko = [
        "mk_prepare_receptor.py", "-i", str(protein_protonated), 
        "-o", str(protein_prep_out), "-p"
    ]
    subprocess.run(cmd_meeko, check=True, capture_output=True, text=True)

def prepare_proteins(input_dir, output_dir, mode):
    """Main workflow to orchestrate the cleaning process."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.debug(f"Created output directory: {output_dir}")
    
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))
    if not pdb_files:
        logging.error(f"No PDB files found in '{input_dir}'. Please check your folder.")
        return {}

    cleaned_paths = []
    if mode == "global":
        cleaned_paths = clean_global_mode(pdb_files, output_dir)
    elif mode == "local":
        cleaned_paths = clean_local_mode(pdb_files, output_dir)
        
    from pathlib import Path
    prepared_results = {}
    
    for cleaned_pdb_path in cleaned_paths:
        protein_path = Path(cleaned_pdb_path)
        protein_dir = protein_path.parent
        protein_base = protein_path.stem

        protein_protonated = protein_dir / f"{protein_base}FH.pdb"
        protein_prep_out = protein_dir / protein_base
        protein_pdbqt = protein_dir / f"{protein_base}.pdbqt"

        logging.debug(f"Preparing receptor {cleaned_pdb_path}...")

        # 1. Adding Hydrogens & Optimizing (REDUCE2)
        run_reduce2(protein_path, protein_protonated)

        # 2. Preparing Receptor (Meeko)
        run_meeko_receptor(protein_protonated, protein_prep_out)
        
        prepared_results[protein_base] = protein_pdbqt
        
    logging.info("Batch cleaning complete!")
    return prepared_results

def run_scrub_ligand(ligand_path, ligand_scrubbed, ph):
    """Runs scrub.py to scrub and protonate the ligand."""
    import subprocess
    from pathlib import Path

    if Path(ligand_scrubbed).exists():
        logging.debug(f"Skipping Scrub: {ligand_scrubbed} already exists.")
        return

    cmd_scrub = [
        "scrub.py", str(ligand_path), "-o", str(ligand_scrubbed), 
        "--ph", str(ph)
    ]
    subprocess.run(cmd_scrub, check=True, capture_output=True, text=True)


def run_meeko_ligand(ligand_scrubbed, ligand_pdbqt):
    """Runs mk_prepare_ligand.py (Meeko) to prepare the ligand."""
    import subprocess
    from pathlib import Path

    if Path(ligand_pdbqt).exists():
        logging.debug(f"Skipping Meeko: {ligand_pdbqt} already exists.")
        return

    cmd_meeko = [
        "mk_prepare_ligand.py", "-i", str(ligand_scrubbed), 
        "-o", str(ligand_pdbqt)
    ]
    subprocess.run(cmd_meeko, check=True, capture_output=True, text=True)

def split_scrubbed_ligand(ligand_scrubbed, output_dir, ligand_base):
    """
    Parses the scrubbed SDF file, finds the lowest energy conformer for each isomer,
    and saves them as separate SDF files.
    """
    from rdkit import Chem
    import json
    from pathlib import Path
    
    suppl = Chem.SDMolSupplier(str(ligand_scrubbed), removeHs=False)
    
    # dictionary: isomerId -> list of (energy, mol)
    isomers = {}
    
    for mol in suppl:
        if mol is None:
            continue
            
        scrub_info_str = mol.GetProp("ScrubInfo") if mol.HasProp("ScrubInfo") else None
        if scrub_info_str:
            scrub_info = json.loads(scrub_info_str)
            isomer_id = scrub_info.get("isomerId", 0)
        else:
            isomer_id = 0
            
        if mol.HasProp("PUBCHEM_MMFF94_ENERGY"):
            energy = float(mol.GetProp("PUBCHEM_MMFF94_ENERGY"))
        else:
            energy = 0.0
                
        if isomer_id not in isomers:
            isomers[isomer_id] = []
        isomers[isomer_id].append((energy, mol))
        
    split_files = []
    out_dir = Path(output_dir)
    
    for isomer_id, confs in isomers.items():
        confs.sort(key=lambda x: x[0])  # sort by energy, lowest first
        best_mol = confs[0][1]
        
        isomer_file = out_dir / f"{ligand_base}_isomer_{isomer_id}.sdf"
        writer = Chem.SDWriter(str(isomer_file))
        writer.write(best_mol)
        writer.close()
        split_files.append(isomer_file)
        
    return split_files

def prepare_ligand(ligand_file: str, ph: float, output_dir: str):
    """
    Prepares the ligand using Molscrub and Meeko.
    Returns the path to the resulting PDBQT file.
    """
    from pathlib import Path

    ligand_path = Path(ligand_file)
    ligand_dir = ligand_path.parent
    ligand_base = ligand_path.stem

    prepared_dir = Path(output_dir)
    prepared_dir.mkdir(parents=True, exist_ok=True)

    ligand_scrubbed = prepared_dir / f"{ligand_base}_scrubbed.sdf"

    logging.debug(f"Preparing ligand {ligand_file}...")

    # 3. Scrubbing & Protonating Ligand (Molscrub)
    run_scrub_ligand(ligand_path, ligand_scrubbed, ph)

    # 3.1 Split Scrubbed Ligand by Isomer and find lowest energy conformer
    isomer_files = split_scrubbed_ligand(ligand_scrubbed, prepared_dir, ligand_base)

    # 4. Preparing Ligand (Meeko)
    pdbqt_files = []
    for isomer_file in isomer_files:
        isomer_base = isomer_file.stem
        ligand_pdbqt = prepared_dir / f"{isomer_base}.pdbqt"
        run_meeko_ligand(isomer_file, ligand_pdbqt)
        pdbqt_files.append(ligand_pdbqt)

    return pdbqt_files

def prepare_ligands(ligand_path, ph, output_dir):
    """
    Prepares all ligands in the given directory.
    Returns a dictionary mapping ligand names to their PDBQT paths.
    """
    from pathlib import Path
    prepared_ligands = {}
    for lig_file in Path(ligand_path).glob("*.sdf"):
        lig_pdbqts = prepare_ligand(str(lig_file), ph, output_dir)
        for lig_pdbqt in lig_pdbqts:
            prepared_ligands[lig_pdbqt.stem] = lig_pdbqt
    return prepared_ligands

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean PDB files by removing water and unwanted HETATMs.")
    parser.add_argument("-i", "--input_dir", type=str, default="protein", help="Input directory containing PDB files.")
    parser.add_argument("-o", "--output_dir", type=str, default="protein-clean", help="Output directory for cleaned PDB files.")
    parser.add_argument("-m", "--mode", type=str, choices=["global", "local"], default="global", help="Elimination mode: 'global' applies to all files, 'local' asks for each file individually.")
    
    args = parser.parse_args()
    
    prepare_proteins(input_dir=args.input_dir, output_dir=args.output_dir, mode=args.mode)