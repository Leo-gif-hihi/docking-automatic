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
        logging.debug(f" -> Saved: {out_filepath}")
    else:
        logging.debug(f" -> Skipped {filename}: No atoms left after cleaning.")

def run_global_mode(pdb_files, output_dir):
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
    

    for filepath, structure in structures:
        clean_and_save_pdb(structure, filepath, output_dir, sel_str)

def run_local_mode(pdb_files, output_dir):
    """Runs the cleaning process in local mode (prompt per file)."""
    logging.info(f"Running in LOCAL mode. Scanning {len(pdb_files)} files individually.")
    logging.info("Tip: You can press Ctrl+C at any prompt to escape and stop processing.\n")
    
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
        clean_and_save_pdb(structure, filepath, output_dir, sel_str)
        logging.debug("")  # Blank line for readability between files if in debug mode

def clean_proteins(input_dir="protein", output_dir="protein-clean", mode="global"):
    """Main workflow to orchestrate the cleaning process."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.debug(f"Created output directory: {output_dir}")
    
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))
    if not pdb_files:
        logging.error(f"No PDB files found in '{input_dir}'. Please check your folder.")
        return

    if mode == "global":
        run_global_mode(pdb_files, output_dir)
    elif mode == "local":
        run_local_mode(pdb_files, output_dir)
        
    logging.info("Batch cleaning complete!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean PDB files by removing water and unwanted HETATMs.")
    parser.add_argument("-i", "--input_dir", type=str, default="protein", help="Input directory containing PDB files.")
    parser.add_argument("-o", "--output_dir", type=str, default="protein-clean", help="Output directory for cleaned PDB files.")
    parser.add_argument("-m", "--mode", type=str, choices=["global", "local"], default="global", help="Elimination mode: 'global' applies to all files, 'local' asks for each file individually.")
    
    args = parser.parse_args()
    
    clean_proteins(input_dir=args.input_dir, output_dir=args.output_dir, mode=args.mode)