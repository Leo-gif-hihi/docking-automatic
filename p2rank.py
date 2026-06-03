import os
import subprocess
import logging
import csv
import shutil
from pathlib import Path

def get_protein_base(protein):
    if protein.endswith('FH.pdb'):
        return protein[:-6]
    elif protein.endswith('.pdb') or protein.endswith('.cif'):
        return protein[:-4]
    return protein

def read_unprocessed_list(file_path):
    """
    Reads and parses the list of unprocessed proteins.
    Returns a list of protein filenames.
    """
    if not os.path.exists(file_path):
        logging.warning(f"Unprocessed file not found: {file_path}")
        return []

    with open(file_path, 'r', encoding='utf-8') as f:
        proteins = [line.strip() for line in f if line.strip()]
        
    return proteins

def create_p2rank_dataset(proteins, prepared_dir, output_ds_path):
    """
    Validates protein files and creates a p2rank dataset (.ds) file.
    Returns True if the dataset was created with at least one entry, False otherwise.
    """
    prepared_path = Path(prepared_dir)
    if not prepared_path.exists():
        logging.error(f"Prepared directory not found: {prepared_dir}")
        return False

    valid_proteins_count = 0
    with open(output_ds_path, 'w', encoding='utf-8') as f:
        for protein in proteins:
            protein_path = prepared_path / protein
            if protein_path.exists():
                f.write(f"{protein_path.absolute()}\n")
                valid_proteins_count += 1
            else:
                logging.warning(f"Protein file missing in prepared directory: {protein_path}")

    if valid_proteins_count == 0:
        logging.info("No valid protein files found to add to the dataset.")
        # Clean up the empty dataset file if no valid proteins were found
        if os.path.exists(output_ds_path):
            os.remove(output_ds_path)
        return False
        
    logging.info(f"Created p2rank dataset file with {valid_proteins_count} entries: {output_ds_path}")
    return True

def run_prank_predict(ds_file_path, output_dir):
    """
    Executes the external prank predict command using the provided dataset.
    """
    prank_path = shutil.which("prank")
    if not prank_path:
        logging.error("Command 'prank' not found in PATH. Please ensure p2rank is installed.")
        return False
        
    # Resolve symlinks because the prank bash script relies on its own location 
    # for resolving the java classpath, which breaks if called via a symlink.
    real_prank_path = os.path.realpath(prank_path)
    
    cmd = [
        real_prank_path,
        "predict",
        str(Path(ds_file_path).absolute()),
        "-o",
        str(Path(output_dir).absolute())
    ]
    
    logging.info(f"Running p2rank command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"p2rank prediction completed successfully. Output saved to: {output_dir}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing p2rank: {e}")
        return False
    except FileNotFoundError:
        logging.error("Command 'prank' not found. Please ensure p2rank is installed and the 'prank' command is in your PATH.")
        return False

def create_box_from_predictions(p2rank_out_dir, proteins, box_dir, prepared_dir, size=30.0, vis_dir=None):
    """
    Extracts the center coordinates of the top 1 pocket from p2rank predictions
    and creates a Vina box configuration file.
    """
    from pocket import generate_pymol_box_script
    box_path = Path(box_dir)
    box_path.mkdir(exist_ok=True, parents=True)
    p2rank_path = Path(p2rank_out_dir)

    for protein in proteins:
        csv_file = p2rank_path / f"{protein}_predictions.csv"
        if not csv_file.exists():
            logging.warning(f"Predictions file not found for {protein}: {csv_file}")
            continue

        try:
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, skipinitialspace=True)
                for row in reader:
                    # Look for rank == 1 or name == pocket1
                    if row.get('rank', '').strip() == '1' or row.get('name', '').strip() == 'pocket1':
                        center_x = float(row['center_x'])
                        center_y = float(row['center_y'])
                        center_z = float(row['center_z'])
                        
                        # Determine base name for the box file
                        protein_base = get_protein_base(protein)

                        box_filename = box_path / f"{protein_base}.box.txt"
                        with open(box_filename, 'w', encoding='utf-8') as bf:
                            bf.write(f"center_x = {center_x:.3f}\n")
                            bf.write(f"center_y = {center_y:.3f}\n")
                            bf.write(f"center_z = {center_z:.3f}\n")
                            bf.write(f"size_x = {size:.3f}\n")
                            bf.write(f"size_y = {size:.3f}\n")
                            bf.write(f"size_z = {size:.3f}\n")
                            bf.write(f"exhaustiveness = 32\n") # Default for size=30^3
                        
                        logging.info(f"Created box file for {protein} at {box_filename}")
                        
                        box_params = {
                            'center_x': center_x, 'center_y': center_y, 'center_z': center_z,
                            'size_x': size, 'size_y': size, 'size_z': size
                        }
                        prepared_path = Path(prepared_dir)
                        protein_file_path = prepared_path / protein
                        
                        pymol_out_dir = vis_dir if vis_dir else box_dir
                        generate_pymol_box_script(
                            protein_file=Path(protein), 
                            box_params=box_params, 
                            output_dir=pymol_out_dir, 
                            suffix="_p2rank",
                            pdb_to_load=protein_file_path.absolute()
                        )
                        break # Only process the top rank
        except Exception as e:
            logging.error(f"Failed to extract predictions for {protein}: {e}")

def process_unprocessed_with_p2rank(unprocessed_file, prepared_dir, box_dir, vis_dir=None):
    """
    Orchestrates the workflow for processing failed proteins using p2rank.
    """
    # 1. Read the list of proteins
    all_proteins = read_unprocessed_list(unprocessed_file)
    if not all_proteins:
        logging.debug("No unprocessed proteins found in the list. Nothing to run for p2rank.")
        return

    # Filter out proteins that already have a box file
    box_path = Path(box_dir)
    proteins = []
    for p in all_proteins:
        protein_base = get_protein_base(p)
        box_file = box_path / f"{protein_base}.box.txt"
        if not box_file.exists():
            proteins.append(p)
        else:
            logging.debug(f"Box file already exists for {p} ({box_file.name}), skipping p2rank.")

    if not proteins:
        logging.info("All unprocessed proteins already have box files. Skipping p2rank workflow.")
        return

    # 2. Create the dataset file
    prepared_path = Path(prepared_dir)
    ds_file = prepared_path / "p2rank_proteins.ds"
    output_dir = prepared_path / "p2rank"
    
    has_valid_dataset = create_p2rank_dataset(proteins, prepared_dir, ds_file)
    if not has_valid_dataset:
        return

    # 3. Execute the prediction
    success = run_prank_predict(ds_file, output_dir)
    
    # 4. Extract predictions and create box files
    if success:
        create_box_from_predictions(output_dir, proteins, box_dir, prepared_dir, vis_dir=vis_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run p2rank on unprocessed proteins")
    parser.add_argument("--unprocessed_list", required=True, help="Path to unprocessed_proteins.txt")
    parser.add_argument("--prepared_dir", required=True, help="Path to the prepared directory containing the FH pdb files")
    parser.add_argument("--box_dir", required=True, help="Path to save the generated box files")
    args = parser.parse_args()
    
    # Configure logging if running as a script
    logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
    
    process_unprocessed_with_p2rank(args.unprocessed_list, args.prepared_dir, args.box_dir)
