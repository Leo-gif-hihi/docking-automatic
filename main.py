import os
import argparse
import subprocess
import logging
import time
from pathlib import Path

from rank import rank_complexes, print_ranking

from pocket import process_pockets

from clean_protein import clean_proteins

from download_protein import download_protein_biopython
from Bio.PDB import PDBList

def parse_args(args=None):
    """Handles argument parsing separately so it can be tested with mock lists of args."""
    parser = argparse.ArgumentParser(description="Automated Docking with AutoDock Vina")
    parser.add_argument("--protein_list", default="protein.txt", help="Text file containing list of protein PDB codes")
    parser.add_argument("--ligand_dir", default="ligand", help="Directory containing ligand SDF files")
    parser.add_argument("--box_dir", default="box", help="Directory containing box TXT files")
    parser.add_argument("--output_dir", default="output", help="Directory for output files")
    parser.add_argument("--cpus", type=int, default=0, help="Number of CPUs to use (default 0 means all CPUs)")
    parser.add_argument("--rank_only", action="store_true", help="Only rank existing results in output_dir without running docking")
    parser.add_argument("--skip_autopoc", action="store_true", help="Skip automatic pocket identification and use existing provided box files")
    parser.add_argument("--identify_pockets_only", action="store_true", help="Only identify pockets and exit; skip docking and ranking")
    parser.add_argument("--skip_clean", action="store_true", help="Skip cleaning protein structures")
    parser.add_argument("--clean_mode", type=str, choices=["global", "local"], default="global", help="Elimination mode for cleaning protein structures")
    return parser.parse_args(args)

def setup_logging(output_dir):
    """Configures logging for console and file."""
    logger = logging.getLogger()
    # Clear any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()
        
    logger.setLevel(logging.DEBUG)

    # 1. File Handler: Always log everything to file
    file_handler = logging.FileHandler(os.path.join(output_dir, "autodock_run.log"))
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_format)
    logger.addHandler(file_handler)

    # 2. Console Handler: Always DEBUG
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

def generate_docking_jobs(protein_path, ligand_path, box_path):
    """
    Generator that creates valid combinations of protein, ligand, and box files.
    Yields: (protein_file, ligand_file, box_file)
    """
    protein_files = list(protein_path.glob("*.pdb"))
    
    if not protein_files:
        logging.error(f"No .pdb files found in {protein_path}.")
        return

    for protein_file in protein_files:
        protein_base = protein_file.stem
        box_file = box_path / f"{protein_base}.box.txt"

        for ligand_file in ligand_path.glob("*.sdf"):
            yield protein_file, ligand_file, box_file

def build_docking_command(script_path, protein_file, ligand_file, box_file, output_dir, cpus):
    """Pure logic function for building the command. Easy to unit test."""
    return [
        "bash",
        str(script_path),
        str(protein_file),
        str(ligand_file),
        str(box_file),
        str(output_dir),
        str(cpus)
    ]

def run_docking(cmd):
    """Isolated subprocess execution. Returns True on success, False on failure."""
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.debug(result.stdout)
        if result.stderr:
            logging.debug(result.stderr)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing command {' '.join(cmd)}")
        if e.stdout:
            logging.debug(f"STDOUT: {e.stdout}")
        if e.stderr:
            logging.debug(f"STDERR: {e.stderr}")
        return False

def main():
    args = parse_args()

    if not args.rank_only:
        existing_dirs = [d for d in [args.output_dir, "protein", "protein-clean"] if os.path.exists(d)]
        if existing_dirs:
            print(f"\n\033[1;33m[WARNING] The following directories already exist: {', '.join(existing_dirs)}\033[0m")
            print("\033[1;33mOld results in these folders can conflict with the current pipeline.\033[0m")
            print(
                "\033[1;33mIf you want to use the previous results and ensure that all files in these directories "
                "are relevant to your project (for example, rerun the workflow after interrupted), "
                "feel free to ignore this warning.\033[0m"
            )
            ans = input("Do you want to delete them before continuing? (y/N): ").strip().lower()
            if ans == 'y':
                import shutil
                for d in existing_dirs:
                    shutil.rmtree(d, ignore_errors=True)
                print("\033[1;32mDirectories deleted.\033[0m\n")
            else:
                print("Continuing without deleting...\n")

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir)

    if args.rank_only:
        step_start = time.time()
        logging.info("\n\033[1;32m[WORKFLOW] Ranking complexes only...\033[0m")
        results = rank_complexes(args.output_dir)
        print_ranking(results, Path(args.output_dir) / "ranking.csv")
        logging.info(f"\033[1;36m[TIME] Step duration: {time.time() - step_start:.2f} seconds\033[0m")
        return

    protein_path = Path("protein")
    protein_clean_dir = "protein-clean"
    
    os.makedirs(protein_path, exist_ok=True)
    try:
        step_start = time.time()
        logging.info("\n\033[1;32m[WORKFLOW] Starting protein download process...\033[0m")
        with open(args.protein_list, 'r') as file:
            protein_codes = [line.strip().upper() for line in file if line.strip()]
            logging.debug(f"Protein codes to download: {protein_codes}")

        pdbl = PDBList(verbose=False)
        for code in protein_codes:
            download_protein_biopython(code, pdbl, str(protein_path))
        logging.info(f"\033[1;36m[TIME] Step duration: {time.time() - step_start:.2f} seconds\033[0m")
    except FileNotFoundError:
        logging.error(f"Error: File '{args.protein_list}' not found.")
        return

    ligand_path = Path(args.ligand_dir)
    box_path = Path(args.box_dir)
    # Identify pockets (either as part of pipeline or as a dedicated identify-only run)
    vis_dir = Path(args.output_dir) / "visualization_pocket"
    os.makedirs(vis_dir, exist_ok=True)
    if args.identify_pockets_only or not args.skip_autopoc:
        step_start = time.time()
        logging.info("\n\033[1;32m[WORKFLOW] Identifying pockets...\033[0m")
        if args.identify_pockets_only and args.skip_autopoc:
            logging.warning("Both --skip_autopoc and --identify_pockets_only provided; ignoring --skip_autopoc and running pocket identification.")
        process_pockets(protein_path, box_path, output_dir=str(vis_dir))
        logging.info(f"\033[1;36m[TIME] Step duration: {time.time() - step_start:.2f} seconds\033[0m")
        if args.identify_pockets_only:
            logging.info("\n\033[1;32m[WORKFLOW] Pocket identification complete. Exiting (identify-only mode).\033[0m")
            return

    if not all(p.exists() for p in (protein_path, ligand_path, box_path)):
        logging.error("Error: One or more input directories (protein, ligand, box) do not exist.")
        return
    
    # Clean proteins
    if not args.skip_clean:
        step_start = time.time()
        logging.info("\n\033[1;32m[WORKFLOW] Cleaning proteins...\033[0m")
        clean_proteins(input_dir=str(protein_path), output_dir=protein_clean_dir, mode=args.clean_mode)
        protein_path = Path(protein_clean_dir)
        logging.info(f"\033[1;36m[TIME] Step duration: {time.time() - step_start:.2f} seconds\033[0m")
    else:
        # If skip_clean is provided, we might still want to use protein_clean_dir if it has files, or protein_path.
        # We'll use protein_path unless protein_clean_dir is specifically needed, but typically skip_clean means we use protein_path directly.
        pass


    script_path = Path(__file__).parent / "vina.sh"

    jobs_list = list(generate_docking_jobs(protein_path, ligand_path, box_path))
    total_jobs = len(jobs_list)
    error_jobs = []
    step_start = time.time()
    logging.info(f"\n\033[1;32m[WORKFLOW] Generated docking jobs. Starting docking process for {total_jobs} combinations...\033[0m")
    for i, (protein_file, ligand_file, box_file) in enumerate(jobs_list, 1):
        logging.info(f"Completed {i}/{total_jobs} complexes")
        logging.debug(f"\n--- Docking {ligand_file.name} to {protein_file.name} ---")
        
        protein_base = protein_file.stem
        ligand_base = ligand_file.stem
        vina_out_dir = Path(args.output_dir) / "vina_output"
        complex_output_dir = vina_out_dir / f"{protein_base}_{ligand_base}"
        os.makedirs(complex_output_dir, exist_ok=True)

        # CHECK IF VINA LOG ALREADY EXISTS
        expected_log_file = complex_output_dir / f"{protein_base}_{ligand_base}_vina.log"
        if expected_log_file.exists():
            logging.debug(f"Log file {expected_log_file.name} already exists. Skipping docking...")
            continue
        
        if not box_file.exists():
            logging.warning(f"Error: Box file {box_file.name} not found. Skipping docking...")
            error_jobs.append(f"{protein_file.name}\t{ligand_file.name}\tMissing Box File")
            continue
            
        cmd = build_docking_command(
            script_path, protein_file, ligand_file, box_file, str(complex_output_dir), args.cpus
        )
        success = run_docking(cmd)
        if not success:
            error_jobs.append(f"{protein_file.name}\t{ligand_file.name}\tDocking Failed")

    if error_jobs:
        error_file = Path(args.output_dir) / "error_complexes.txt"
        with open(error_file, "w") as ef:
            ef.write("Protein\tLigand\tError_Type\n")
            ef.write("\n".join(error_jobs) + "\n")
        logging.warning(f"\nRecorded {len(error_jobs)} failed docking jobs in {error_file}")

    logging.info(f"\033[1;36m[TIME] Step duration: {time.time() - step_start:.2f} seconds\033[0m")

    # Rank complexes after all docking jobs are complete
    step_start = time.time()
    logging.info("\n\033[1;32m[WORKFLOW] Docking complete. Generating ranking...\033[0m")
    results = rank_complexes(args.output_dir)
    print_ranking(results, Path(args.output_dir) / "ranking.csv")
    logging.info(f"\033[1;36m[TIME] Step duration: {time.time() - step_start:.2f} seconds\033[0m")

if __name__ == "__main__":
    main()
