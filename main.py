import os
import argparse
import subprocess
import logging
import time
from pathlib import Path

from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from logger_utils import log_step, console

from rank import rank_complexes, print_ranking, generate_complexes

from pocket import process_pockets

from prepare import prepare_proteins, prepare_ligands

from download_protein import download_proteins

from p2rank import process_unprocessed_with_p2rank

def parse_args(args=None):
    """Handles argument parsing separately so it can be tested with mock lists of args."""
    parser = argparse.ArgumentParser(description="Automated Docking with AutoDock Vina")
    parser.add_argument("--protein_list", default="protein.txt", help="Text file containing list of protein PDB codes")
    parser.add_argument("--ligand_dir", default="ligand", help="Directory containing ligand SDF files")
    parser.add_argument("--box_dir", default=None, help="Directory containing box TXT files (defaults: box_{protein_list})")
    parser.add_argument("--output_dir", default=None, help="Directory for output files (default: output_{protein_list}_{ligand_dir})")
    parser.add_argument("--cpus", type=int, default=0, help="Number of CPUs to use (default 0 means all CPUs)")
    parser.add_argument("--skip_autopoc", action="store_true", help="Skip automatic pocket identification and use existing provided box files")
    parser.add_argument("--dock_all_pockets", action="store_true", help="Generate box files for all pockets and dock ligands into all of them. Default is to only use the highest-scoring pocket.")
    parser.add_argument("--identify_pockets_only", action="store_true", help="Only identify pockets and exit; skip docking and ranking")
    parser.add_argument("--clean_mode", type=str, choices=["auto","global", "local"], default="auto", help="Elimination mode for cleaning protein structures")
    parser.add_argument("--ph", type=float, default=7.4, help="pH value to prepare ligands (default: 7.4)")
    parser.add_argument("--skip_cofactor", action="store_true", default=False, help="Delete all HETATM records instead of only water (default: keep cofactors, only delete water)")
    parser.add_argument("--generate_isomers", action="store_true", help="Generate acid-base and tautomer isomers during ligand preparation (default is to skip)")
    parser.add_argument("--skip_minimization", action="store_true", help="Skip energy minimization step")
    parser.add_argument("--num_runs", type=int, default=3, help="Number of independent docking runs per complex (default: 3)")
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

    # 2. Console Handler: WARNING level for terminal (UI is handled by log_step)
    console_handler = RichHandler(rich_tracebacks=True, show_time=False, show_path=False)
    console_handler.setLevel(logging.WARNING)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

def generate_docking_jobs(prepared_proteins, prepared_ligands, box_path, num_runs=3):
    """
    Generator that creates valid combinations of protein, ligand, and box files.
    Yields: (protein_base, ligand_base, box_file, run_index)
    """
    if not prepared_proteins:
        logging.error(f"No prepared proteins found.")
        return

    for protein_base in prepared_proteins:
        box_files = list(box_path.glob(f"{protein_base}.box.txt")) + list(box_path.glob(f"{protein_base}_pocket_*.box.txt"))
        if not box_files:
            logging.warning(f"No box files found for {protein_base} in {box_path}")
            continue

        for box_file in box_files:
            for ligand_base in prepared_ligands:
                for run_index in range(1, num_runs + 1):
                    yield protein_base, ligand_base, box_file, run_index

def run_docking_pipeline(protein_pdbqt, ligand_pdbqt, box_file, output_dir, protein_base, ligand_base, run_index, cpus):
    """Runs the docking pipeline for a single pair (already prepared)."""
    try:
        out_pdbqt = Path(output_dir) / f"{protein_base}_{ligand_base}_vina_out.pdbqt"
        out_sdf = Path(output_dir) / f"{protein_base}_{ligand_base}_vina_out.sdf"
        out_log = Path(output_dir) / f"{protein_base}_{ligand_base}_vina.log"

        # 5. Running AutoDock Vina
        cmd_vina = [
            "vina", "--receptor", str(protein_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--config", str(box_file),
            "--out", str(out_pdbqt)
        ]
        if cpus > 0:
            cmd_vina.extend(["--cpu", str(cpus)])
            
        logging.debug(f"Running Vina: {' '.join(cmd_vina)}")
        
        with open(out_log, "w") as log_file:
            subprocess.run(cmd_vina, check=True, stdout=log_file, stderr=subprocess.STDOUT)

        # 6. Exporting Results to SDF (Meeko)
        cmd_export = ["mk_export.py", str(out_pdbqt), "-s", str(out_sdf)]
        subprocess.run(cmd_export, check=True, capture_output=True, text=True)

        return True
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing command: {e.cmd}")
        if e.stdout:
            logging.debug(f"STDOUT: {e.stdout}")
        if e.stderr:
            logging.debug(f"STDERR: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error during docking pipeline: {e}")
        return False

def main():
    args = parse_args()

    # Dynamically set protein_path, protein_clean_dir, output_dir, and box_dir based on protein_list and ligand_dir
    protein_list_name = Path(args.protein_list).stem
    ligand_dir_name = Path(args.ligand_dir).stem
    protein_path = Path(protein_list_name)
    protein_clean_dir = f"{protein_list_name}_prepared"
    ligand_prepared_dir = f"{ligand_dir_name}_prepared"

    if args.output_dir is None:
        args.output_dir = f"output_{protein_list_name}_{ligand_dir_name}"

    if args.box_dir is None:
        args.box_dir = f"box_{protein_list_name}"

    existing_dirs = [d for d in [args.output_dir, str(protein_path), protein_clean_dir, ligand_prepared_dir] if os.path.exists(d)]
    if existing_dirs:
        print()
        log_step("WARNING", f"The following directories already exist: {', '.join(existing_dirs)}", color="yellow")
        log_step("WARNING", "Old results in these folders can conflict with the current pipeline.", color="yellow")
        log_step(
            "WARNING",
            "If you want to use the previous results and ensure that all files in these directories "
            "are relevant to your project (for example, rerun the workflow after interrupted), "
            "feel free to ignore this warning.", color="yellow"
        )
        ans = input("Do you want to delete them before continuing? (y/N): ").strip().lower()
        if ans == 'y':
            import shutil
            for d in existing_dirs:
                shutil.rmtree(d, ignore_errors=True)
            log_step("INFO", "Directories deleted.")
            print()
        else:
            print("Continuing without deleting...\n")

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir)

    os.makedirs(protein_path, exist_ok=True)
    try:
        step_start = time.time()
        print()
        log_step("WORKFLOW", "Starting protein download process...")
        download_proteins(args.protein_list, str(protein_path))
        log_step("TIME", f"Step duration: {time.time() - step_start:.2f} seconds", color="cyan")
    except FileNotFoundError:
        return

    ligand_path = Path(args.ligand_dir)
    box_path = Path(args.box_dir)
    if not all(p.exists() for p in (protein_path, ligand_path)):
        logging.error("Error: One or more input directories (protein, ligand) do not exist.")
        return
    
    # Clean and prepare proteins
    prepared_proteins = {}
    step_start = time.time()
    print()
    log_step("WORKFLOW", "Preparing proteins...")
    prepared_proteins = prepare_proteins(input_dir=str(protein_path), output_dir=protein_clean_dir, mode=args.clean_mode, skip_cofactor=args.skip_cofactor, skip_minimization=args.skip_minimization)
    
    protein_clean_path = Path(protein_clean_dir)
    log_step("TIME", f"Step duration: {time.time() - step_start:.2f} seconds", color="cyan")

    # Prepare ligands
    step_start = time.time()
    print()
    log_step("WORKFLOW", "Preparing ligands...")
    
    prepared_ligands = prepare_ligands(ligand_path, args.ph, ligand_prepared_dir, generate_isomers=args.generate_isomers)
        
    log_step("TIME", f"Step duration: {time.time() - step_start:.2f} seconds", color="cyan")

    # Identify pockets (either as part of pipeline or as a dedicated identify-only run)
    vis_dir = Path(args.output_dir) / "visualization_pocket"
    os.makedirs(vis_dir, exist_ok=True)
    if args.identify_pockets_only or not args.skip_autopoc:
        step_start = time.time()
        print()
        log_step("WORKFLOW", "Identifying pockets...")
        if args.identify_pockets_only and args.skip_autopoc:
            logging.warning("Both --skip_autopoc and --identify_pockets_only provided; ignoring --skip_autopoc and running pocket identification.")
        unprocessed_list = process_pockets(protein_clean_path, box_path, output_dir=str(vis_dir), dock_all_pockets=args.dock_all_pockets)
        
        if args.dock_all_pockets:
            print()
            log_step("INTERACTIVE", f"Pocket identification complete. Check {vis_dir}/pocket_reliability.csv.", color="cyan")
            eliminated = input("Enter pocket IDs to eliminate separated by commas (e.g. 1, 3, 5), or press Enter to keep all: ").strip()
            if eliminated:
                eliminated_ids = [pid.strip() for pid in eliminated.split(',') if pid.strip()]
                for pid in eliminated_ids:
                    for box_file in box_path.glob(f"*_pocket_{pid}.box.txt"):
                        logging.info(f"Eliminating pocket: {box_file.name}")
                        try:
                            box_file.unlink()
                        except FileNotFoundError:
                            pass

                # Fallback to p2rank for proteins that have no pockets left
                if unprocessed_list is None:
                    unprocessed_list = str(protein_path / "unprocessed_proteins.txt")
                
                existing_unprocessed = set()
                if os.path.exists(unprocessed_list):
                    with open(unprocessed_list, "r", encoding="utf-8") as f:
                        existing_unprocessed = set(line.strip() for line in f if line.strip())

                for protein_file in protein_clean_path.glob("*FH.cif"):
                    protein_base = protein_file.stem[:-2] if protein_file.stem.endswith('FH') else protein_file.stem
                    remaining_boxes = list(box_path.glob(f"{protein_base}*.box.txt"))
                    if not remaining_boxes and protein_file.name not in existing_unprocessed:
                        with open(unprocessed_list, "a", encoding="utf-8") as f:
                            f.write(f"{protein_file.name}\n")
                        existing_unprocessed.add(protein_file.name)
                        logging.info(f"Protein {protein_base} has no pockets left after elimination. Added to unprocessed list for p2rank fallback.")

        # Backup pocket identification using p2rank
        if unprocessed_list and os.path.exists(unprocessed_list):
            print()
            log_step("WORKFLOW", "Running p2rank as a backup for unprocessed proteins...")
            process_unprocessed_with_p2rank(str(unprocessed_list), str(protein_clean_path), str(box_path), vis_dir=str(vis_dir))
            
        log_step("TIME", f"Step duration: {time.time() - step_start:.2f} seconds", color="cyan")
        if args.identify_pockets_only:
            print()
            log_step("WORKFLOW", "Workflow complete. Exiting (identify-only mode).")
            return

    if not box_path.exists() and not args.skip_autopoc:
        logging.error("Error: Box directory does not exist after pocket identification.")
        return

    jobs_list = list(generate_docking_jobs(prepared_proteins, prepared_ligands, box_path, args.num_runs))
    total_jobs = len(jobs_list)
    error_jobs = []
    step_start = time.time()
    print()
    log_step("WORKFLOW", f"Generated docking jobs. Starting docking process for {total_jobs} combinations...")
    
    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TextColumn("•"),
        TimeElapsedColumn(),
        TextColumn("•"),
        TimeRemainingColumn(),
        expand=True
    ) as progress:
        docking_task = progress.add_task("[cyan]Starting docking...", total=total_jobs)
        for i, (protein_base, ligand_base, box_file, run_index) in enumerate(jobs_list, 1):
            try:
                protein_pocket_base = box_file.name.replace(".box.txt", "")
                progress.update(docking_task, description=f"[cyan]Docking [bold]{protein_pocket_base}[/bold] & [bold]{ligand_base}[/bold] Run {run_index} ({i}/{total_jobs})")
                logging.debug(f"\n--- Docking {ligand_base} to {protein_pocket_base} (Run {run_index}) ---")
                
                vina_out_dir = Path(args.output_dir) / "vina_output" / f"run_{run_index}"
                complex_output_dir = vina_out_dir / f"{protein_pocket_base}_{ligand_base}"
                os.makedirs(complex_output_dir, exist_ok=True)
        
                # CHECK IF VINA LOG ALREADY EXISTS AND IS COMPLETE
                expected_log_file = complex_output_dir / f"{protein_pocket_base}_{ligand_base}_vina.log"
                if expected_log_file.exists():
                    with open(expected_log_file, "r") as log_file:
                        log_content = log_file.read()
                        if log_content.count("*") >= 51:
                            logging.debug(f"Log file {expected_log_file.name} already exists and docking is complete. Skipping docking...")
                            continue
                        else:
                            logging.debug(f"Log file {expected_log_file.name} exists but docking is incomplete. Proceeding with docking...")
                    
                if not box_file.exists():
                    logging.warning(f"Error: Box file {box_file.name} not found. Skipping docking...")
                    error_jobs.append(f"{protein_base}\t{ligand_base}\tRun {run_index}\tMissing Box File")
                    continue
                    
                protein_pdbqt = prepared_proteins.get(protein_base)
                ligand_pdbqt = prepared_ligands.get(ligand_base)
                
                if not protein_pdbqt or not ligand_pdbqt:
                    logging.error(f"Error: Missing prepared files for {protein_base} or {ligand_base}. Skipping...")
                    error_jobs.append(f"{protein_base}\t{ligand_base}\tRun {run_index}\tPreparation Failed")
                    continue
                    
                success = run_docking_pipeline(
                    protein_pdbqt, ligand_pdbqt, box_file, str(complex_output_dir), 
                    protein_pocket_base, ligand_base, run_index, args.cpus
                )
                if not success:
                    error_jobs.append(f"{protein_base}\t{ligand_base}\tRun {run_index}\tDocking Failed")
            finally:
                progress.advance(docking_task)

    if error_jobs:
        error_file = Path(args.output_dir) / "error_complexes.txt"
        with open(error_file, "w") as ef:
            ef.write("Protein\tLigand\tRun\tError_Type\n")
            ef.write("\n".join(error_jobs) + "\n")
        logging.warning(f"\nRecorded {len(error_jobs)} failed docking jobs in {error_file}")

    log_step("TIME", f"Step duration: {time.time() - step_start:.2f} seconds", color="cyan")

    # Rank complexes after all docking jobs are complete
    step_start = time.time()
    print()
    log_step("WORKFLOW", "Docking complete. Generating ranking...")
    results = rank_complexes(args.output_dir, list(prepared_ligands.keys()) if prepared_ligands else None)
    curated_results = print_ranking(results, Path(args.output_dir) / "ranking.csv")
    
    log_step("WORKFLOW", "Generating complex files for top results...")
    if curated_results:
        generate_complexes(curated_results, args.output_dir, protein_clean_dir, display_limit=20)
    else:
        logging.warning("No valid curated results to generate complexes from.")
    
    log_step("TIME", f"Step duration: {time.time() - step_start:.2f} seconds", color="cyan")

if __name__ == "__main__":
    main()
