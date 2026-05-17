import os
import argparse
import subprocess
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
    parser.add_argument("--skip_clean", action="store_true", help="Skip cleaning protein structures")
    parser.add_argument("--clean_mode", type=str, choices=["global", "local"], default="global", help="Elimination mode for cleaning protein structures")
    return parser.parse_args(args)

def generate_docking_jobs(protein_path, ligand_path, box_path):
    """
    Generator that creates valid combinations of protein, ligand, and box files.
    Yields: (protein_file, ligand_file, box_file)
    """
    protein_files = list(protein_path.glob("*.pdb"))
    
    if not protein_files:
        print(f"No .pdb files found in {protein_path}.")
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
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command {' '.join(cmd)}: {e}")
        return False

def main():
    args = parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    if args.rank_only:
        results = rank_complexes(args.output_dir)
        print_ranking(results, Path(args.output_dir) / "ranking.csv")
        return

    protein_path = Path("protein")
    protein_clean_dir = "protein-clean"
    
    os.makedirs(protein_path, exist_ok=True)
    try:
        with open(args.protein_list, 'r') as file:
            protein_codes = [line.strip().upper() for line in file if line.strip()]
            print(f"Protein codes to download: {protein_codes}")

        pdbl = PDBList(verbose=False)
        for code in protein_codes:
            download_protein_biopython(code, pdbl, str(protein_path))
    except FileNotFoundError:
        print(f"Error: File '{args.protein_list}' not found.")
        return

    ligand_path = Path(args.ligand_dir)
    box_path = Path(args.box_dir)
    # Identify pockets
    vis_dir = Path(args.output_dir) / "visualization_pocket"
    os.makedirs(vis_dir, exist_ok=True)
    if not args.skip_autopoc:
        process_pockets(protein_path, box_path, output_dir=str(vis_dir))

    if not all(p.exists() for p in (protein_path, ligand_path, box_path)):
        print("Error: One or more input directories (protein, ligand, box) do not exist.")
        return
    
    # Clean proteins
    if not args.skip_clean:
        clean_proteins(input_dir=str(protein_path), output_dir=protein_clean_dir, mode=args.clean_mode)
        protein_path = Path(protein_clean_dir)
    else:
        # If skip_clean is provided, we might still want to use protein_clean_dir if it has files, or protein_path.
        # We'll use protein_path unless protein_clean_dir is specifically needed, but typically skip_clean means we use protein_path directly.
        pass


    script_path = Path(__file__).parent / "vina.sh"

    jobs = generate_docking_jobs(protein_path, ligand_path, box_path)
    error_jobs = []
    
    for protein_file, ligand_file, box_file in jobs:
        print(f"\n--- Docking {ligand_file.name} to {protein_file.name} ---")
        
        protein_base = protein_file.stem
        ligand_base = ligand_file.stem
        vina_out_dir = Path(args.output_dir) / "vina_output"
        complex_output_dir = vina_out_dir / f"{protein_base}_{ligand_base}"
        os.makedirs(complex_output_dir, exist_ok=True)
        
        if not box_file.exists():
            print(f"Error: Box file {box_file.name} not found. Skipping docking...")
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
        print(f"\nRecorded {len(error_jobs)} failed docking jobs in {error_file}")

    # Rank complexes after all docking jobs are complete
    results = rank_complexes(args.output_dir)
    print_ranking(results, Path(args.output_dir) / "ranking.csv")

if __name__ == "__main__":
    main()
