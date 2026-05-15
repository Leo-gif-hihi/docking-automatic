import os
import argparse
import subprocess
from pathlib import Path

from rank import rank_complexes, print_ranking

from pocket import process_pockets

def parse_args(args=None):
    """Handles argument parsing separately so it can be tested with mock lists of args."""
    parser = argparse.ArgumentParser(description="Automated Docking with AutoDock Vina")
    parser.add_argument("--protein_dir", default="protein", help="Directory containing protein files")
    parser.add_argument("--ligand_dir", default="ligand", help="Directory containing ligand SDF files")
    parser.add_argument("--box_dir", default="box", help="Directory containing box TXT files")
    parser.add_argument("--output_dir", default="output", help="Directory for output files")
    parser.add_argument("--cpus", type=int, default=0, help="Number of CPUs to use (default 0 means all CPUs)")
    parser.add_argument("--rank_only", action="store_true", help="Only rank existing results in output_dir without running docking")
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

        if not box_file.exists():
            print(f"Warning: Box file {box_file} not found for {protein_file.name}. Skipping...")
            continue

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
    """Isolated subprocess execution. Easy to mock during tests."""
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command {' '.join(cmd)}: {e}")

def main():
    args = parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    if args.rank_only:
        results = rank_complexes(args.output_dir)
        print_ranking(results, Path(args.output_dir) / "ranking.csv")
        return

    protein_path = Path(args.protein_dir)
    ligand_path = Path(args.ligand_dir)
    box_path = Path(args.box_dir)

    if not all(p.exists() for p in (protein_path, ligand_path, box_path)):
        print("Error: One or more input directories (protein, ligand, box) do not exist.")
        return

    # Phase 1: Identify pockets
    process_pockets(protein_path)

    script_path = Path(__file__).parent / "vina.sh"

    jobs = generate_docking_jobs(protein_path, ligand_path, box_path)
    
    for protein_file, ligand_file, box_file in jobs:
        print(f"\n--- Docking {ligand_file.name} to {protein_file.name} ---")
        
        cmd = build_docking_command(
            script_path, protein_file, ligand_file, box_file, args.output_dir, args.cpus
        )
        run_docking(cmd)

    # Rank complexes after all docking jobs are complete
    results = rank_complexes(args.output_dir)
    print_ranking(results, Path(args.output_dir) / "ranking.csv")

if __name__ == "__main__":
    main()
