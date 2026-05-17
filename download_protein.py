import os
from Bio.PDB import PDBList

def download_protein_biopython(protein_code, pdbl, download_dir):
    """
    Downloads a single PDB file using Biopython's PDBList, 
    skipping the download if the file already exists.
    """
    # 1. Define what the final expected filename will look like
    expected_filename = os.path.join(download_dir, f"{protein_code.lower()}.pdb")
    
    # 2. Check if this file already exists in your target folder
    if os.path.exists(expected_filename):
        print(f"Protein {protein_code} already exists. Skipping download.")
        return # Exit the function early without downloading

    try:
        # Fetch the standard PDB file (changed overwrite to False)
        filename = pdbl.retrieve_pdb_file(
            pdb_code=protein_code,
            pdir=download_dir,
            file_format='pdb',
            overwrite=False 
        )
        
        # Rename "pdbXXXX.ent" to "XXXX.pdb"
        if filename and filename.endswith(".ent"):
            os.rename(filename, expected_filename)
            print(f"Successfully downloaded and saved: {expected_filename}")
        else:
            print(f"Successfully downloaded {protein_code} to {filename}")

    except Exception as e:
        print(f"Failed to download protein {protein_code}. Error: {e}")

def main():
    cwd = os.getcwd()
    protein_file_path = os.path.join(cwd, "protein.txt")
    
    # Define and create the "protein" directory
    protein_dir = os.path.join(cwd, "protein")
    os.makedirs(protein_dir, exist_ok=True)
    
    print(f"Current working directory: {cwd}")
    print(f"Protein file path: {protein_file_path}")

    try:
            # Read the PDB codes from the file
            with open(protein_file_path, 'r') as file:
                protein_codes = [line.strip().upper() for line in file if line.strip()]
                print(f"Protein codes to download: {protein_codes}")

            # Initialize Biopython's PDB downloader tool
            pdbl = PDBList(verbose=False)

            # Loop through the list sequentially
            for code in protein_codes:
                download_protein_biopython(code, pdbl, protein_dir)

    except FileNotFoundError:
        print(f"Error: File 'protein.txt' not found at {protein_file_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()