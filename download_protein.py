import os
import logging
from Bio.PDB import PDBList
import time

def download_protein_biopython(protein_code, pdbl, download_dir):
    """
    Downloads a single PDB file using Biopython's PDBList, 
    skipping the download if the file already exists.
    """
    # 1. Define what the final expected filename will look like
    expected_filename = os.path.join(download_dir, f"{protein_code.lower()}.pdb")
    
    # 2. Check if this file already exists in your target folder
    if os.path.exists(expected_filename):
        logging.debug(f"Protein {protein_code} already exists. Skipping download.")
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
            logging.debug(f"Successfully downloaded and saved: {expected_filename}")
        else:
            logging.debug(f"Successfully downloaded {protein_code} to {filename}")

    except Exception as e:
        logging.error(f"Failed to download protein {protein_code}. Error: {e}")

def download_proteins(protein_list_file, download_dir):
    """
    Reads a list of PDB codes from a file and downloads them sequentially.
    Retries up to 3 times for each protein in case of errors, with a 10-second delay between attempts.
    """
    try:
        with open(protein_list_file, 'r') as file:
            protein_codes = [line.strip().upper() for line in file if line.strip()]
            logging.debug(f"Protein codes to download: {protein_codes}")

        # Initialize Biopython's PDB downloader tool
        pdbl = PDBList(verbose=False)

        # Loop through the list sequentially
        for code in protein_codes:
            attempts = 0
            while attempts < 3:
                try:
                    download_protein_biopython(code, pdbl, download_dir)
                    break  # Exit the retry loop if successful
                except Exception as e:
                    attempts += 1
                    logging.warning(f"Attempt {attempts} failed for protein {code}. Error: {e}")
                    if attempts < 3:
                        time.sleep(10)  # Wait 10 seconds before retrying
                    else:
                        logging.error(f"Failed to download protein {code} after 3 attempts.")

    except FileNotFoundError:
        logging.error(f"Error: File '{protein_list_file}' not found.")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        raise

def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    cwd = os.getcwd()
    protein_file_path = os.path.join(cwd, "protein.txt")
    
    # Define and create the "protein" directory
    protein_dir = os.path.join(cwd, "protein")
    os.makedirs(protein_dir, exist_ok=True)
    
    logging.debug(f"Current working directory: {cwd}")
    logging.debug(f"Protein file path: {protein_file_path}")

    try:
        download_proteins(protein_file_path, protein_dir)
    except Exception:
        pass

if __name__ == "__main__":
    main()