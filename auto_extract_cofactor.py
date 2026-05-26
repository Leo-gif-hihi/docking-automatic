import os
import urllib.request
import urllib.parse
import re
import logging
import gzip
import shutil
import time

def download_chebi_accession_file(filepath="database_accession.tsv"):
    """
    Downloads the compressed database_accession.tsv.gz file from ChEBI's HTTPS server,
    extracts it, and saves it as a standard TSV.
    
    Args:
        filepath (str): The local path where the extracted file should be saved.
    """
    time.sleep(1)  # Sleep to avoid overwhelming the server with requests
    # Updated to the new HTTPS path and .gz filename
    url = "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/database_accession.tsv.gz"
    gz_filepath = filepath + ".gz"

    if not os.path.exists(filepath):
        for attempt in range(3):
            try:
                # 1. Download the compressed file
                logging.info(f"Downloading {gz_filepath} from {url}... (Attempt {attempt+1}/3)")
                urllib.request.urlretrieve(url, gz_filepath)
                logging.info(f"Successfully downloaded {gz_filepath}. Extracting...")
                
                # 2. Extract the .gz file to the target filepath
                with gzip.open(gz_filepath, 'rb') as f_in:
                    with open(filepath, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        
                logging.info(f"Successfully extracted to {filepath}.")
                
                # 3. Clean up the compressed .gz file to save disk space
                os.remove(gz_filepath)
                break  # Exit the retry loop on success
                
            except Exception as e:
                logging.error(f"Failed to download or extract {url}: {e}")
                if attempt < 2:
                    time.sleep(10)  # Wait before retrying
    else:
        logging.info(f"File {filepath} already exists. Skipping download.")

def get_uniprot_cofactor(uniprot_id):
    """
    Extracts cofactor information for a target UniProt ID using the UniProt REST API
    and returns a list of unique ChEBI IDs.
    
    Args:
        uniprot_id (str): The UniProt Accession ID (e.g., 'V5NC32').
        
    Returns:
        list: A list of unique ChEBI IDs (e.g., ['CHEBI:18420', 'CHEBI:49883']).
    """
    time.sleep(1)
    url = f"https://rest.uniprot.org/uniprotkb/search?query=accession_id:{urllib.parse.quote(uniprot_id)}&format=tsv&fields=accession,protein_name,cc_cofactor"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req) as response:
                data = response.read().decode('utf-8')
                
            lines = data.strip().split('\n')
            if len(lines) > 1:
                # The first line is the header, the second line has the data
                columns = lines[1].split('\t')
                # cc_cofactor is the 3rd field requested, so index 2
                if len(columns) >= 3:
                    cc_cofactor = columns[2].strip()
                    chebi_ids = re.findall(r'ChEBI:(CHEBI:\d+)', cc_cofactor)
                    return list(set(chebi_ids))
            return []  # Success but no cofactors
        except Exception as e:
            logging.error(f"Failed to fetch cofactor info for {uniprot_id} (Attempt {attempt+1}/3): {e}")
            if attempt < 2:
                time.sleep(5)
        
    return []

def map_chebi_to_pdb(chebi_ids, filepath="database_accession.tsv"):
    """
    Maps a list of ChEBI IDs to their corresponding PDB chemical component IDs.
    
    Args:
        chebi_ids (list): List of ChEBI IDs (e.g., ['CHEBI:18420']).
        filepath (str): Path to the database_accession.tsv file.
        
    Returns:
        dict: A dictionary mapping ChEBI IDs to lists of mapped PDB IDs.
    """
    if not os.path.exists(filepath):
        logging.error(f"Mapping file not found: {filepath}")
        return {}
        
    # Create a mapping from numeric ChEBI ID to full ChEBI ID (e.g., '18420' -> 'CHEBI:18420')
    numeric_to_chebi = {}
    for cid in chebi_ids:
        match = re.search(r'\d+', cid)
        if match:
            numeric_to_chebi[match.group()] = cid
            
    pdb_mapping = {cid: set() for cid in chebi_ids}
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                columns = line.strip().split('\t')
                if len(columns) >= 6:
                    compound_id = columns[1]
                    accession_number = columns[2]
                    source_id = columns[5]
                    
                    # source_id '63' represents PDBeChem (PDB IDs)
                    if source_id == '63' and compound_id in numeric_to_chebi:
                        chebi_key = numeric_to_chebi[compound_id]
                        pdb_mapping[chebi_key].add(accession_number)
    except Exception as e:
        logging.error(f"Error reading {filepath}: {e}")
        
    return {k: list(v) for k, v in pdb_mapping.items() if v}

def get_pdb_cofactors_for_uniprot(uniprot_id, filepath="database_accession.tsv"):
    """
    Orchestrates the overall process: downloads the database mapping file (if needed),
    extracts ChEBI cofactors from UniProt, and maps them to PDB IDs.
    
    Args:
        uniprot_id (str): The UniProt Accession ID.
        filepath (str): Path to the database mapping file.
        
    Returns:
        dict: The resulting mapping of ChEBI IDs to PDB IDs.
    """
    logging.debug(f"Starting execution for UniProt ID: {uniprot_id}")
    
    # 1. Download
    download_chebi_accession_file(filepath)
    
    # 2. Get ChEBI IDs
    chebi_ids = get_uniprot_cofactor(uniprot_id)
    if not chebi_ids:
        logging.debug(f"No cofactors found for {uniprot_id}.")
        return {}
        
    logging.debug(f"Found ChEBI cofactors: {chebi_ids}")
    
    # 3. Map to PDB
    pdb_mapping = map_chebi_to_pdb(chebi_ids, filepath)
    logging.debug(f"Mapped PDB cofactors: {pdb_mapping}")
    
    return pdb_mapping
