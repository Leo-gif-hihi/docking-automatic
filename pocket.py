import time
import os
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio import Align
import requests
import json
import numpy as np
from sklearn.cluster import DBSCAN

def extract_uniprot_ids_from_pdb(pdb_file):
    """
    Parses a PDB file to extract UniProt IDs for each chain from DBREF records.
    Returns a dictionary mapping chain ID to UniProt ID.
    """
    chain_to_uniprot = {}
    
    with open(pdb_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("DBREF"):
                # Example DBREF line
                # DBREF  2X9A A    1   373  UNP    P00519   ABL1_HUMAN       9   381
                parts = line.split()
                if len(parts) >= 7 and parts[5] == "UNP":
                    chain_id = parts[2]
                    uniprot_id = parts[6]
                    chain_to_uniprot[chain_id] = uniprot_id
                    
    return chain_to_uniprot

def fetch_biolip_data(uniprot_id, max_retries=3):
    """
    Fetches empirical binding data from the BioLiP API for a given UniProt ID.
    Implements exponential backoff to handle rate limits in the form of HTML responses.
    Caches responses locally to avoid redundant API calls.
    """
    cache_dir = "biolip_cache"
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, f"{uniprot_id}.txt")
    
    if os.path.exists(cache_file):
        print(f"Loading BioLiP data from cache for UniProt ID: {uniprot_id}")
        with open(cache_file, 'r', encoding='utf-8') as f:
            return f.read()

    url = f"https://aideepmed.com/BioLiP/qsearch.cgi?uniprot={uniprot_id}&outfmt=txt"
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
    }

    for attempt in range(max_retries):
        try:
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            
            text = response.text
            #print(text)  # Debug print to check the response content
            # Check for the HTML rate limit trap or explicitly "Too many requests"
            if "Too many requests" in text or text.strip().startswith("<html") or "<html>" in text.lower():
                print(f"Rate limited by BioLiP for UniProt ID {uniprot_id}. Attempt {attempt + 1}/{max_retries}.")
            else:
                with open(cache_file, 'w', encoding='utf-8') as f:
                    f.write(text)
                return text
                
        except requests.exceptions.RequestException as e:
            print(f"Request error fetching BioLiP data for UniProt ID {uniprot_id}: {e}")
            
        if attempt < max_retries - 1:
            sleep_time = 30
            print(f"Sleeping for {sleep_time} seconds before retrying...")
            time.sleep(sleep_time)
            
    print(f"Failed to fetch BioLiP data for UniProt ID {uniprot_id} after {max_retries} attempts.")
    return None

def parse_biolip_tsv(tsv_text):
    """
    Parses BioLiP TSV data to extract binding residues (Index 8) and protein sequence (Index 20).
    Returns a list of dictionaries with parsed data.
    """
    parsed_entries = []
    lines = tsv_text.strip().split('\n')
    for line in lines:
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) > 20:
            # Extract binding residues and split by space into a list
            binding_residues_str = parts[8].strip()
            binding_residues = binding_residues_str.split() if binding_residues_str else []
            
            # Extract protein sequence
            protein_sequence = parts[20].strip()
            
            parsed_entries.append({
                'binding_residues': binding_residues,
                'sequence': protein_sequence
            })
            
    return parsed_entries

def extract_sequence_from_pdb_atoms(pdb_file):
    """
    Extracts the amino acid sequence for each chain physically present in the PDB file 
    based solely on ATOM records (avoiding missing loops from SEQRES).
    Returns a dict mapping chain ID to a list of tuples: (1-letter-AA, ResSeq).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    chain_sequences = {}

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            seq_data = []
            for residue in chain:
                # Filter out hetero/water residues; pure amino acids have " " as hetero flag
                if residue.get_id()[0] != " ":
                    continue
                # Ensure the residue has a CA atom (truly part of the backbone)
                if 'CA' in residue:
                    resname = residue.get_resname().strip().upper()
                    resseq = residue.get_id()[1]
                    # Map 3-letter to 1-letter code; use 'X' for unknown
                    aa_1_letter = protein_letters_3to1.get(resname, 'X')
                    seq_data.append((aa_1_letter, resseq))
            
            if seq_data:
                chain_sequences[chain_id] = seq_data
                
        # Only process the first model
        break
        
    return chain_sequences

def align_sequences(biolip_sequence, pdb_atom_sequence):
    """
    Performs a Needleman-Wunsch global alignment between the BioLiP sequence 
    and the PDB ATOM sequence using sensible conservative penalties.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(biolip_sequence, pdb_atom_sequence)
    if not alignments:
        return None
        
    best_alignment = alignments[0]
    return best_alignment

def calculate_sequence_identity(best_alignment):
    """
    Calculates the sequence identity (matches / total alignment length).
    Returns the percentage.
    """
    seq1, seq2 = best_alignment[0], best_alignment[1]
    
    matches = 0
    total_length = len(seq1)
    
    for char1, char2 in zip(seq1, seq2):
        if char1 == char2 and char1 != '-':
            matches += 1
                
    if total_length == 0:
        return 0.0
        
    return (matches / total_length) * 100.0

def map_binding_residues(best_alignment, binding_residues, pdb_seq_data):
    """
    Maps BioLiP interacting residues to PDB resseq numbers using the global alignment.
    Returns a list of valid mapped resseq numbers in the PDB file.
    """
    # In Biopython PairwiseAligner, the alignment object behaves like a tuple of strings (with gap characters)
    biolip_aligned = best_alignment[0]
    pdb_aligned = best_alignment[1]
    
    biolip_to_pdb_map = {}
    biolip_pos = 0
    pdb_pos = 0
    
    # Traverse the alignment to build the index cross-reference map
    for char_b, char_p in zip(biolip_aligned, pdb_aligned):
        if char_b != '-':
            if char_p != '-':
                biolip_to_pdb_map[biolip_pos] = pdb_pos
            biolip_pos += 1
        if char_p != '-':
            pdb_pos += 1

    mapped_resseqs = []
    
    for res_str in binding_residues:
        if not res_str or len(res_str) < 2:
            continue
            
        req_aa = res_str[0] # Amino acid 1-letter code
        try:
            # BioLiP uses 1-based indexing, convert to 0-based
            req_idx = int(res_str[1:]) - 1
        except ValueError:
            print(f"Skipping malformed BioLiP residue string: {res_str}")
            continue
            
        if req_idx in biolip_to_pdb_map:
            mapped_pdb_idx = biolip_to_pdb_map[req_idx]
            pdb_aa, resseq = pdb_seq_data[mapped_pdb_idx]
            
            # Validation check for mutation between pure UniProt sequence and PDB crystal homologous chain
            if pdb_aa != req_aa:
                print(f"Mutation warning at BioLiP {res_str}: PDB structure has {pdb_aa}. Keeping coordinate.")
            
            mapped_resseqs.append(resseq)
        else:
            print(f"Missing loop warning: BioLiP residue {res_str} maps to a gap in the PDB structure. Dropping.")
            
    return mapped_resseqs

def extract_active_coordinates(pdb_file, active_residues):
    """
    Extracts 3D coordinates (x,y,z) for the Alpha-Carbon of high-scoring residues.
    Returns coords (Nx3 NumPy array), scores (Nx1 NumPy array), and corresponding residue mappings.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    coords = []
    scores = []
    residue_ids = []
    
    for model in structure:
        for chain in model:
            cid = chain.get_id()
            for res in chain:
                resseq = res.get_id()[1]
                if (cid, resseq) in active_residues:
                    if 'CA' in res:
                        # Extract the Alpha-Carbon 3D coordinate vector
                        coord = res['CA'].get_coord()
                        coords.append(coord)
                        scores.append(active_residues[(cid, resseq)])
                        residue_ids.append((cid, resseq))
        break # Only process the first model
        
    return np.array(coords), np.array(scores), residue_ids

def cluster_and_select_pocket(coords, scores, residue_ids, eps=10.0, min_samples=3):
    """
    Clusters 3D coordinates using DBSCAN to identify distinct binding pockets.
    Resolves multiple clusters by summing BioLiP occurrence scores and taking the maximum.
    """
    if len(coords) == 0:
        return None
        
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    labels = dbscan.fit_predict(coords)
    
    unique_labels = set(labels)
    best_cluster = -1
    max_score = -1
    
    for label in unique_labels:
        if label == -1:
            # Skip noise points identified by DBSCAN
            continue
            
        # Sum the BioLiP occurrence scores of all residues within this cluster
        cluster_mask = (labels == label)
        cluster_score = np.sum(scores[cluster_mask])
        print(f"  -> DBSCAN Cluster {label}: score {cluster_score}, size {np.sum(cluster_mask)} atoms")
        
        # Select cluster with the highest occurrence score as the primary active site
        if cluster_score > max_score:
            max_score = cluster_score
            best_cluster = label
            
    if best_cluster == -1:
        print("  -> DBSCAN failed to find any dense clusters. Check eps/min_samples.")
        return None
        
    print(f"*** Primary Active Site Resolved: Cluster {best_cluster} with Total Score {max_score} ***")
    
    best_coords = coords[labels == best_cluster]
    best_residue_ids = [residue_ids[i] for i in range(len(labels)) if labels[i] == best_cluster]
    
    return best_coords, best_residue_ids

def calculate_bounding_box(coords, padding=5.0):
    """
    Calculates the spatial center and dimensions of the selected pocket in 3D space,
    adding the user-defined padding.
    """
    if len(coords) == 0:
        return None
        
    min_x, min_y, min_z = np.min(coords, axis=0)
    max_x, max_y, max_z = np.max(coords, axis=0)
    
    center_x = (max_x + min_x) / 2
    center_y = (max_y + min_y) / 2
    center_z = (max_z + min_z) / 2
    
    size_x = (max_x - min_x) + padding
    size_y = (max_y - min_y) + padding
    size_z = (max_z - min_z) + padding
    
    return {
        'center_x': center_x, 'center_y': center_y, 'center_z': center_z,
        'size_x': size_x, 'size_y': size_y, 'size_z': size_z
    }

def collect_interactive_data(chain_to_uniprot):
    """
    Deduplicates UniProt IDs and fetches BioLiP data for each unique ID.
    Returns mapping of UniProt ID to BioLiP data, and chains grouped by UniProt ID.
    """
    unique_uniprot_ids = set(chain_to_uniprot.values())
    uniprot_to_chains = {u: [] for u in unique_uniprot_ids}
    
    for chain, uniprot in chain_to_uniprot.items():
        uniprot_to_chains[uniprot].append(chain)
        
    biolip_data = {}
    for uniprot_id in unique_uniprot_ids:
        print(f"Fetching BioLiP data for UniProt ID: {uniprot_id}")
        
        # Check if it will hit cache before fetching to know if we need to sleep
        is_cached = os.path.exists(os.path.join("biolip_cache", f"{uniprot_id}.txt"))
        
        data = fetch_biolip_data(uniprot_id)
        if data:
            parsed_data = parse_biolip_tsv(data)
            biolip_data[uniprot_id] = parsed_data
            
        if not is_cached:
            sleep_time = 10
            print(f"Sleeping for {sleep_time} seconds before next API request...")
            time.sleep(sleep_time)
            
    return biolip_data, uniprot_to_chains

def calculate_volume_and_exhaustiveness(box_params):
    """
    Calculates grid box volume and determines Vina exhaustiveness scaling.
    """
    volume = box_params['size_x'] * box_params['size_y'] * box_params['size_z']
    exhaustiveness = 32 if volume <= 27000 else 64
    return volume, exhaustiveness

def write_vina_box_file(protein_file, box_path, box_params, exhaustiveness):
    """
    Writes the Vina configuration box file.
    """
    import os
    os.makedirs(box_path, exist_ok=True)
    box_file = box_path / f"{protein_file.stem}.box.txt"
    
    with open(box_file, 'w', encoding='utf-8') as f:
        f.write(f"center_x = {box_params['center_x']:.3f}\n")
        f.write(f"center_y = {box_params['center_y']:.3f}\n")
        f.write(f"center_z = {box_params['center_z']:.3f}\n")
        f.write(f"size_x = {box_params['size_x']:.3f}\n")
        f.write(f"size_y = {box_params['size_y']:.3f}\n")
        f.write(f"size_z = {box_params['size_z']:.3f}\n")
        f.write(f"exhaustiveness = {exhaustiveness}\n")
        
    print(f"Wrote Vina box configuration to {box_file}")

def generate_pymol_box_script(protein_file, box_params, output_dir="output"):
    """
    Generates a PyMOL script (.pml) to visualize the generated heatmap 
    PDB alongside the calculated Vina binding box.
    """
    import os
    from pathlib import Path
    out_path = Path(output_dir)
    out_path.mkdir(exist_ok=True, parents=True)
    
    pml_file = out_path / f"{protein_file.stem}_visualize.pml"
    heatmap_file = f"{protein_file.stem}_heatmap.pdb"
    
    # Calculate box corners based on center and size
    min_x = box_params['center_x'] - (box_params['size_x'] / 2)
    max_x = box_params['center_x'] + (box_params['size_x'] / 2)
    min_y = box_params['center_y'] - (box_params['size_y'] / 2)
    max_y = box_params['center_y'] + (box_params['size_y'] / 2)
    min_z = box_params['center_z'] - (box_params['size_z'] / 2)
    max_z = box_params['center_z'] + (box_params['size_z'] / 2)

    script_content = f"""
# Load the heatmap structure
load {heatmap_file}
hide all
show cartoon

# Color the protein based on the mapped occurrence scores stored in B-factor
spectrum b, white_red, minimum=0, maximum=100

# Script to draw the 3D Vina box using CGO
python
from pymol.cgo import *
from pymol import cmd

box = [
    BEGIN, LINES,
    COLOR, 0.0, 1.0, 0.0, # Green Box
    
    # Bottom Face
    VERTEX, {min_x}, {min_y}, {min_z},
    VERTEX, {max_x}, {min_y}, {min_z},
    VERTEX, {min_x}, {max_y}, {min_z},
    VERTEX, {max_x}, {max_y}, {min_z},
    VERTEX, {min_x}, {min_y}, {min_z},
    VERTEX, {min_x}, {max_y}, {min_z},
    VERTEX, {max_x}, {min_y}, {min_z},
    VERTEX, {max_x}, {max_y}, {min_z},

    # Top Face
    VERTEX, {min_x}, {min_y}, {max_z},
    VERTEX, {max_x}, {min_y}, {max_z},
    VERTEX, {min_x}, {max_y}, {max_z},
    VERTEX, {max_x}, {max_y}, {max_z},
    VERTEX, {min_x}, {min_y}, {max_z},
    VERTEX, {min_x}, {max_y}, {max_z},
    VERTEX, {max_x}, {min_y}, {max_z},
    VERTEX, {max_x}, {max_y}, {max_z},

    # Vertical Pillars
    VERTEX, {min_x}, {min_y}, {min_z},
    VERTEX, {min_x}, {min_y}, {max_z},
    VERTEX, {max_x}, {min_y}, {min_z},
    VERTEX, {max_x}, {min_y}, {max_z},
    VERTEX, {min_x}, {max_y}, {min_z},
    VERTEX, {min_x}, {max_y}, {max_z},
    VERTEX, {max_x}, {max_y}, {min_z},
    VERTEX, {max_x}, {max_y}, {max_z},
    END
]

cmd.load_cgo(box, "docking_box")
python end

center docking_box
zoom docking_box, 15
"""
    with open(pml_file, 'w', encoding='utf-8') as f:
        f.write(script_content.strip())
    print(f"Generated PyMOL Visualization Script: {pml_file}")

def generate_heatmap_pdb(protein_file, active_residues, output_dir="output"):
    """
    Step 12: Writes BioLiP occurrence scores into the PDB B-factor column.
    Allows users to visually verify the pipeline's decisions in 3D (e.g., in PyMOL).
    """
    import os
    from pathlib import Path
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_file)
    
    max_score = max(active_residues.values()) if active_residues else 0
    
    for model in structure:
        for chain in model:
            cid = chain.get_id()
            for res in chain:
                resseq = res.get_id()[1]
                count = active_residues.get((cid, resseq), 0)
                # Normalize the occurrence score to 0-100 for the B-factor column
                norm_score = (count / max_score) * 100.0 if max_score > 0 else 0.0
                
                for atom in res:
                    atom.set_bfactor(norm_score)
                    
    out_path = Path(output_dir)
    out_path.mkdir(exist_ok=True, parents=True)
    
    io = PDBIO()
    io.set_structure(structure)
    heatmap_file = out_path / f"{protein_file.stem}_heatmap.pdb"
    io.save(str(heatmap_file))
    print(f"Generated Heatmap PDB for visual verification: {heatmap_file}")

def process_pockets(protein_path, box_path, output_dir="output"):
    """Phase 1: Identify pockets by querying BioLiP data for UniProt IDs extracted from PDB."""
    protein_files = list(protein_path.glob("*.pdb"))
    if not protein_files:
        print(f"No .pdb files found in {protein_path} for pocket identification.")
        return

    for protein_file in protein_files:
        print(f"\n--- Processing {protein_file.name} for pocket identification ---")
        chain_to_uniprot = extract_uniprot_ids_from_pdb(protein_file)
        
        if not chain_to_uniprot:
            print(f"No UniProt mappings found in DBREF records for {protein_file.name}.")
            continue
            
        print(f"Extracted UniProt mappings: {chain_to_uniprot}")
        
        # Extract sequences physically present in the ATOM records to avoid missing loops
        chain_sequences = extract_sequence_from_pdb_atoms(protein_file)
        print(f"Extracted physical sequences for chains: {list(chain_sequences.keys())}")
        
        # Step 7: Initialize occurrence heatmap for all physical residues
        occurrence_heatmap = {}
        for cid, seq_data in chain_sequences.items():
            for _, resseq in seq_data:
                occurrence_heatmap[(cid, resseq)] = 0
        
        biolip_data, uniprot_to_chains = collect_interactive_data(chain_to_uniprot)
        
        print(f"Grouped chains by UniProt ID: {uniprot_to_chains}")
        
        # Phase 1, Step 5: Global Alignment
        alignment_cache = {}
        for uniprot_id, chains in uniprot_to_chains.items():
            if uniprot_id not in biolip_data or not biolip_data[uniprot_id]:
                continue
                
            entries = biolip_data[uniprot_id]
            for chain_id in chains:
                if chain_id not in chain_sequences:
                    continue
                    
                # The PDB ATOM sequence is kept as a list of tuples (AA, resseq)
                # We need to construct the pure AA string for alignment
                pdb_seq_data = chain_sequences[chain_id]
                pdb_atom_sequence_str = "".join([aa for aa, resseq in pdb_seq_data])
                
                # We might have multiple BioLiP entries per UniProt ID (different ligands)
                for i, entry in enumerate(entries):
                    biolip_sequence = entry['sequence']
                    
                    cache_key = (biolip_sequence, pdb_atom_sequence_str)
                    
                    if cache_key in alignment_cache:
                        print(f"Aligning Chain {chain_id} (UniProt {uniprot_id}) to BioLiP Record {i+1} (using cached alignment)...")
                        best_alignment = alignment_cache[cache_key]
                    else:
                        print(f"Aligning Chain {chain_id} (UniProt {uniprot_id}) to BioLiP Record {i+1}...")
                        best_alignment = align_sequences(biolip_sequence, pdb_atom_sequence_str)
                        alignment_cache[cache_key] = best_alignment
                    
                    if best_alignment:
                        print(f"Alignment Score: {best_alignment.score}")
                        identity_pct = calculate_sequence_identity(best_alignment)
                        print(f"Sequence Identity: {identity_pct:.2f}%")
                        
                        if identity_pct <= 80.0:
                            print(f"Sequence identity {identity_pct:.2f}% is not above 80%, skipping...")
                            continue

                        mapped_resseqs = map_binding_residues(best_alignment, entry['binding_residues'], pdb_seq_data)
                        print(f"Successfully mapped {len(mapped_resseqs)} valid 3D coordinates for BioLiP Record {i+1}: {mapped_resseqs}")
                        
                        # Step 7: Aggregation
                        for resseq in mapped_resseqs:
                            occurrence_heatmap[(chain_id, resseq)] += 1

        # Step 7: Merge the final score arrays into a single, unified pool
        active_residues = {k: count for k, count in occurrence_heatmap.items() if count > 0}
        print(f"Aggregated {len(active_residues)} active residues across all chains for {protein_file.name}")
        
        # Step 12: Generate Heatmap PDB
        if active_residues:
            generate_heatmap_pdb(protein_file, active_residues, output_dir=output_dir)
        
        # Step 8: 3D Spatial Clustering (DBSCAN)
        if active_residues:
            coords, scores, residue_ids = extract_active_coordinates(protein_file, active_residues)
            print(f"Executing DBSCAN Clustering on {len(coords)} accumulated 3D spatial points...")
            primary_pocket = cluster_and_select_pocket(coords, scores, residue_ids)
            
            if primary_pocket:
                best_coords, best_residues = primary_pocket
                print(f"Final Primary Binding Pocket consists of {len(best_residues)} residues: {best_residues}")
                
                # Step 9: Bounding Box Calculation
                box_params = calculate_bounding_box(best_coords, padding=5.0)
                print(f"Calculated Vina Grid Box: Center({box_params['center_x']:.2f}, {box_params['center_y']:.2f}, {box_params['center_z']:.2f}) | Dimensions({box_params['size_x']:.2f}, {box_params['size_y']:.2f}, {box_params['size_z']:.2f})")
                
                # Step 10: Volume Check & Dynamic Exhaustiveness
                volume, exhaustiveness = calculate_volume_and_exhaustiveness(box_params)
                print(f"Box Volume: {volume:.2f} Å³ -> Scaled Exhaustiveness: {exhaustiveness}")
                write_vina_box_file(protein_file, box_path, box_params, exhaustiveness)
                
                # NEW CALL: Generate Visual Script
                generate_pymol_box_script(protein_file, box_params, output_dir=output_dir)
                
            else:
                print("Could not resolve a primary binding pocket via clustering.")
        else:
            print("No active residues found to cluster.")
        # At this point, active_residues contains the pooled binding residues ready for 3D clustering

