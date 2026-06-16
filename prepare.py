import os
import glob
import argparse
import sys
import logging
from prody import parseMMCIF, writeMMCIF, confProDy
from logger_utils import log_step, console

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
            writeMMCIF(out_filepath, clean_selection)
            logging.debug(f" -> Saved: {out_filepath}")
            return out_filepath
    else:
        logging.debug(f" -> Skipped {filename}: No atoms left after cleaning.")
        return None

def get_suggested_cofactors(filepath):
    """Helper to query UniProt mapping logic and return suggested cofactors."""
    from pocket import extract_uniprot_ids_from_cif
    from auto_extract_cofactor import get_pdb_cofactors_for_uniprot
    
    chain_to_uniprot = extract_uniprot_ids_from_cif(filepath)
    essential_cofactors = set()
    uniprot_to_pdb_mapping = {}
    for chain_id, uniprot_id in chain_to_uniprot.items():
        if uniprot_id not in uniprot_to_pdb_mapping:
            uniprot_to_pdb_mapping[uniprot_id] = get_pdb_cofactors_for_uniprot(uniprot_id)
        
        pdb_mapping = uniprot_to_pdb_mapping[uniprot_id]
        for chebi, pdb_list in pdb_mapping.items():
            for pdb_id in pdb_list:
                essential_cofactors.add(pdb_id.upper())
    return essential_cofactors

def generate_file_stats(filepath, file_hetatms, to_eliminate, suggested=None):
    """Helper to generate stats dictionary for a processed file."""
    filename = os.path.basename(filepath)
    if suggested is None:
        suggested = get_suggested_cofactors(filepath)
        
    removed = [h for h in file_hetatms if h in to_eliminate]
    kept = [h for h in file_hetatms if h not in to_eliminate]
    
    return {
        "Protein": filename,
        "Kept Ligands": "; ".join(kept),
        "Removed Ligands": "; ".join(removed),
        "Suggested Cofactors": "; ".join(suggested)
    }

def clean_global_mode(pdb_files, output_dir):
    """Runs the cleaning process in global mode (one prompt for all files)."""
    logging.info(f"Scanning {len(pdb_files)} files for HETATM residues...")
    all_hetatms = set()
    structures = []
    
    # Pre-parse and collect all HETATMs
    for filepath in pdb_files:
        structure = parseMMCIF(filepath)
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
    
    cleaned_paths = []
    file_stats = {}
    for filepath, structure in structures:
        cleaned_pdb_path = clean_and_save_pdb(structure, filepath, output_dir, sel_str)
        if cleaned_pdb_path:
            cleaned_paths.append(cleaned_pdb_path)
            
            filename = os.path.basename(filepath)
            file_hetatms = get_hetatms(structure)
            file_stats[filename] = generate_file_stats(filepath, file_hetatms, to_eliminate)
            
    return cleaned_paths, file_stats

def clean_local_mode(pdb_files, output_dir):
    """Runs the cleaning process in local mode (prompt per file)."""
    logging.info(f"Running in LOCAL mode. Scanning {len(pdb_files)} files individually.")
    logging.info("Tip: You can press Ctrl+C at any prompt to escape and stop processing.\n")
    
    cleaned_paths = []
    file_stats = {}
    for filepath in pdb_files:
        filename = os.path.basename(filepath)
        structure = parseMMCIF(filepath)
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
        cleaned_pdb_path = clean_and_save_pdb(structure, filepath, output_dir, sel_str)
        if cleaned_pdb_path:
            cleaned_paths.append(cleaned_pdb_path)
            
            file_stats[filename] = generate_file_stats(filepath, file_hetatms, to_eliminate)
        logging.debug("")  # Blank line for readability between files if in debug mode
        
    return cleaned_paths, file_stats

def clean_auto_mode(pdb_files, output_dir):
    """Runs the cleaning process in auto mode (keeps essential cofactors automatically)."""
    logging.debug(f"Running in AUTO mode. Scanning {len(pdb_files)} files individually.")
    
    cleaned_paths = []
    file_stats = {}
    for filepath in pdb_files:
        filename = os.path.basename(filepath)
        structure = parseMMCIF(filepath)
        if structure is None:
            continue
            
        file_hetatms = get_hetatms(structure)
        logging.debug(f"--- {filename} ---")
        
        suggested = set()
        to_eliminate = []
        if not file_hetatms:
            logging.debug("No HETATM residues found. Just cleaning water.")
            sel_str = get_selection_string(to_eliminate, file_hetatms)
        else:
            suggested = get_suggested_cofactors(filepath)
            
            for hetatm in file_hetatms:
                if hetatm.upper() not in suggested:
                    to_eliminate.append(hetatm)
            
            kept = suggested & set([h.upper() for h in file_hetatms])
            kept_str = ", ".join(kept) if kept else "none"
            action_desc = "ALL HETATM residues in this file" if set(to_eliminate) == file_hetatms else ", ".join(to_eliminate)
            logging.debug(f"Action: Eliminating {action_desc} (Kept {kept_str} as essential cofactors).")
            
            sel_str = get_selection_string(to_eliminate, file_hetatms)
            
        cleaned_pdb_path = clean_and_save_pdb(structure, filepath, output_dir, sel_str)
        if cleaned_pdb_path:
            cleaned_paths.append(cleaned_pdb_path)
            
            file_stats[filename] = generate_file_stats(filepath, file_hetatms, to_eliminate, suggested)
        logging.debug("")
        
    return cleaned_paths, file_stats

def run_pdbfixer(protein_path, fixed_protein_path):
    """Uses PDBFixer to rebuild missing heavy atoms in the protein structure."""
    import logging
    from pathlib import Path
    
    if Path(fixed_protein_path).exists():
        logging.debug(f"Skipping PDBFixer: {fixed_protein_path} already exists.")
        return True
        
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBxFile
    except ImportError:
        logging.error("PDBFixer/OpenMM is required to rebuild missing atoms but is not installed.")
        return False

    try:
        logging.info(f"Running PDBFixer on {protein_path} to rebuild missing atoms...")
        fixer = PDBFixer(filename=str(protein_path))
        
        # Find missing residues and atoms, then add them
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        
        # Save the repaired structure
        with open(fixed_protein_path, 'w') as f:
            PDBxFile.writeFile(fixer.topology, fixer.positions, f)
            
        logging.debug(f"Saved fixed structure to {fixed_protein_path}")
        return True
    except Exception as e:
        logging.error(f"PDBFixer failed for {protein_path}: {e}")
        return False

def run_reduce2(protein_path, protein_protonated):
    """Runs mmtbx.reduce2 to add hydrogens and optimize the structure."""
    import subprocess
    from pathlib import Path

    if Path(protein_protonated).exists():
        logging.debug(f"Skipping REDUCE2: {protein_protonated} already exists.")
        return True

    cmd_reduce2 = [
        "mmtbx.reduce2", str(protein_path), "approach=add", 
        "add_flip_movers=True", f"output.filename={protein_protonated}"
    ]
    try:
        subprocess.run(cmd_reduce2, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"REDUCE2 failed for {protein_path}:\n{e.stderr}")
        return False


def run_openmm_minimization(protein_protonated, freeze_backbone=True):
    """Minimizes the structure using OpenMM, optionally freezing the backbone."""
    try:
        from openmm.app import PDBxFile, ForceField, Simulation, HBonds, NoCutoff
        from openmm import LangevinMiddleIntegrator
        from openmm.unit import kelvin, picosecond, kilojoules_per_mole, nanometer
    except ImportError:
        logging.error("OpenMM is required to minimize structures but is not installed.")
        return False

    try:
        logging.info(f"Running OpenMM minimization on {protein_protonated}")
        pdb = PDBxFile(str(protein_protonated))
        
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml', 'implicit/gbn2.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None)
        
        # --- Freeze the backbone to preserve the crystal structure ---
        if freeze_backbone:
            logging.debug("Freezing backbone atoms (N, CA, C, O) to prevent structural drift...")
            for atom in pdb.topology.atoms():
                if atom.name in ['N', 'CA', 'C', 'O']: # Identify backbone atoms
                    system.setParticleMass(atom.index, 0.0) # Mass of 0 freezes the atom
        # -----------------------------------------------------------------

        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picosecond)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        
        logging.debug("Minimizing energy...")
        simulation.minimizeEnergy()
        
        positions = simulation.context.getState(getPositions=True).getPositions()
        
        # Overwrite the protein_protonated file with the minimized structure
        with open(protein_protonated, 'w') as f:
            PDBxFile.writeFile(simulation.topology, positions, f)
            
        logging.debug(f"Saved minimized structure to {protein_protonated}")
        return True
    except Exception as e:
        logging.error(f"OpenMM minimization failed: {e}")
        return False

def fix_cif_indentation(cif_path):
    """ProDy (used by Meeko) fails to parse mmCIF files if ATOM/HETATM lines have leading spaces.
       This function strips leading spaces from these lines."""
    with open(cif_path, 'r') as f:
        lines = f.readlines()
    with open(cif_path, 'w') as f:
        for line in lines:
            if line.lstrip().startswith('ATOM') or line.lstrip().startswith('HETATM'):
                f.write(line.lstrip())
            else:
                f.write(line)

def run_meeko_receptor(protein_protonated, protein_prep_out, protein_pdbqt):
    """Runs mk_prepare_receptor.py (Meeko) to prepare the receptor."""
    import subprocess
    from pathlib import Path

    if Path(protein_pdbqt).exists():
        logging.debug(f"Skipping Meeko: {protein_pdbqt} already exists.")
        return True
        
    fix_cif_indentation(protein_protonated)

    cmd_meeko = [
        "mk_prepare_receptor.py", "-i", str(protein_protonated), 
        "-o", str(protein_prep_out), "-p", "-a"
    ]
    try:
        subprocess.run(cmd_meeko, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Meeko failed to prepare {protein_protonated}:\n{e.stderr}")
        return False

def _extract_cif_loop(lines, loop_prefix):
    """Extracts a loop block from CIF lines based on a prefix (e.g., '_entity.')."""
    loop_lines = []
    in_loop = False
    in_kv = False
    in_multiline = False
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        if in_multiline:
            loop_lines.append(line.rstrip() + '\n')
            if line.startswith(';'):
                in_multiline = False
            continue
            
        if stripped == "loop_":
            in_kv = False
            j = i + 1
            is_target = False
            while j < len(lines) and lines[j].strip().startswith("_"):
                if lines[j].strip().startswith(loop_prefix):
                    is_target = True
                    break
                j += 1
            
            if is_target:
                in_loop = True
                loop_lines.append(line.rstrip() + '\n')
                continue
                
        if in_loop:
            if stripped == "loop_" or (stripped.startswith("_") and not stripped.startswith(loop_prefix)):
                break
            if line.startswith(';'):
                in_multiline = True
            if stripped != "#":
                loop_lines.append(line.rstrip() + '\n')
        else:
            if stripped.startswith(loop_prefix):
                in_kv = True
                loop_lines.append(line.rstrip() + '\n')
            elif in_kv:
                if stripped.startswith("_") or stripped == "loop_" or (stripped == "#" and not line.startswith(";")):
                    in_kv = False
                else:
                    if line.startswith(';'):
                        in_multiline = True
                    loop_lines.append(line.rstrip() + '\n')
                
    return loop_lines

def _get_atom_site_cols(lines):
    """Returns the start index, end index, and column names of the _atom_site loop."""
    start_idx, end_idx = -1, -1
    col_names = []
    in_header, in_data = False, False
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "loop_":
            if i + 1 < len(lines) and lines[i+1].strip().startswith("_atom_site."):
                in_header = True
                start_idx = i
                continue
                
        if in_header:
            if stripped.startswith("_atom_site."):
                col_names.append(stripped)
            elif stripped.startswith("ATOM") or stripped.startswith("HETATM"):
                in_header = False
                in_data = True
                
        if in_data:
            if stripped == "loop_" or stripped == "#" or (stripped.startswith("_") and not stripped.startswith("_atom_site.")):
                end_idx = i
                break
                
    if in_data and end_idx == -1:
        end_idx = len(lines)
        
    return start_idx, end_idx, col_names

def _get_asym_to_entity_from_atom_site(lines, col_names, start_idx, end_idx):
    """Extracts asym_to_entity mapping from an _atom_site loop."""
    mapping = {}
    try:
        asym_idx = col_names.index("_atom_site.label_asym_id")
        entity_idx = col_names.index("_atom_site.label_entity_id")
    except ValueError:
        return mapping
        
    data_start = start_idx + 1 + len(col_names)
    for i in range(data_start, end_idx):
        tokens = lines[i].strip().split()
        if len(tokens) > max(asym_idx, entity_idx):
            asym = tokens[asym_idx]
            ent = tokens[entity_idx]
            if asym not in mapping:
                mapping[asym] = ent
    return mapping

def _infer_asym_to_entity(proc_lines, col_names, start_idx, end_idx):
    """Infers asym_to_entity mapping based on comp_id (HOH/WAT vs others)."""
    try:
        asym_idx = col_names.index("_atom_site.label_asym_id")
        comp_idx = col_names.index("_atom_site.label_comp_id")
        seq_idx = col_names.index("_atom_site.label_seq_id")
    except ValueError:
        return {}, []
        
    chain_info = {}
    data_start = start_idx + 1 + len(col_names)
    for i in range(data_start, end_idx):
        line = proc_lines[i].strip()
        if not line or line.startswith("#"): continue
        tokens = line.split()
        if len(tokens) <= max(asym_idx, comp_idx, seq_idx): continue
            
        asym = tokens[asym_idx]
        seq = tokens[seq_idx]
        comp = tokens[comp_idx]
        
        if asym not in chain_info:
            chain_info[asym] = {'has_seq': False, 'is_water': False}
        if seq not in ('.', '?'):
            chain_info[asym]['has_seq'] = True
        if comp in ('HOH', 'WAT'):
            chain_info[asym]['is_water'] = True
            
    asym_to_entity = {}
    entity_lines = ["loop_\n", "_entity.id\n", "_entity.type\n"]
    current_ent_id = 1
    
    for asym, info in chain_info.items():
        ent_type = 'polymer' if info['has_seq'] else ('water' if info['is_water'] else 'non-polymer')
        asym_to_entity[asym] = str(current_ent_id)
        entity_lines.append(f"{current_ent_id} {ent_type}\n")
        current_ent_id += 1
        
    return asym_to_entity, entity_lines

def _update_struct_asym_loop(proc_lines, asym_to_entity):
    """Replaces or creates _struct_asym loop with entity_id."""
    start_idx, end_idx = -1, -1
    for i, line in enumerate(proc_lines):
        if line.strip() == "loop_" and i + 1 < len(proc_lines) and proc_lines[i+1].strip().startswith("_struct_asym."):
            start_idx = i
            j = i + 1
            while j < len(proc_lines):
                s = proc_lines[j].strip()
                if s == "#" or s == "loop_" or (s.startswith("_") and not s.startswith("_struct_asym.")):
                    end_idx = j
                    break
                j += 1
            if end_idx == -1: end_idx = len(proc_lines)
            break
            
    if start_idx != -1:
        new_loop = ["loop_\n", "  _struct_asym.id\n", "  _struct_asym.entity_id\n"]
        for asym, ent in asym_to_entity.items():
            new_loop.append(f"   {asym} {ent}\n")
        proc_lines[start_idx:end_idx] = new_loop

def _update_atom_site_entity_id(proc_lines, col_names, start_idx, end_idx, asym_to_entity):
    """Updates _atom_site.label_entity_id in proc_lines."""
    try:
        asym_idx = col_names.index("_atom_site.label_asym_id")
        entity_idx = col_names.index("_atom_site.label_entity_id")
    except ValueError:
        return # Cannot update without these columns
        
    data_start = start_idx + 1 + len(col_names)
    for i in range(data_start, end_idx):
        line = proc_lines[i].strip()
        if not line or line.startswith("#"): continue
        tokens = line.split()
        if len(tokens) > max(asym_idx, entity_idx):
            asym = tokens[asym_idx]
            if asym in asym_to_entity:
                tokens[entity_idx] = asym_to_entity[asym]
                proc_lines[i] = " ".join(tokens) + "\n"

def _get_asym_mappings_from_atom_site(lines, col_names, start_idx, end_idx):
    """Extracts mapping of label_asym_id to auth_asym_id from the _atom_site loop."""
    mapping = {}
    try:
        label_idx = col_names.index("_atom_site.label_asym_id")
        auth_idx = col_names.index("_atom_site.auth_asym_id")
    except ValueError:
        return mapping
        
    comp_idx = col_names.index("_atom_site.label_comp_id") if "_atom_site.label_comp_id" in col_names else -1
    
    # Standard amino acids to ensure we're looking at a protein chain
    protein_res = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
                   "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
        
    data_start = start_idx + 1 + len(col_names)
    for i in range(data_start, end_idx):
        line = lines[i].strip()
        if not line or line.startswith("#"): continue
        tokens = line.split()
        if len(tokens) > max(label_idx, auth_idx):
            label = tokens[label_idx]
            auth = tokens[auth_idx]
            
            # Only map if the residue is a standard protein amino acid
            if comp_idx != -1 and len(tokens) > comp_idx:
                comp_id = tokens[comp_idx]
                if comp_id not in protein_res:
                    continue
                    
            if label not in mapping:
                mapping[label] = auth
    return mapping

def _update_struct_ref_seq_pdbx_strand_id(struct_ref_seq_lines, orig_label_to_auth, prody_label_to_auth):
    """Updates _struct_ref_seq.pdbx_strand_id to match the final label_asym_id."""
    if not struct_ref_seq_lines:
        return
        
    orig_auth_to_label = {}
    for label, auth in orig_label_to_auth.items():
        if auth not in orig_auth_to_label:
            orig_auth_to_label[auth] = label
            
    prody_auth_to_label = {}
    for label, auth in prody_label_to_auth.items():
        if auth not in prody_auth_to_label:
            prody_auth_to_label[auth] = label
    
    col_names = []
    header_idx = -1
    for i, line in enumerate(struct_ref_seq_lines):
        stripped = line.strip()
        if stripped.startswith("_struct_ref_seq."):
            col_names.append(stripped)
        elif not stripped.startswith("_") and stripped != "loop_" and stripped != "#":
            header_idx = i
            break
            
    if header_idx != -1 and "_struct_ref_seq.pdbx_strand_id" in col_names:
        target_idx = col_names.index("_struct_ref_seq.pdbx_strand_id")
        for i in range(header_idx, len(struct_ref_seq_lines)):
            line = struct_ref_seq_lines[i].strip()
            if not line or line == "#": continue
            tokens = line.split()
            if len(tokens) > target_idx:
                orig_auth = tokens[target_idx]
                orig_label = orig_auth_to_label.get(orig_auth, orig_auth)
                final_label = prody_auth_to_label.get(orig_label, orig_label)
                tokens[target_idx] = final_label
                struct_ref_seq_lines[i] = " ".join(tokens) + "\n"
    else:
        for i, line in enumerate(struct_ref_seq_lines):
            tokens = line.strip().split()
            if len(tokens) >= 2 and tokens[0] == "_struct_ref_seq.pdbx_strand_id":
                orig_auth = tokens[1]
                orig_label = orig_auth_to_label.get(orig_auth, orig_auth)
                final_label = prody_auth_to_label.get(orig_label, orig_label)
                struct_ref_seq_lines[i] = f"_struct_ref_seq.pdbx_strand_id {final_label}\n"

def _get_asym_id_mapping(orig_label_to_auth, prody_label_to_auth):
    """Generates a mapping from original label_asym_id to processed label_asym_id."""
    label_mapping = {}
    if orig_label_to_auth and prody_label_to_auth:
        for prody_label, prody_auth in prody_label_to_auth.items():
            # In ProDy's output, the auth_asym_id column corresponds to the original label_asym_id.
            orig_label = prody_auth
            if orig_label in orig_label_to_auth:
                label_mapping[orig_label] = prody_label
    return label_mapping

def restore_cif_entity_metadata(original_cif_path, prody_cif_path, final_cif_path):
    """
    Reads the _entity block from the original CIF and appends it 
    to the final CIF file. Also restores the _atom_site.label_entity_id
    and _struct_asym block using mappings cross-checked from the ProDy output.
    Returns a dictionary mapping original label_asym_id to processed label_asym_id.
    """
    if not os.path.exists(original_cif_path) or not os.path.exists(prody_cif_path) or not os.path.exists(final_cif_path): return {}
    
    with open(original_cif_path, 'r', encoding='utf-8') as f: orig_lines = f.readlines()
    with open(prody_cif_path, 'r', encoding='utf-8') as f: prody_lines = f.readlines()
    with open(final_cif_path, 'r', encoding='utf-8') as f: final_lines = f.readlines()

    # 1. Extract _entity and original mapping
    entity_lines = _extract_cif_loop(orig_lines, "_entity.")
    orig_start, orig_end, orig_cols = _get_atom_site_cols(orig_lines)
    orig_asym_to_entity = _get_asym_to_entity_from_atom_site(orig_lines, orig_cols, orig_start, orig_end)
    orig_label_to_auth = _get_asym_mappings_from_atom_site(orig_lines, orig_cols, orig_start, orig_end)

    # 2. Parse ProDy file's _atom_site for label -> auth mapping
    prody_start, prody_end, prody_cols = _get_atom_site_cols(prody_lines)
    prody_label_to_auth = _get_asym_mappings_from_atom_site(prody_lines, prody_cols, prody_start, prody_end)
    
    # 3. Create mapping for final label_asym_id to entity_id
    asym_to_entity = {}
    if prody_label_to_auth and orig_asym_to_entity:
        for prody_label, prody_auth in prody_label_to_auth.items():
            if prody_auth not in orig_asym_to_entity:
                raise ValueError(f"auth_asym_id '{prody_auth}' in the ProDy output does not match any label_asym_id in the original CIF file.")
            asym_to_entity[prody_label] = orig_asym_to_entity[prody_auth]
            if prody_label != prody_auth:
                logging.debug(f"ProDy reassigned label_asym_id from {prody_auth} to {prody_label}. Using the new label.")

    # 4. Fallback inference if missing
    final_start, final_end, final_cols = _get_atom_site_cols(final_lines)
    if not asym_to_entity or not entity_lines:
        logging.warning(f"Inferring entity metadata for {final_cif_path}")
        asym_to_entity, entity_lines = _infer_asym_to_entity(final_lines, final_cols, final_start, final_end)

    # 5. Update final lines
    if asym_to_entity and final_cols:
        _update_atom_site_entity_id(final_lines, final_cols, final_start, final_end, asym_to_entity)
        _update_struct_asym_loop(final_lines, asym_to_entity)
        
    if entity_lines:
        final_lines.extend(["\n#\n"] + entity_lines + ["#\n"])

    struct_ref_lines = _extract_cif_loop(orig_lines, "_struct_ref.")
    if struct_ref_lines:
        final_lines.extend(["\n#\n"] + struct_ref_lines + ["#\n"])
        
    struct_ref_seq_lines = _extract_cif_loop(orig_lines, "_struct_ref_seq.")
    if struct_ref_seq_lines:
        _update_struct_ref_seq_pdbx_strand_id(struct_ref_seq_lines, orig_label_to_auth, prody_label_to_auth)
        final_lines.extend(["\n#\n"] + struct_ref_seq_lines + ["#\n"])
        
    with open(final_cif_path, 'w', encoding='utf-8') as f:
        f.writelines(final_lines)

    return _get_asym_id_mapping(orig_label_to_auth, prody_label_to_auth)


def create_fh_zip_archive(phase1_results, output_dir):
    """Creates a zip archive of the generated FH files for easy cloud upload."""
    import zipfile
    import os
    import logging
    
    zip_path = os.path.join(output_dir, "fh_files.zip")
    try:
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for item in phase1_results:
                _, protein_protonated, _, _ = item
                if os.path.exists(protein_protonated):
                    zipf.write(protein_protonated, os.path.basename(str(protein_protonated)))
        log_step(None, f"Created archive {zip_path} containing FH files for easy cloud upload.", color="white")
    except Exception as e:
        logging.error(f"Failed to create {zip_path}: {e}")


def convert_pdb_to_cif(input_dir):
    """Converts any .pdb files in the input directory to .cif format."""
    from prody import parsePDB, writeMMCIF
    pdb_files_in_dir = glob.glob(os.path.join(input_dir, "*.pdb"))
    if pdb_files_in_dir:
        logging.info(f"Found {len(pdb_files_in_dir)} .pdb file(s) in '{input_dir}'. Converting to .cif format...")
        for pdb_path in pdb_files_in_dir:
            try:
                base_name = os.path.splitext(os.path.basename(pdb_path))[0]
                cif_out_path = os.path.join(input_dir, f"{base_name}.cif")
                if not os.path.exists(cif_out_path):
                    logging.debug(f"Parsing {pdb_path}...")
                    struct = parsePDB(pdb_path)
                    if struct is not None:
                        writeMMCIF(cif_out_path, struct)
                        logging.debug(f"Converted {pdb_path} to {cif_out_path}")
                    else:
                        logging.warning(f"Failed to parse {pdb_path}")
            except Exception as e:
                logging.error(f"Error converting {pdb_path} to CIF: {e}")

def prepare_proteins(input_dir, output_dir, mode, skip_cofactor=False, skip_minimization=False):
    """Main workflow to orchestrate the cleaning process."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.debug(f"Created output directory: {output_dir}")
    
    # Check for PDB files and convert them to CIF first
    convert_pdb_to_cif(input_dir)

    pdb_files = glob.glob(os.path.join(input_dir, "*.cif"))
        
    if not pdb_files:
        logging.error(f"No CIF files found in '{input_dir}'. Please check your folder.")
        return {}
    
    # Partition files: check if they already exist in the output directory
    files_to_clean = []
    already_cleaned_paths = []
    
    for filepath in pdb_files:
        filename = os.path.basename(filepath)
        expected_out_path = os.path.join(output_dir, filename)
        
        if os.path.exists(expected_out_path):
            logging.debug(f"Skipping cleaning for {filename}: already exists in output directory.")
            already_cleaned_paths.append(expected_out_path)
        else:
            files_to_clean.append(filepath)

    # Only run the cleaning modes on files that actually need it
    newly_cleaned_paths = []
    file_stats = {}
    if files_to_clean:
        if skip_cofactor:
            # Skip all mode logic: only remove water, keep all HETATM cofactors
            logging.debug("skip_cofactor=False: Removing only water from all files.")
            for filepath in files_to_clean:
                structure = parseMMCIF(filepath)
                if structure is None:
                    continue
                cleaned_pdb_path = clean_and_save_pdb(structure, filepath, output_dir, "protein or (not water)")
                if cleaned_pdb_path:
                    newly_cleaned_paths.append(cleaned_pdb_path)
                    
                    filename = os.path.basename(filepath)
                    file_hetatms = get_hetatms(structure)
                    file_stats[filename] = generate_file_stats(filepath, file_hetatms, [])
        elif mode == "global":
            newly_cleaned_paths, file_stats = clean_global_mode(files_to_clean, output_dir)
        elif mode == "local":
            newly_cleaned_paths, file_stats = clean_local_mode(files_to_clean, output_dir)
        elif mode == "auto":
            newly_cleaned_paths, file_stats = clean_auto_mode(files_to_clean, output_dir)
            
    # --- Manage CSV Stats ---
    if file_stats:
        import csv
        csv_path = os.path.join(output_dir, "prepare_protein_summary.csv")
        fieldnames = ["Protein", "Kept Ligands", "Removed Ligands", "Suggested Cofactors"]
        existing_stats = {}
        if os.path.exists(csv_path):
            with open(csv_path, "r", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    existing_stats[row["Protein"]] = row
                    
        # Update with new stats
        existing_stats.update(file_stats)
        
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row_key in sorted(existing_stats.keys()):
                writer.writerow(existing_stats[row_key])

    # Combine lists so downstream tasks process all files
    all_cleaned_paths = already_cleaned_paths + newly_cleaned_paths
        
    from pathlib import Path
    prepared_results = {}
    
    # Phase 1: Rebuild and REDUCE2
    phase1_results = []
    for cleaned_pdb_path in all_cleaned_paths:
        protein_path = Path(cleaned_pdb_path)
        protein_dir = protein_path.parent
        protein_base = protein_path.stem

        # Define a new path for the fixed CIF
        protein_fixed = protein_dir / f"{protein_base}_fixed.cif" 
        protein_protonated = protein_dir / f"{protein_base}FH.cif"
        protein_prep_out = protein_dir / protein_base
        protein_pdbqt = protein_dir / f"{protein_base}.pdbqt"

        logging.debug(f"Preparing receptor {cleaned_pdb_path}...")

        # 0. Rebuild Missing Heavy Atoms (PDBFixer)
        if not run_pdbfixer(protein_path, protein_fixed):
            logging.error(f"Failed to fix missing atoms for {protein_base}. Skipping.")
            continue

        # 1. Adding Hydrogens & Optimizing (REDUCE2)
        # Note: We now pass 'protein_fixed' instead of 'protein_path'
        if not run_reduce2(protein_fixed, protein_protonated):
            logging.error(f"Failed to protonate {protein_base}. Skipping.")
            continue
            
        phase1_results.append((protein_base, protein_protonated, protein_prep_out, protein_pdbqt))

    # Phase 1.5: Minimization and User Prompt
    phase15_results = []
    if skip_minimization and phase1_results:
        print()
        log_step("INTERACTIVE", "--skip_minimization is provided.", color="yellow")
        log_step("INTERACTIVE", "The pipeline has stopped to allow manual minimization of the generated FH files (*FH.cif).", color="yellow")
        
        create_fh_zip_archive(phase1_results, output_dir)

        while True:
            ans = input("Have you manually minimized the FH files in the cloud server and replaced the local files? (y/n): ").strip().lower()
            if ans == 'y':
                log_step(None, "Continuing processing...", color="white")
                phase15_results = phase1_results
                break
            else:
                log_step("WARNING", "Please minimize the FH files and replace them in the folder before continuing.", color="yellow")
    else:
        for item in phase1_results:
            protein_base, protein_protonated, protein_prep_out, protein_pdbqt = item
            if not run_openmm_minimization(protein_protonated):
                logging.error(f"Failed to minimize {protein_base}. Skipping.")
                continue
            phase15_results.append(item)

    # Restore _entity metadata from original CIF to the protonated FH file
    new_label_mappings = {}
    for item in phase15_results:
        protein_base, protein_protonated, protein_prep_out, protein_pdbqt = item
        original_cif_path = os.path.join(input_dir, f"{protein_base}.cif")
        prody_cif_path = os.path.join(output_dir, f"{protein_base}.cif")
        mapping = restore_cif_entity_metadata(original_cif_path, prody_cif_path, protein_protonated)
        if mapping:
            new_label_mappings[protein_base] = mapping

    if new_label_mappings:
        import json
        mapping_file = os.path.join(output_dir, "asym_id_mapping.json")
        existing_mappings = {}
        if os.path.exists(mapping_file):
            try:
                with open(mapping_file, "r") as f:
                    existing_mappings = json.load(f)
            except json.JSONDecodeError:
                pass
        existing_mappings.update(new_label_mappings)
        with open(mapping_file, "w") as f:
            json.dump(existing_mappings, f, indent=4)

    # Phase 2: Meeko
    for item in phase15_results:
        protein_base, protein_protonated, protein_prep_out, protein_pdbqt = item
        
        # 2. Preparing Receptor (Meeko)
        if run_meeko_receptor(protein_protonated, protein_prep_out, protein_pdbqt):
            prepared_results[protein_base] = protein_pdbqt
        else:
            logging.error(f"Failed to prepare receptor for {protein_base}. Skipping.")
        
    logging.info("Batch cleaning complete!")
    return prepared_results

def run_scrub_ligand(ligand_path, ligand_scrubbed, ph, generate_isomers):
    """Runs scrub.py to scrub and protonate the ligand."""
    import subprocess
    from pathlib import Path

    if Path(ligand_scrubbed).exists():
        logging.debug(f"Skipping Scrub: {ligand_scrubbed} already exists.")
        return

    cmd_scrub = [
        "scrub.py", str(ligand_path), "-o", str(ligand_scrubbed), 
        "--ph", str(ph)
    ]
    if not generate_isomers:
        cmd_scrub.extend(["--skip_acidbase", "--skip_tautomers"])
    subprocess.run(cmd_scrub, check=True, capture_output=True, text=True)


def run_meeko_ligand(ligand_scrubbed, ligand_pdbqt):
    """Runs mk_prepare_ligand.py (Meeko) to prepare the ligand."""
    import subprocess
    from pathlib import Path
    import logging

    if Path(ligand_pdbqt).exists():
        logging.debug(f"Skipping Meeko: {ligand_pdbqt} already exists.")
        return True

    cmd_meeko = [
        "mk_prepare_ligand.py", "-i", str(ligand_scrubbed), 
        "-o", str(ligand_pdbqt), "--bad_charge_ok"
    ]
    try:
        subprocess.run(cmd_meeko, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Meeko failed to prepare {ligand_scrubbed}:\n{e.stderr}")
        return False

def split_scrubbed_ligand(ligand_scrubbed, output_dir, ligand_base):
    """
    Parses the scrubbed SDF file, finds the lowest energy conformer for each isomer,
    and saves them as separate SDF files.
    """
    from rdkit import Chem
    import json
    from pathlib import Path
    
    suppl = Chem.SDMolSupplier(str(ligand_scrubbed), removeHs=False)
    
    # dictionary: isomerId -> list of (energy, mol)
    isomers = {}
    
    for mol in suppl:
        if mol is None:
            continue
            
        scrub_info_str = mol.GetProp("ScrubInfo") if mol.HasProp("ScrubInfo") else None
        if scrub_info_str:
            scrub_info = json.loads(scrub_info_str)
            isomer_id = scrub_info.get("isomerId", 0)
        else:
            isomer_id = 0
            
        if mol.HasProp("PUBCHEM_MMFF94_ENERGY"):
            energy = float(mol.GetProp("PUBCHEM_MMFF94_ENERGY"))
        else:
            energy = 0.0
                
        if isomer_id not in isomers:
            isomers[isomer_id] = []
        isomers[isomer_id].append((energy, mol))
        
    split_files = []
    out_dir = Path(output_dir)
    
    for isomer_id, confs in isomers.items():
        confs.sort(key=lambda x: x[0])  # sort by energy, lowest first
        best_mol = confs[0][1]
        
        isomer_file = out_dir / f"{ligand_base}_isomer_{isomer_id}.sdf"
        writer = Chem.SDWriter(str(isomer_file))
        writer.write(best_mol)
        writer.close()
        split_files.append(isomer_file)
        
    return split_files

def preprocess_ligand(ligand_path, output_file):
    """
    Pre-processes the ligand:
    - Chooses the largest fragment if there are multiple.
    - Sanitizes the molecule using RDKit.
    - Relaxes the molecule using UFFOptimizeMolecule.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import logging

    try:
        suppl = Chem.SDMolSupplier(str(ligand_path), removeHs=False, sanitize=False)
        mol = None
        for m in suppl:
            if m is not None:
                mol = m
                break
    except Exception as e:
        logging.error(f"Error opening molecule from {ligand_path}: {e}")
        return False
            
    if mol is None:
        logging.error(f"Could not read molecule from {ligand_path}")
        return False

    # Get fragments
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) > 1:
        logging.warning(f"File {ligand_path.name} contains {len(frags)} fragments. Choosing the largest one.")
        largest_frag = max(frags, key=lambda x: x.GetNumHeavyAtoms())
        mol = largest_frag

    # Sanitize
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        logging.warning(f"Sanitization failed for {ligand_path.name}: {e}")
        return False

    try:
        writer = Chem.SDWriter(str(output_file))
        writer.write(mol)
        writer.close()
    except Exception as e:
        logging.error(f"Error writing preprocessed molecule to {output_file}: {e}")
        return False

    return True

def prepare_ligand(ligand_file: str, ph: float, output_dir: str, generate_isomers: bool):
    """
    Prepares the ligand using Molscrub and Meeko.
    Returns the path to the resulting PDBQT file.
    """
    from pathlib import Path

    ligand_path = Path(ligand_file)
    ligand_dir = ligand_path.parent
    ligand_base = ligand_path.stem

    prepared_dir = Path(output_dir)
    prepared_dir.mkdir(parents=True, exist_ok=True)

    ligand_preprocessed = prepared_dir / f"{ligand_base}_preprocessed.sdf"
    ligand_scrubbed = prepared_dir / f"{ligand_base}_scrubbed.sdf"

    logging.debug(f"Preparing ligand {ligand_file}...")

    # 2.5 Preprocessing Ligand (Largest fragment, Sanitize, relax)
    if not preprocess_ligand(ligand_path, ligand_preprocessed):
        logging.error(f"Preprocessing failed for {ligand_base}. Skipping.")
        return []

    # 3. Scrubbing & Protonating Ligand (Molscrub)
    run_scrub_ligand(ligand_preprocessed, ligand_scrubbed, ph, generate_isomers)

    # 3.1 Split Scrubbed Ligand by Isomer and find lowest energy conformer
    isomer_files = split_scrubbed_ligand(ligand_scrubbed, prepared_dir, ligand_base)

    # 4. Preparing Ligand (Meeko)
    pdbqt_files = []
    for isomer_file in isomer_files:
        isomer_base = isomer_file.stem
        ligand_pdbqt = prepared_dir / f"{isomer_base}.pdbqt"
        if run_meeko_ligand(isomer_file, ligand_pdbqt):
            pdbqt_files.append(ligand_pdbqt)

    return pdbqt_files

def prepare_ligands(ligand_path, ph, output_dir, generate_isomers=False):
    """
    Prepares all ligands in the given directory.
    Returns a dictionary mapping ligand names to their PDBQT paths.
    """
    from pathlib import Path
    prepared_ligands = {}
    for lig_file in Path(ligand_path).glob("*.sdf"):
        lig_pdbqts = prepare_ligand(str(lig_file), ph, output_dir, generate_isomers)
        for lig_pdbqt in lig_pdbqts:
            prepared_ligands[lig_pdbqt.stem] = lig_pdbqt
    return prepared_ligands

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean PDB files by removing water and unwanted HETATMs.")
    parser.add_argument("-i", "--input_dir", type=str, default="protein", help="Input directory containing PDB files.")
    parser.add_argument("-o", "--output_dir", type=str, default="protein-clean", help="Output directory for cleaned PDB files.")
    parser.add_argument("-m", "--mode", type=str, choices=["global", "local", "auto"], default="global", help="Elimination mode: 'global' applies to all files, 'local' asks for each file individually, 'auto' automatically keeps essential cofactors based on UniProt DB.")
    
    args = parser.parse_args()
    
    prepare_proteins(input_dir=args.input_dir, output_dir=args.output_dir, mode=args.mode)