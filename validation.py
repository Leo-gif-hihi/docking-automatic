import os
import sys
import glob
import time
import shutil
import logging
import csv
from pathlib import Path

from prody import parsePDB, parseMMCIF, writePDB, writeMMCIF
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from logger_utils import log_step
from prepare import prepare_proteins, prepare_ligands
from download_protein import download_proteins


def calculate_rmsd(ref_sdf, docked_sdf):
    """Calculates the best RMSD between the reference ligand and docked poses."""
    try:
        ref_supplier = Chem.SDMolSupplier(str(ref_sdf))
        ref_mol = ref_supplier[0] if len(ref_supplier) > 0 else None
        
        if not ref_mol:
            logging.error(f"Could not read reference SDF: {ref_sdf}")
            return None
            
        # Try to sanitize, sometimes PDB to SDF lacks proper sanitization
        try:
            Chem.SanitizeMol(ref_mol)
        except:
            pass

        docked_supplier = Chem.SDMolSupplier(str(docked_sdf))
        
        # Only evaluate the top-ranked pose (index 0)
        top_pose = docked_supplier[0] if len(docked_supplier) > 0 else None
        
        if not top_pose:
            logging.error(f"Could not read docked poses from: {docked_sdf}")
            return None
            
        try:
            Chem.SanitizeMol(top_pose)
        except:
            pass
            
        try:
            # Use GetBestRMS to account for automorphisms/symmetry on the top pose
            rmsd = rdMolAlign.GetBestRMS(ref_mol, top_pose)
            return rmsd
        except Exception as e:
            logging.error(f"RMSD calculation failed for top pose: {e}")
            return None
    except Exception as e:
        logging.error(f"Error calculating RMSD: {e}")
        return None

def generate_rmsd_plot(csv_path, output_dir):
    """Generates a grouped bar chart of RMSD values across validation runs."""
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        df = pd.read_csv(csv_path)
        if df.empty:
            return
        # Combine Complex and Ligand to create a unique identifier for each validation target
        df['Target'] = df['Complex'] + "_" + df['Ligand'].str.replace('_isomer_0', '', regex=False)
            
        plt.figure(figsize=(max(10, len(df['Target'].unique()) * 2), 6))
        sns.set_theme(style="whitegrid")
        
        # Create a box plot to show the distribution of RMSD across runs for each target
        ax = sns.boxplot(
            data=df,
            x='Target',
            y='RMSD',
            hue='Target',
            palette='viridis',
            legend=False
        )
        
        # Add a horizontal line at 2.0 Å threshold
        plt.axhline(y=2.0, color='red', linestyle='--', linewidth=2, label='2.0 Å Threshold')
        
        plt.title('Validation RMSD Distribution per Target', fontsize=16)
        plt.xlabel('Target', fontsize=14)
        plt.ylabel('RMSD (Å)', fontsize=14)
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        plot_path = Path(output_dir) / 'validation_visualization' / 'rmsd_plot.png'
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        from logger_utils import log_step
        log_step("VALIDATION", f"RMSD box plot saved to {plot_path}", color="cyan")
        
    except ImportError:
        logging.warning("Could not generate RMSD plot: pandas, matplotlib, or seaborn is not installed.")
    except Exception as e:
        logging.warning(f"Failed to generate RMSD plot: {e}")

def run_validation_pipeline(args):
    """Runs the validation pipeline."""
    from main import generate_docking_jobs, run_docking_pipeline
    
    input_path = Path(args.protein_input)
    
    if input_path.is_dir():
        protein_path = input_path
    else:
        protein_path = Path(input_path.stem)
        if not protein_path.exists():
            os.makedirs(protein_path, exist_ok=True)
            log_step("WORKFLOW", f"Downloading proteins listed in {input_path}...")
            download_proteins(str(input_path), str(protein_path))
            
    complex_files = list(protein_path.glob("*.pdb")) + list(protein_path.glob("*.cif"))
        
    if not complex_files:
        logging.error(f"No valid complex files found in {protein_path}")
        return
        
    if not args.output_dir:
        output_base = protein_path.stem
        args.output_dir = f"validation_output_{output_base}"
        
    os.makedirs(args.output_dir, exist_ok=True)
    
    temp_protein_dir = Path(args.output_dir) / "validation_protein_temp"
    temp_ligand_dir = Path(args.output_dir) / "validation_ligand_temp"
    box_dir = Path(args.output_dir) / "validation_box"
    val_vis_dir = Path(args.output_dir) / "validation_visualization"
    
    os.makedirs(temp_protein_dir, exist_ok=True)
    os.makedirs(temp_ligand_dir, exist_ok=True)
    os.makedirs(box_dir, exist_ok=True)
    os.makedirs(val_vis_dir, exist_ok=True)
    
    positive_control_map = {}
    original_ligand_sdfs = {}
    validation_context = {}
    
    rmsd_results = []
    
    for complex_file in complex_files:
        log_step("VALIDATION", f"Processing {complex_file.name}...")
        
        # Parse the complex
        try:
            if complex_file.suffix.lower() == '.cif':
                structure = parseMMCIF(str(complex_file))
            else:
                structure = parsePDB(str(complex_file))
        except Exception as e:
            logging.error(f"Failed to parse {complex_file}: {e}")
            continue
            
        if not structure:
            logging.error(f"Parsed structure is empty for {complex_file}")
            continue
            
        # Find potential ligands (HETATM excluding water/ions)
        hetatms = structure.select('not protein and not water and not ion')
        if not hetatms:
            logging.warning(f"No non-water/ion HETATMs found in {complex_file.name}. Skipping.")
            continue
            
        # Get unique (chain, resnum, resname) instances
        ligand_instances = list(set(zip(hetatms.getChids(), hetatms.getResnums(), hetatms.getResnames())))
        # Sort them for deterministic ordering
        ligand_instances.sort(key=lambda x: (x[0], x[1], x[2]))
        
        if not ligand_instances:
            logging.warning(f"No valid ligand residues found in {complex_file.name}. Skipping.")
            continue
            
        target_instances = []
        if len(ligand_instances) > 1:
            log_step("INTERACTIVE", f"Multiple ligand instances found in {complex_file.name}:")
            print("0: Validate ALL ligands")
            instance_strs = [f"{i+1}: Chain {chain}, Res {resnum} ({resname})" for i, (chain, resnum, resname) in enumerate(ligand_instances)]
            for s in instance_strs:
                print(s)
                
            log_step("WARNING", "If you select a cofactor as a target, it will be completely removed from the protein structure during validation.", color="yellow")
            
            while True:
                user_input = input("Enter the number(s) of the ligand(s) you want to validate, separated by commas (default 0 for ALL): ").strip()
                if not user_input:
                    user_input = "0"
                
                try:
                    selected_indices = [int(x.strip()) for x in user_input.split(',')]
                    
                    if 0 in selected_indices:
                        target_instances = ligand_instances
                        break
                    
                    valid_selection = True
                    temp_targets = []
                    for idx in selected_indices:
                        if 1 <= idx <= len(ligand_instances):
                            temp_targets.append(ligand_instances[idx - 1])
                        else:
                            valid_selection = False
                            break
                            
                    if valid_selection and temp_targets:
                        target_instances = temp_targets
                        break
                    else:
                        print(f"Invalid selection. Please enter numbers between 0 and {len(ligand_instances)}.")
                except ValueError:
                    print("Please enter valid numbers separated by commas (e.g., 1, 3).")
        else:
            target_instances = ligand_instances
            log_step("VALIDATION", f"Auto-selected ligand: Chain {target_instances[0][0]}, Res {target_instances[0][1]} ({target_instances[0][2]})")
            
        complex_base = complex_file.stem
        
        # Build exclusion string for all targets to prepare protein ONCE
        exclusion_clauses = []
        for instance in target_instances:
            c, r, n = instance
            exclusion_clauses.append(f"(chain {c} and resnum {r} and resname {n})")
        exclusion_str = " or ".join(exclusion_clauses)
        
        protein_sel = structure.select(f'not ({exclusion_str})')
        protein_cif_path = temp_protein_dir / f"{complex_base}.cif"
        writeMMCIF(str(protein_cif_path), protein_sel)
        log_step("VALIDATION", f"Extracted protein (all targets removed) saved to {protein_cif_path}")
        
        for target_instance in target_instances:
            target_chain, target_resnum, target_resname = target_instance
                
            # Extract ligand
            ligand_sel = structure.select(f'chain {target_chain} and resnum {target_resnum} and resname {target_resname}')
            ligand_base_name = f"{target_resname}_{target_chain}_{target_resnum}"
            ligand_pdb_path = temp_ligand_dir / f"{ligand_base_name}.pdb"
            writePDB(str(ligand_pdb_path), ligand_sel)
            
            # Convert ligand PDB to SDF using RDKit
            ligand_sdf_path = temp_ligand_dir / f"{ligand_base_name}.sdf"
            try:
                mol = Chem.MolFromPDBFile(str(ligand_pdb_path), sanitize=False)
                if mol:
                    writer = Chem.SDWriter(str(ligand_sdf_path))
                    writer.write(mol)
                    writer.close()
                    log_step("VALIDATION", f"Extracted ligand saved to {ligand_sdf_path}")
                else:
                    logging.error(f"RDKit failed to parse the extracted ligand PDB for {ligand_base_name}.")
                    continue
            except Exception as e:
                logging.error(f"Error converting ligand PDB to SDF for {ligand_base_name}: {e}")
                continue
                
            # Name the protein specifically for this ligand (we will copy the prepared protein later)
            specific_protein_base = f"{complex_base}_{ligand_base_name}"
            
            box_file = box_dir / f"{specific_protein_base}.box.txt"
            
            custom_box_file = None
            if getattr(args, 'box_dir', None):
                potential_custom_box = Path(args.box_dir) / f"{specific_protein_base}.box.txt"
                
                if potential_custom_box.exists():
                    custom_box_file = potential_custom_box
                    
            if custom_box_file:
                shutil.copy(custom_box_file, box_file)
                log_step("VALIDATION", f"Using custom box file: {custom_box_file}")
            else:
                # Get box size from user
                log_step("INTERACTIVE", f"We need the size of the box around {ligand_base_name}.")
                box_size_input = input("Enter box size (default 20x20x20): ").strip()
                
                size_x = size_y = size_z = 20.0
                if box_size_input:
                    parts = box_size_input.replace('x', ' ').replace(',', ' ').split()
                    try:
                        if len(parts) == 1:
                            size_x = size_y = size_z = float(parts[0])
                        elif len(parts) == 3:
                            size_x, size_y, size_z = map(float, parts)
                        else:
                            log_step("WARNING", "Invalid input format. Using default 20x20x20.", color="yellow")
                    except ValueError:
                        log_step("WARNING", "Could not parse numbers. Using default 20x20x20.", color="yellow")
                        
                if size_x <= 0 or size_y <= 0 or size_z <= 0:
                    log_step("WARNING", "Box size must be greater than 0. Falling back to default 20x20x20.", color="yellow")
                    size_x = size_y = size_z = 20.0
                    
                # Calculate ligand center
                coords = ligand_sel.getCoords()
                center_x, center_y, center_z = coords.mean(axis=0)
                
                # Write box file
                with open(box_file, 'w') as f:
                    f.write(f"center_x = {center_x:.3f}\n")
                    f.write(f"center_y = {center_y:.3f}\n")
                    f.write(f"center_z = {center_z:.3f}\n")
                    f.write(f"size_x = {size_x:.3f}\n")
                    f.write(f"size_y = {size_y:.3f}\n")
                    f.write(f"size_z = {size_z:.3f}\n")
                    f.write("exhaustiveness = 8\n") # Default
                    
                log_step("VALIDATION", f"Created box file: {box_file}")
            
            # Update positive control map and original ligand sdfs
            if ligand_base_name.lower() not in positive_control_map:
                positive_control_map[ligand_base_name.lower()] = []
            positive_control_map[ligand_base_name.lower()].append(specific_protein_base.lower())
            
            original_ligand_sdfs[ligand_base_name] = ligand_sdf_path
            
            validation_context[specific_protein_base] = {
                'complex_path': str(complex_file.resolve()),
                'resname': target_resname,
                'chain': target_chain,
                'resnum': target_resnum
            }
        
    # Load flex_res_map if it exists in args
    flex_res_map = {}
    if getattr(args, 'flex_res', None) and os.path.exists(args.flex_res):
        try:
            with open(args.flex_res, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    prot = row.get('protein', '').strip()
                    flex = row.get('flex_res', '').strip()
                    if prot and flex:
                        flex_res_map[prot.lower()] = flex
        except Exception as e:
            logging.error(f"Error reading flexible residues file in validation: {e}")

    # Prepare Protein
    protein_clean_dir = Path(args.output_dir) / "validation_protein_prepared"
    prepared_proteins = prepare_proteins(
        input_dir=str(temp_protein_dir),
        output_dir=str(protein_clean_dir),
        mode=args.clean_mode,
        skip_cofactor=args.skip_cofactor,
        skip_minimization=args.skip_minimization,
        flex_res_map=flex_res_map
    )
    
    # Duplicate prepared proteins for each specific ligand target so they have distinct boxes and jobs
    new_prepared_proteins = {}
    for specific_protein_base, context in validation_context.items():
        complex_base = Path(context['complex_path']).stem
        if complex_base in prepared_proteins:
            orig_pdbqt_dict = prepared_proteins[complex_base]
            new_pdbqt_dict = {}
            
            if orig_pdbqt_dict.get("rigid"):
                orig_rigid = Path(orig_pdbqt_dict["rigid"])
                new_rigid = orig_rigid.parent / (f"{specific_protein_base}_rigid.pdbqt" if orig_pdbqt_dict.get("flex") else f"{specific_protein_base}.pdbqt")
                shutil.copy(orig_rigid, new_rigid)
                new_pdbqt_dict["rigid"] = new_rigid
                
            if orig_pdbqt_dict.get("flex"):
                orig_flex = Path(orig_pdbqt_dict["flex"])
                new_flex = orig_flex.parent / f"{specific_protein_base}_flex.pdbqt"
                shutil.copy(orig_flex, new_flex)
                new_pdbqt_dict["flex"] = new_flex
                
            new_prepared_proteins[specific_protein_base] = new_pdbqt_dict
    
    prepared_proteins = new_prepared_proteins
    
    # Prepare Ligand
    ligand_prepared_dir = Path(args.output_dir) / "validation_ligand_prepared"
    prepared_ligands = prepare_ligands(
        ligand_path=str(temp_ligand_dir),
        ph=args.ph,
        output_dir=str(ligand_prepared_dir),
        generate_isomers=args.generate_isomers
    )
    
    # Docking
    jobs_list = list(generate_docking_jobs(prepared_proteins, prepared_ligands, box_dir, args.num_runs, positive_control_map))
    
    # Filter jobs to only include the specific ligands we targeted in this validation run 
    # (prevents leftover ligands from previous runs in the output dir from being docked)
    valid_ligand_bases = set(original_ligand_sdfs.keys())
    filtered_jobs = []
    for job in jobs_list:
        ligand_base = job[1]
        base_ligand = ligand_base.split("_isomer_")[0] if "_isomer_" in ligand_base else ligand_base
        if base_ligand in valid_ligand_bases:
            filtered_jobs.append(job)
            
    jobs_list = filtered_jobs
    
    total_jobs = len(jobs_list)
    if total_jobs == 0:
        logging.error("No valid docking jobs generated.")
        return
        
    log_step("VALIDATION", f"Starting docking for {len(complex_files)} complexes...")
    for i, (protein_base, ligand_base, box_file, run_index) in enumerate(jobs_list, 1):
        protein_pdbqt_dict = prepared_proteins.get(protein_base)
        ligand_pdbqt = prepared_ligands.get(ligand_base)
        
        base_ligand = ligand_base.split("_isomer_")[0] if "_isomer_" in ligand_base else ligand_base
        suffix = f"_{base_ligand}"
        actual_complex_base = protein_base[:-len(suffix)] if protein_base.endswith(suffix) else protein_base
        
        vina_out_dir = Path(args.output_dir) / "vina_output" / f"run_{run_index}"
        complex_output_dir = vina_out_dir / f"{actual_complex_base}_{ligand_base}"
        os.makedirs(complex_output_dir, exist_ok=True)
        
        success = run_docking_pipeline(
            protein_pdbqt_dict, ligand_pdbqt, box_file, str(complex_output_dir),
            actual_complex_base, ligand_base, run_index, args.cpus, args.exhaustiveness
        )
        
        if success:
            # Calculate RMSD
            docked_sdf = complex_output_dir / f"{actual_complex_base}_{ligand_base}_vina_out.sdf"
            
            base_ligand = ligand_base.split("_isomer_")[0] if "_isomer_" in ligand_base else ligand_base
            ref_ligand_sdf_path = original_ligand_sdfs.get(base_ligand)
            
            if not ref_ligand_sdf_path:
                logging.error(f"Could not find reference ligand SDF for {base_ligand}")
                continue
                
            rmsd = calculate_rmsd(ref_ligand_sdf_path, docked_sdf)
            if rmsd is not None:
                log_step("VALIDATION", f"[{protein_base}] Run {run_index} RMSD: {rmsd:.3f} Å", color="green")
                rmsd_results.append({
                    'Complex': actual_complex_base,
                    'Ligand': ligand_base,
                    'Run': run_index,
                    'RMSD': rmsd
                })
            else:
                log_step("WARNING", f"Failed to calculate RMSD for Run {run_index}.", color="yellow")
                
            # Generate PyMOL script
            context = validation_context.get(protein_base)
            if context and docked_sdf.exists():
                pml_script_path = val_vis_dir / f"visualize_{actual_complex_base}_{ligand_base}_run{run_index}.pml"
                
                complex_name = Path(context['complex_path']).stem
                docked_name = docked_sdf.stem
                resn = context['resname']
                chain = context['chain']
                resi = context['resnum']
                
                # Extract top pose to a separate file for PyMOL
                top_pose_sdf = complex_output_dir / f"{actual_complex_base}_{ligand_base}_top_pose.sdf"
                try:
                    supplier = Chem.SDMolSupplier(str(docked_sdf))
                    if len(supplier) > 0 and supplier[0] is not None:
                        writer = Chem.SDWriter(str(top_pose_sdf))
                        writer.write(supplier[0])
                        writer.close()
                except Exception as e:
                    logging.warning(f"Failed to extract top pose for PyMOL visualization: {e}")
                    
                with open(pml_script_path, "w") as f:
                    f.write(f"load {context['complex_path']}, {complex_name}\n")
                    f.write(f"load {ref_ligand_sdf_path.resolve()}, ref_ligand\n")
                    
                    # Load the newly extracted single-pose SDF
                    if top_pose_sdf.exists():
                        docked_name = top_pose_sdf.stem
                        f.write(f"load {top_pose_sdf.resolve()}, {docked_name}\n")
                    else:
                        docked_name = docked_sdf.stem
                        f.write(f"load {docked_sdf.resolve()}, {docked_name}\n")
                    
                    f.write(f"hide everything\n")
                    f.write(f"show cartoon, {complex_name}\n")
                    f.write(f"color green, {complex_name}\n")
                    
                    f.write(f"show sticks, ref_ligand\n")
                    f.write(f"color yellow, ref_ligand\n")
                    
                    f.write(f"show sticks, {docked_name}\n")
                    f.write(f"color cyan, {docked_name}\n")
                    
                    f.write(f"center {docked_name}\n")
                    f.write(f"zoom {docked_name}, 10\n")
                    
                    if rmsd is not None:
                        f.write(f"print('Python RDKit RMSD: {rmsd:.3f} A')\n")

    # Save RMSD results
    if rmsd_results:
        csv_path = Path(args.output_dir) / "validation_rmsd.csv"
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['Complex', 'Ligand', 'Run', 'RMSD'])
            writer.writeheader()
            writer.writerows(rmsd_results)
        log_step("VALIDATION", f"Validation RMSD results saved to {csv_path}", color="green")
        
        # Generate the RMSD visualization plot
        generate_rmsd_plot(csv_path, args.output_dir)
        
    log_step("VALIDATION", f"PyMOL visualization scripts have been saved to the '{val_vis_dir}' directory.", color="cyan")
