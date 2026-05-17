import os

# Dynamically get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Look for the input file in that same directory
input_file = os.path.join(script_dir, "ligand", "PubChem_search_records.sdf")

with open(input_file, "r", encoding="utf-8") as f:
    current_molecule = []
    
    for line in f:
        current_molecule.append(line)
        
        # '$$$$' indicates the end of a molecular record in an SDF file
        if line.strip() == "$$$$":
            # PubChem puts the CID on the very first line of the block
            cid_name = current_molecule[0].strip()
            
            if not cid_name:
                cid_name = "unknown_compound"
                
            # Clean up filename for safety
            filename = "".join(c for c in cid_name if c.isalnum() or c in ("_", "-"))
            
            # Set the output path directly to the ligand directory
            output_path = os.path.join(script_dir, "ligand", f"{filename}.sdf")
            
            # Save the individual compound
            with open(output_path, "w", encoding="utf-8") as out_f:
                out_f.writelines(current_molecule)
            
            # Reset for the next molecule
            current_molecule = []

print(f"Done! Separate files have been saved directly to: {script_dir}")