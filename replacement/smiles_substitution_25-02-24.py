# %%

import os
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
from pathlib import Path
import argparse
import pandas as pd

# %%

def perform_substitution(molecule, match, replacement_smiles, replacement_bond_atom_idx):
    # Convert the replacement SMILES to a molecule
    replacement_mol = Chem.MolFromSmiles(replacement_smiles)
    if not replacement_mol:
        raise ValueError("Invalid SMILES for replacement")

    # Create editable molecule
    emol = Chem.EditableMol(molecule)

    # Get atoms for the part to replace
    atoms_to_replace = set(match)

    # Find bonds to break in the original molecule
    bonds_to_break = []
    for bond in molecule.GetBonds():
        if bond.GetBeginAtomIdx() in atoms_to_replace and bond.GetEndAtomIdx() not in atoms_to_replace:
            bonds_to_break.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        elif bond.GetEndAtomIdx() in atoms_to_replace and bond.GetBeginAtomIdx() not in atoms_to_replace:
            bonds_to_break.append((bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()))

    # Remove atoms to replace and create a mapping of old indices to new indices
    old_to_new_indices = {}
    atoms_deleted = 0
    
    new_idx = 0
    for old_idx in range(molecule.GetNumAtoms()):
        if old_idx not in atoms_to_replace:
            old_to_new_indices[old_idx] = new_idx
            new_idx = new_idx+1
            
    for old_idx in sorted(atoms_to_replace, reverse=True):
        emol.RemoveAtom(old_idx)

    # Add the atoms and bonds from the replacement molecule
    new_atom_indices = {}  # Map from atom index in replacement_mol to atom index in emol
    for atom in replacement_mol.GetAtoms():
        new_idx = emol.AddAtom(atom)
        new_atom_indices[atom.GetIdx()] = new_idx

    for bond in replacement_mol.GetBonds():
        emol.AddBond(new_atom_indices[bond.GetBeginAtomIdx()],
                     new_atom_indices[bond.GetEndAtomIdx()],
                     bond.GetBondType())

    # Reconnect bonds to the replacement molecule
    for bond in bonds_to_break:
        original_atom_idx = bond[1]
        # Map the original index to the new index
        new_atom_idx = old_to_new_indices.get(original_atom_idx, original_atom_idx)

        # Check if the bond formation is chemically valid
        if new_atom_idx == new_atom_indices[replacement_bond_atom_idx]:
            print("Attempting to form a self-bond. Skipping.")
            continue

        # Add bond
        emol.AddBond(new_atom_idx, new_atom_indices[replacement_bond_atom_idx], rdchem.BondType.SINGLE)  # Adjust the bond type as needed

    # Get the final modified molecule
    modified_mol = emol.GetMol()

    # Optionally: sanitize the molecule, compute 2D coordinates for visualization
    Chem.SanitizeMol(modified_mol)
    AllChem.Compute2DCoords(modified_mol)

    return modified_mol



def substitute_smiles(input_smiles, part_to_replace_smarts, replacement_smiles, replacement_bond_atom_idx):
  
    # Read the input molecule from a .mol file
    mol = Chem.MolFromSmiles(input_smiles)    
    if not mol:
        raise ValueError("Invalid input molecule file")

    # Find the substructures to replace
    query = Chem.MolFromSmarts(part_to_replace_smarts)
    matches = mol.GetSubstructMatches(query)

    # Filter matches based on the number of outgoing bonds
    filtered_matches = []
    for match in matches:
        num_outgoing_bonds = sum(1 for bond in mol.GetBonds() if (bond.GetBeginAtomIdx() in match and bond.GetEndAtomIdx() not in match) 
                                                           or (bond.GetEndAtomIdx() in match and bond.GetBeginAtomIdx() not in match))
        if num_outgoing_bonds == 1:  # Only allow substitution if there's exactly one outgoing bond
            filtered_matches.append(match)

    if not filtered_matches:
        # No matching substructure with the correct number of outgoing bonds found
        return []

    modified_smiles = []
    for i, match in enumerate(filtered_matches):
        try:
            # Perform the substitution
            modified_mol = perform_substitution(mol, match, replacement_smiles, replacement_bond_atom_idx)

            # Convert the modified molecule to SMILES
            modified_smiles_str = Chem.MolToSmiles(modified_mol)
            modified_smiles.append(modified_smiles_str)
        except Exception as e:
            print(f"Error during substitution: {e}")

    return modified_smiles


def write_file(name, df_cont):
    file = open(name, "w")
    for i in df_cont:
        file.writelines(i["smiles_replaced"] + "    " + i["id"] + "\n")
    file.close()
        
    
# %%

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Replace parts of molecules specified in an Excel file.')
    parser.add_argument('-i', type=str, required=True, help='Input Excel file with list of SMILES')
    parser.add_argument('-o', type=str, required=True, help='Output Excel file with list of generated SMILES')
    parser.add_argument('--from', dest='from_smarts', type=str, required=True, help='SMARTS pattern for replacement')
    parser.add_argument('--to', dest='to_smiles', type=str, required=True, help='SMILES to replace to')
    args = parser.parse_args()

    # Read the input Excel file
    input_df = pd.read_excel(args.i)
    if 'id' not in input_df.columns or 'smiles_original' not in input_df.columns:
        raise ValueError("Input Excel must contain 'id' and 'smiles_original' columns")

    # Prepare the output DataFrame
    output_df = []

    # Iterate over each row in the input DataFrame
    for idx, row in input_df.iterrows():
        original_smiles = row['smiles_original']
        molecule_id = row['id']
        # Perform the substitution
        replaced_smiles_list = substitute_smiles(original_smiles, args.from_smarts, args.to_smiles, 0)  # Assuming bond atom index 0

        # Append results to the output DataFrame
        for replaced_smiles in replaced_smiles_list:
            output_df.append({'id': molecule_id, 'smiles_original': original_smiles, 'smiles_replaced': replaced_smiles})
            
            
    
    #output_df = pd.DataFrame(output_df)
    print(f"{len(output_df)} replacements were made")

    # Write the output DataFrame to an Excel file
    #output_df.to_csv(args.o, index=False)
    write_file(args.o, output_df)
    print(f'Results saved to {args.o}')

if __name__ == "__main__":
    main()

