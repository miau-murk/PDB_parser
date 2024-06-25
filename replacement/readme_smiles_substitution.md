# Molecule Substitution Script

This script allows for the substitution of parts of molecules defined by SMILES strings in an Excel file. The user specifies a part of the molecule to replace (using a SMARTS pattern), what to replace it with (using a SMILES string), and the script performs the substitution, saving the results in a new Excel file.

## Pre-requisites

Create and activate conda environment
```
conda create --name chem
conda activate chem
```

within the environment install libraries:
```
conda install -c conda-forge rdkit
pip install openpyxl
```

## Usage

1.  Prepare your input Excel file (input.xlsx) with two columns:

- id: Unique identifier for each molecule.
- smiles_original: SMILES string of the original molecule.

2. Run the script from the command line, providing the necessary arguments:

```bash
python smiles_substitution.py -i input.xlsx -o output.xlsx --from "<from_smarts>" --to "<to_smiles>"
```

Where:
-i: Path to the input Excel file.
-o: Path to the output Excel file.
--from: SMARTS pattern of the part of the molecule you want to replace.
--to: SMILES string of the part you want to introduce in place of the part defined by --from.

3. Check the output Excel file (output.xlsx) with three columns:

id: Copied from the input file.
smiles_original: Copied from the input file.
smiles_replaced: SMILES string of the molecule after substitution.

## Notes

- The script assumes the atom in the replacement molecule used for bonding is at index 0. This can be modified in the script if needed.
- The chemical validity of the substitutions is not guaranteed and should be verified separately.


## Examples:

### Pyridine to F-phenyl

need to generate all 3 possible forms of F-phenyl (o, m and p-substituents)

```
python smiles_substitution.py -i input.xlsx -o 4Fphenyl.xlsx --from c1ccccn1 --to c1ccc(F)cc1
python smiles_substitution.py -i input.xlsx -o 3Fphenyl.xlsx --from c1ccccn1 --to c1cc(F)ccc1
python smiles_substitution.py -i input.xlsx -o 2Fphenyl.xlsx --from c1ccccn1 --to c1c(F)cccc1
```

Note: no atom mapping could be specified at this time, so position of C-F might be different from position of N in the ring.

### NO2 -> CF3

```
python smiles_substitution.py -i input.xlsx -o nitroCF3.xlsx --from [NX3+](=O)[O-] --to C(F)(F)F
```

### Br-> CF3

```
python smiles_substitution.py -i input.xlsx -o BrCF3.xlsx --from Br --to C(F)(F)F
```