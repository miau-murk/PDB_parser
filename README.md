# PDB_parser
Analysis and search of structures in the PDB

This program was created to search for bioisosters in the Protein Data Bank. 
Bioisosters are molecules that differ in one or more atoms or functional groups.

<!--Installation-->
## Installation (Linux)
You must have installed pandas, openpyxl, RDKit. At this stage, this project is not a complete ready-made program.
These are several directories that contain scripts that perform independent tasks, for example, creating sdf files for further analysis.

**Description of directories:**

1. exe_program.zip

The program is designed to generate bioisosters for specified mol files.
You can read more details in the readme file inside the directory.

2. replacement

The folder contains scripts with which to create a bioisoster for each mol file inside the sdf.
This project uses a locally installed database of PCB ligands - ```Components-pub.sdf```
```bioster_replacement.ipynb``` slices ```Components-pub.sdf``` into separate mol files and alternately creates mol files of bioisosters.

The program is used to work with F -- Phenyl substitutions smiles_substitution.py
You can read more about this program in the ```readme_smiles_substitution.md``` file.

The script ``converter.py `` is used to convert ``Components-pub.sdf`` to smiles format.

3. search

```bioster_search.py``` is necessary to determine the identity of the mol files. The algorithm calculates Morgan types based on information about valences, stereochemistry and aromaticity of molecules encrypted in mol files. This code is imported into ```bioster_search.py``` which compares the original ```Components-pub.sdf``` database and the created file with bioisosters. As a result, we get a list of potential bioisosters. 
The ```parser.ipynb``` script for each of the found bioisosters parses information about protein complexes from the PDB, taking into account the selection criteria (for example, matching pH, method, absence of HEM, etc.)

4. sources_file

Local databases

5. created_files

Files created as a result of scripts.

6. result_tables

Summary tables of bioisosters and protein data.

<!--description of commits-->
## description of commits
| Commit   | Description                                                     |
|----------|-----------------------------------------------------------------|
|   	   |                                                                 |
|   	   |                                                                 |
|   	   |                                                                 |
|   	   |                                                                 |
|   	   |                                                                 |
