# sap_redesign

A tool for resurfacing proteins to increase their solubility.
Requires numpy and pyrosetta to be installed.

## Installation

TODO

## Usage
To use this you must provide at minimum the input filename(s) and a per 
residue score cutoff. Phil found a cutoff of above 0.3 seems to pick off any 
hydrophobics that are even a little bit exposed, if you are a little careful,
0.6 seems to pick off most exposed hydrophobics while still leaving boundary 
regions for the most part. If you just want to clean the absolute worst 
offenders, 1.0 will do that.

Redesign a folder full of pdbs, picking the worst 20 residues by SAP score:
```bash
~/sap_redesign/sap_redesign/redesign.py --pdbs *.pdb --worst_n 20 --redesign_above 0.3
```
Redesign a single pdb, without touching HBNet polars or PRO and GLY positions,
and designing only 10 residues at a time to increase speed:
```bash
~/sap_redesign/sap_redesign/redesign.py --in:file:silent in.silent  --redesign_above 0.6 --lock_PG --lock_HNQST --chunk

```
Redesign a all pdbs in a silent file, removing only the worst residues,
favoring mutation and penalizing ARG:
```bash
~/sap_redesign/sap_redesign/redesign.py --in:file:silent in.silent  --redesign_above 1.0 --encourage_mutation --penalize_ARG
```

## Contributing
Pull requests are welcome. 

## License
[MIT](https://choosealicense.com/licenses/mit/)
