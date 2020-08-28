# sap_redesign

A tool for resurfacing proteins to increase their solubility.
Requires numpy and pyrosetta to be installed.

## Installation

TODO

## Usage

Redesign a folder full of pdbs, picking the worst 20 residues by SAP score:
```bash
python redesign.py -worst_n 20 --pdbs *.pdb
```
Redesign a single pdb, picking the worst 40 residues:
```bash
python redesign.py -worst_n 40 --pdbs A_REALLY_STICKY_PROTEIN.pdb 
```
Redesign a all pdbs in a silent file:
```bash
python redesign.py -in:file:silent in.silent
```

## Contributing
Pull requests are welcome. 

## License
[MIT](https://choosealicense.com/licenses/mit/)
