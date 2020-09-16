# sap_redesign

A tool for resurfacing proteins to increase their solubility.
Requires numpy and pyrosetta to be installed.

## Installation

TODO

## Usage

Redesign a folder full of pdbs, picking the worst 20 residues by SAP score:
```bash
~/sap_redesign/sap_redesign/redesign.py -worst_n 20 -redesign_above 0.3 --pdbs *.pdb
```
Redesign a single pdb, picking the worst 40 residues:
```bash
TODO
```
Redesign a all pdbs in a silent file:
```bash
TODO
```

## Contributing
Pull requests are welcome. 

## License
[MIT](https://choosealicense.com/licenses/mit/)
