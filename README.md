# Run Boltz-1/2 Model With Arbitrary Non-Canonical Amino Acids

A jupyter notebook which walks you through injecting a new residue into the CCD cache and then using
that new residue as a modification in a boltz yaml file.
All you need for the input is a smiles string representing the molecule attached to an amino acid.

Runs with any python environment that can run Boltz inference.

## Changes
Forked from [https://github.com/benf549/boltz-generalized-covalent-modification](https://github.com/benf549/boltz-generalized-covalent-modification). Made it executable with `covalent_inference.py`. Also supports Boltz-2.

### Register CCD
Simply register a new molecule to CCD from SMILES string. Can be used for covalent ligands etc.

```
python src/add_ccd_entry.py \
    "C([C@@H]1[C@H]([C@@H]([C@@H]([C@H](O1)OC[C@@H](C(=O)O)N)O)O)O)O" \
    "SEMAN"
```

### Register CCD aligned to amino acid
Register a new molecule to CCD from SMILES string while aligning to existing amino acid structure. Can be used for modifications.
```
python src/covalent_inference.py \
    "C[C@H]([C@@H](C(=O)O)N)O[C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O" \
    "THR" \
    "THMAN"
```
