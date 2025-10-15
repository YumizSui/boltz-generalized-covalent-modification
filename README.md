# Run Boltz-1 Model With Arbitrary Non-Canonical Amino Acids

A jupyter notebook which walks you through injecting a new residue into the CCD cache and then using
that new residue as a modification in a boltz yaml file.
All you need for the input is a smiles string representing the molecule attached to an amino acid.

Runs with any python environment that can run Boltz inference.

## Changes
Forked from [https://github.com/benf549/boltz-generalized-covalent-modification](https://github.com/benf549/boltz-generalized-covalent-modification). Made it executable with `covalent_inference.py`. Also supports Boltz-2.

### Example 1
```
python src/covalent_inference.py \
    "C[C@H]([C@@H](C(=O)O)N)O[C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O" \
    "THR" \
    "THMAN"
```

### Example 2
```
python src/covalent_inference.py \
    "C([C@@H]1[C@H]([C@@H]([C@@H]([C@H](O1)OC[C@@H](C(=O)O)N)O)O)O)O" \
    "SER" \
    "SEMAN"
```
