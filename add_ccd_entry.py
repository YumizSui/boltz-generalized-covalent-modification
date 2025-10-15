#!/usr/bin/env python
# Add CCD Entry to Boltz Cache
#
# Creates a Chemical Component Dictionary entry from SMILES.
# Simply processes the input SMILES and adds it to the Boltz CCD cache.

import os
import pickle
from pathlib import Path
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem

# THIS IS ABSOLUTELY NECESSARY TO MAKE SURE THAT ALL PROPERTIES ARE PICKLED OR ELSE BOLTZ WILL CRASH.
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

def main():
    parser = argparse.ArgumentParser(description='Add a custom molecule to Boltz CCD cache from SMILES')
    parser.add_argument('smiles', type=str, help='SMILES string of the molecule')
    parser.add_argument('output_ccd_code', type=str, help='New CCD code for the molecule (first 5 letters are written to the CIF file as the resname)')
    parser.add_argument('--force', '-f', action='store_true', help='Force overwrite of existing CCD code')
    parser.add_argument('--cache_dir', type=str, default=os.path.expanduser("~/.boltz"), help='Boltz cache directory')
    parser.add_argument('--export_sdf', action='store_true', help='Export the molecule as an SDF file')
    args = parser.parse_args()

    smiles_str = args.smiles
    output_ccd_code = args.output_ccd_code
    cache_dir = Path(args.cache_dir)

    # Load the Boltz Cache Chemical Component Dictionary
    mol_dir = cache_dir / 'mols'
    ccd_path = cache_dir / 'ccd.pkl'
    with ccd_path.open("rb") as file:
        ccd = pickle.load(file)  # noqa: S301
    if ccd.get(output_ccd_code) is not None and not args.force:
        raise ValueError(f"CCD code {output_ccd_code} already exists in the CCD cache. Use --force to overwrite.")
    print(f"Loaded Boltz CCD cache from {ccd_path}")

    # Parse the SMILES string into a molecule
    smiles_mol = Chem.MolFromSmiles(smiles_str)
    if smiles_mol is None:
        raise ValueError("Invalid SMILES string provided")

    # Remove hydrogens for atom naming
    mol_no_h = AllChem.RemoveHs(smiles_mol)

    # Add the metadata boltz expects to the atoms
    for idx, atom in enumerate(mol_no_h.GetAtoms()):
        # Simple atom naming: ElementSymbol + Index
        default_name = f'{atom.GetSymbol()}{str(idx)}'.upper()

        # Set atom properties
        atom.SetProp('name', default_name)
        atom.SetProp('alt_name', default_name)
        atom.SetBoolProp('leaving_atom', False)

    # Add hydrogens back
    mol_with_h = AllChem.AddHs(mol_no_h)

    # Generate a simple conformer
    embed_success = AllChem.EmbedMolecule(mol_with_h)
    if embed_success == -1:
        print("Warning: Could not embed molecule, trying with random coordinates...")
        AllChem.EmbedMolecule(mol_with_h, useRandomCoords=True)

    optimize_success = AllChem.UFFOptimizeMolecule(mol_with_h)
    if optimize_success == -1:
        print("Warning: UFF optimization failed, structure may not be ideal")

    # Set conformer properties
    for c in mol_with_h.GetConformers():
        c.SetProp('name', 'Ideal')

    # Verify atom properties
    print("\nAtom properties:")
    for atom in mol_with_h.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        props = {k: v for k, v in atom.GetPropsAsDict().items() if k not in ['__computedProps']}
        print(f"Atom {atom.GetIdx()}: {atom.GetSymbol()} - {props}")

    # Save the molecule to the CCD cache
    ccd[output_ccd_code] = mol_with_h

    with ccd_path.open("wb") as file:
        pickle.dump(ccd, file)

    output_mol_path = mol_dir/f"{output_ccd_code}.pkl"
    if output_mol_path.exists() and not args.force:
        raise ValueError(f"Output molecule {output_ccd_code} already exists. Use --force to overwrite.")

    with output_mol_path.open("wb") as file:
        mol_no_hs_final = AllChem.RemoveAllHs(mol_with_h)
        pickle.dump(mol_no_hs_final, file)

    if args.export_sdf:
        with Chem.SDWriter(f"{output_ccd_code}.sdf") as writer:
            writer.write(mol_with_h)

    print(f"\nSuccessfully added {output_ccd_code} to the Boltz CCD cache")

if __name__ == "__main__":
    main()
