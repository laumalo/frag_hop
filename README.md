
# FragHop: Fragments Hopping


The `FragHop` (Fragments Hopping) is a Python package capable of replacing fragments in a ligand in a protein-ligand complex in order to generate a similar hit molecule.

The current supported fragment replacements are:
* Terminal Fragments

## Examples

### Protein-ligand complex with covalent interactions


```bash
python -m frag_hop.main data/TestingFiles/complexes/covalent.pdb data/TestingFiles/fragments/frag1.pdb -c1 C13-N7 -c2 C1-H4 --covalent
```

### Protein-ligand complex with non-covalent interactions
