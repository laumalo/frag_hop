
# FragHop: Fragments Hopping


The `FragHop` (Fragments Hopping) is a Python package capable of replacing fragments in a ligand in a protein-ligand complex in order to generate a similar hit molecule.

The current supported fragment replacements are:
* Terminal Fragments (covalent and non covalent)
* Scaffolds (only non covalent)

## Examples

### Terminal fragments

#### Protein-ligand complex with covalent interactions


```bash
python -m frag_hop.main data/TestingFiles/complexes/test_covalent.pdb data/TestingFiles/fragments/frag1.pdb -c data/TestingFiles/test_configurations.conf --covalent
```

### Protein-ligand complex with non-covalent interactions

```bash
python -m frag_hop.main data/TestingFiles/complexes/test_noncovalent.pdb data/TestingFiles/fragments/frag2.pdb -c data/TestingFiles/test_configurations.conf
```
# Scaffolds
### Protein-ligand complex with non-covalent interactions
```bash
python -m frag_hop.main data/TestingFiles/complexes/test_scaffold.pdb data/TestingFiles/scaffolds/scaffold_2FJP.pdb -c data/TestingFiles/test_configurations.conf --core
```
