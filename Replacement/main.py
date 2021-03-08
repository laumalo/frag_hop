from replace_fragments import PrepareSystem, Replacer
from utils import PDBTools
import os

# Arguments (Example 1)
complex_pdb = \
    '/Users/laura/LAURA/BSC/repos/frag-replacement/TestingFiles/noncovalent_ligand/complex_laura.pdb'

fragment_original_pdb = \
    '/Users/laura/LAURA/BSC/repos/frag-replacement/TestingFiles/noncovalent_ligand/fragment.pdb'

bond_atoms = [['C4', 'H8'], ['C13', 'N2']]


# -------
PDBHandler = PDBTools()

os.makedirs('out', exist_ok=True)
PDBHandler.extract_chain(chain_id='L',
                         input_file=complex_pdb,
                         output_file='out/ligand.pdb')


PrepareSystem(initial_complex=complex_pdb,
              fragment=fragment_original_pdb, ligand_resnum=900,
              bond_atoms=bond_atoms)

Replacer(ligand_pdb='out/ligand.pdb',
         fragment_pdb='out/frag_prepared.pdb',
         bond_atoms=bond_atoms)

