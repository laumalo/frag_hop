import logging
import sys

logging.basicConfig(format="%(message)s", level=logging.INFO,
                    stream=sys.stdout)

conf_file = '/Users/laura/LAURA/BSC/repos/frag-hop2/frag_hop/data/TestingFiles/test_configurations.conf'
# Example covalent terminal fragment replacement
path_complex = '/Users/laura/LAURA/BSC/repos/frag-hop2/frag_hop/data/TestingFiles/complexes/test_covalent.pdb'

path_fragment = '/Users/laura/LAURA/BSC/repos/frag-hop2/frag_hop/data/TestingFiles/fragments/frag2.pdb'


"""
from frag_hop.helpers.hopping import run_covalent_frag_replacement
run_covalent_frag_replacement(complex_pdb = path_complex,
                              fragment_pdb = path_fragment,
                              conf_file = conf_file,
                              output = 'test_out1',
                              resname = 'GRW')"""


#Example non-covalent terminal fragment replacement
path_complex = '/Users/laura/LAURA/BSC/repos/frag-hop2/frag_hop/data/TestingFiles/complexes/test_noncovalent.pdb'

path_fragment = '/Users/laura/LAURA/BSC/repos/frag-hop2/frag_hop/data/TestingFiles/fragments/frag1.pdb'


from frag_hop.helpers.hopping import run_frag_replacement
run_frag_replacement(complex_pdb = path_complex,
                              fragment_pdb = path_fragment,
                              conf_file = conf_file,
                              output = 'test_out2',
                              chain_id = 'L')

"""
# Example non-covalent scaffold replacement
path_complex = '/Users/laura/LAURA/BSC/repos/frag-hop2/frag_hop/data/TestingFiles/complexes/test_scaffold.pdb'

path_scaffold = '/Users/laura/LAURA/BSC/repos/frag-hop2/frag_hop/data/TestingFiles/scaffolds/scaffold_2FJP.pdb'

from frag_hop.helpers.hopping import run_core_replacement
run_core_replacement(complex_pdb = path_complex,
                                 scaffold_pdb = path_scaffold,
                                 conf_file = conf_file,
                                 output = 'test_out3',
                                 chain_id = 'L')
"""
