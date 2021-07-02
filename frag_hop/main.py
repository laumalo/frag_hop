"""
This module is designed to run FragHop through the command-line.
"""

__author__ = "Laura Malo"
__license__ = "GPL"
__maintainer__ = "Laura Malo"
__email__ = "laura.maloroset@bsc.es"


# General imports
import argparse as ap
import logging
import sys

logging.basicConfig(format="%(message)s", level=logging.INFO,
                    stream=sys.stdout)

def parse_args():
    """
    It parses the command-line arguments.
    Parameters
    ----------
    args : list[str]
        List of command-line arguments to parse
    Returns
    -------
    parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """
    parser = ap.ArgumentParser(description=
                               "FRAGHOP: PERFORM FRAGMENT REPLACEMENT IN A " +
                               "PROTEIN-LIGAND COMPLEX")
    parser.add_argument("complex_pdb", type=str,
                        help="Path to the protein-ligand complex.")

    parser.add_argument("fragment_pdb", type=str,
                        help="Path to the hit fragment.")

    parser.add_argument("-c", "--conf_file", type=str,
                        help="Path to the configuration file.")

    parser.add_argument("-o", "--output", type=str, default='out',
                        help="Path to the output folder. Default: out.")

    parser.add_argument("--covalent",
                        dest="covalent",
                        action='store_true',
                        help="Perform the fragment replacement in a covalent" +
                        " complex.")

    parser.add_argument("--core",
                        dest="core",
                        action='store_true',
                        help="Perform a scaffold replacement.")

    parser.set_defaults(covalent=False)
    parser.set_defaults(core=False)

    parsed_args = parser.parse_args()
    return parsed_args



def main(args):
    """
    It reads the command-line arguments and runs the fragment replacement
    protocol.

    Parameters
    ----------
    args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user

    Examples
    --------

    From the command-line:

    >>> python -m frag_hop.main data/TestingFiles/complexes/covalent.pdb data/TestingFiles/fragments/frag1.pdb -c1 C13-N7 -c2 C1-H4 --covalent

    """
    CHAIN_ID = 'L' # Default parameters for non-covalent ligands
    RESNAME_COV = 'GRW' # Default parameters for covalent ligands

    type_rep = 'scaffold' if args.core else 'terminal'
    type_ligand = 'covalent' if args.covalent else 'non-covalent'

    from frag_hop.utils import parse_conf_file
    dict_connectivities = parse_conf_file(file_path=args.conf_file)
    BONDS_FRAG = dict_connectivities.get(args.fragment_pdb)
    BONDS_LIG = dict_connectivities.get(args.complex_pdb)

    logging.info('-' * 75)
    logging.info('Fragment replacement for a protein-ligand complex')
    logging.info('-' * 75)
    logging.info(' - Input information:')
    logging.info('    - Protein-ligand complex PDB: %s',args.complex_pdb)
    logging.info('    - Hit fragment PDB: %s', args.fragment_pdb)
    logging.info('    - Type of fragment replacement: %s', type_rep)
    logging.info('    - Type of ligand: %s', type_ligand)
    logging.info('    - Connectivity remaining ligand: %s', BONDS_LIG)
    logging.info('    - Connectivity fragment: %s', BONDS_FRAG)
    logging.info('-' * 75)

    # Scaffold hopping
    if args.core:
        if args.covalent:
            raise NotImplementedError
        else:
            from frag_hop.helpers.hopping import run_core_replacement
            run_core_replacement(complex_pdb = args.complex_pdb,
                                 scaffold_pdb = args.fragment_pdb,
                                 conf_file = args.conf_file,
                                 output = args.output,
                                 chain_id = CHAIN_ID)
    # Terminal fragment hopping
    else:
        if args.covalent:
           from frag_hop.helpers.hopping import run_covalent_frag_replacement
           run_covalent_frag_replacement(complex_pdb = args.complex_pdb,
                                         fragment_pdb = args.fragment_pdb,
                                         conf_file = args.conf_file,
                                         output = args.output,
                                         resname = RESNAME_COV)
        else:
            from frag_hop.helpers.hopping import run_frag_replacement
            run_frag_replacement(complex_pdb = args.complex_pdb,
                                 fragment_pdb = args.fragment_pdb,
                                 conf_file = args.conf_file,
                                 output = args.output,
                                 chain_id = CHAIN_ID)
    logging.info('-' * 75)

if __name__ == '__main__':
    args = parse_args()
    main(args)
