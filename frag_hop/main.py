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
import os

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
    parser = ap.ArgumentParser(description="FRAGHOP: PERFORM FRAGMENT  " +
                               "REPLACEMENT IN A PROTEIN-LIGAND COMPLEX")
    parser.add_argument("complex_pdb", type=str,
                        help="Path to the protein-ligand complex.")

    parser.add_argument("fragment_pdb", type=str,
                        help="Path to the hit fragment.")

    parser.add_argument("-c", "--conf_file", type=str,
                        help="Path to the configuration file.")

    parser.add_argument("-o", "--output", type=str, default='out',
                        help="Path to the output folder. Default: out.")

    parser.add_argument("-p", "--pele", type=str, default=None,
                        help="Complete path to PELE serial. Default: None")

    parser.add_argument("-ct", "--control_file", type=str, default=None,
                        help="Path to control file for the PELE exploration.")

    parser.add_argument("-s", "--screening", type=str, default='basic',
                        help="Screening method.")

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

    >>> python -m frag_hop.main data/TestingFiles/complexes/test_covalent.pdb data/TestingFiles/fragments/frag1.pdb -c data/TestingFiles/test_configurations.conf --covalent

    >>> python -m frag_hop.main data/TestingFiles/complexes/test_noncovalent.pdb data/TestingFiles/fragments/frag2.pdb -c data/TestingFiles/test_configurations.conf

    >>> python -m frag_hop.main data/TestingFiles/complexes/test_scaffold.pdb data/TestingFiles/scaffolds/scaffold_2FJP.pdb -c data/TestingFiles/test_configurations.conf --core
    """

    # Type of hopping
    type_rep = 'scaffold' if args.core else 'terminal'
    type_ligand = 'covalent' if args.covalent else 'non-covalent'

    is_screening = os.path.isdir(args.fragment_pdb)

    from frag_hop.utils import parse_conf_file
    dict_connectivities = parse_conf_file(file_path=args.conf_file)
    BONDS_LIG = dict_connectivities.get(args.complex_pdb)
    if not is_screening:
        BONDS_FRAG = dict_connectivities.get(args.fragment_pdb)

    logging.info('-' * 75)
    logging.info('Fragment hopping for a protein-ligand complex')
    logging.info('-' * 75)
    logging.info(' - Input information:')
    logging.info('    - Protein-ligand complex PDB: %s', args.complex_pdb)
    if not is_screening:
        logging.info('    - Hit fragment PDB: %s', args.fragment_pdb)
    else:
        logging.info('    - Fragment library: %s', args.fragment_pdb)
    logging.info('    - Type of fragment replacement: %s', type_rep)
    logging.info('    - Type of ligand: %s', type_ligand)
    logging.info('    - Connectivity remaining ligand: %s', BONDS_LIG[0])
    if not is_screening:
        logging.info('    - Connectivity fragment: %s', BONDS_FRAG[0])
        logging.info(' - You are going to perform a single fragment hopping.')
    else:
        logging.info(' - You are going to perform a library screening.')
    logging.info('-' * 75)

    # Library screening
    if is_screening:
        from frag_hop.helpers.screening import Screening
        screen_library = Screening(complex_pdb=args.complex_pdb,
                                   lib_path=args.fragment_pdb,
                                   conf_file=args.conf_file)
        screen_library.run(type_screening=args.screening,
                           CONTROL_FILE=args.control_file,
                           PELE_exc=args.pele)

    # Simple hopping
    else:
        from frag_hop.helpers.hopping import HoppingSelector
        hopping_method = HoppingSelector(core=args.core,
                                         covalent=args.covalent)
        hopping_method.run(complex_pdb=args.complex_pdb,
                           hit_pdb=args.fragments_pdb,
                           conf_file=args.conf_file,
                           id_lig='GRW' if args.covalent else 'L',
                           output=args.output)
        logging.info('-' * 75)


if __name__ == '__main__':
    args = parse_args()
    main(args)
