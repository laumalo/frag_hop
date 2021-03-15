# General imports
import os
import argparse as ap

# Local imports
from replace_fragments import PrepareFragment, Replacer
from utils import perform_residue_substitution, extract_residue


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
                               "FRAGHOP: PERFORM FRAGMENT REPLACEMENT")
    parser.add_argument("complex_pdb", type=str,
                        help="The path to the complex PDB.")

    parser.add_argument("fragment_pdb", type=str,
                        help="The path to the fragment PDB.")

    parser.add_argument("-c1", "--connectivity1", type=str,
                        help="Connection from the ligand to the fragment.")

    parser.add_argument("-c2", "--connectivity2", type=str,
                        help="Connection from the fragment to the ligand.")

    parser.add_argument("-o", "--output", type=str, default='out',
                        help="Path to the output folder. Default: out.")

    parsed_args = parser.parse_args()

    return parsed_args


def main(args):
    """
    It reads the command-line arguments and runs the fragment replacement
    protocol for a covalent ligand.

    Parameters
    ----------
    args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user

    Examples
    --------
    From the command-line:
    >>> python Replacement/main.py TestingFiles/covalent_ligand/covalent_covid_scaffold.pdb TestingFiles/covalent_ligand/S1p_6367.pdb -c1 C13-N7 -c2 C1-H4
    """

    OUTPUT_FOLDER = os.path.join(args.output, 'out_rep')
    RESNAME = 'GRW'
    BOND_ATOMS = [args.connectivity2.split('-'), args.connectivity1.split('-')]

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    extract_residue(input_pdb=args.complex_pdb,
                    resname=RESNAME,
                    output_pdb=os.path.join(OUTPUT_FOLDER, 'RES.pdb'))

    PrepareFragment(initial_complex=args.complex_pdb,
                    fragment=args.fragment_pdb,
                    bond_atoms=BOND_ATOMS,
                    out_folder=OUTPUT_FOLDER,
                    resname=RESNAME)

    Replacer(ligand_pdb=os.path.join(OUTPUT_FOLDER, 'RES.pdb'),
             fragment_pdb=os.path.join(OUTPUT_FOLDER, 'frag_prepared.pdb'),
             bond_atoms=BOND_ATOMS,
             out_folder=OUTPUT_FOLDER)

    perform_residue_substitution(complex_pdb=args.complex_pdb,
                                 new_residue_pdb=os.path.join(OUTPUT_FOLDER,
                                                              'merged.pdb'),
                                 new_complex_pdb=os.path.join(OUTPUT_FOLDER,
                                                              'complex_merged.pdb'),
                                 resname=RESNAME)


if __name__ == '__main__':
    args = parse_args()
    main(args)
