from replace_fragments import PrepareSystem, Replacer
from utils import perform_complex_mutation, extract_residue
import os
import argparse as ap


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
    parser = ap.ArgumentParser(description="FRAGMENTS EXTRACTION")
    parser.add_argument("complex_pdb", type=str,
                        help="The path to the complex PDB.")

    parser.add_argument("fragment_original_pdb", type=str,
                        help="The path to thefragment PDB.")

    parser.add_argument("-c1", "--connectivity1", type=str,
                        help="Connection from the ligand to the fragment.")

    parser.add_argument("-c2", "--connectivity2", type=str,
                        help="Connection from the fragment to the ligand.")

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
    >>> python Replacement/main.py TestingFiles/covalent_ligand/covalent_covid_scaffold.pdb TestingFiles/covalent_ligand/S1p_6367.pdb
    -c1 C13-N7 -c2 C1-H4
    """
    bond_atoms = [args.connectivity2.split('-'), args.connectivity1.split('-')]

    os.makedirs('out', exist_ok=True)
    extract_residue(pdb_file=args.complex_pdb,
                    resname='GRW')

    PrepareSystem(initial_complex=args.complex_pdb,
                  fragment=args.fragment_original_pdb,
                  bond_atoms=bond_atoms)

    Replacer(ligand_pdb='out/RES.pdb',
             fragment_pdb='out/frag_prepared.pdb',
             bond_atoms=bond_atoms)

    perform_complex_mutation(complex_pdb=args.complex_pdb,
                             new_residue_pdb='out/merged.pdb',
                             resname='GRW')


if __name__ == '__main__':
    args = parse_args()
    main(args)
