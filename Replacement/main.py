# General imports
import os
import argparse as ap

# Local imports
from replace_fragments import PrepareFragment, Replacer
from utils import perform_residue_substitution, extract_residue, run_preprocess


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
                        help="Path to the complex PDB.")

    parser.add_argument("fragment_pdb", type=str,
                        help="Path to the fragment PDB.")

    parser.add_argument("-c1", "--connectivity1", type=str,
                        help="Connection from the ligand to the fragment.")

    parser.add_argument("-c2", "--connectivity2", type=str,
                        help="Connection from the fragment to the ligand.")

    parser.add_argument("-o", "--output", type=str, default='out',
                        help="Path to the output folder. Default: out.")

    parser.add_argument("--covalent",
                        dest="covalent",
                        action='store_true',
                        help="Perform the fragment replacement in a covalent" +
                        " complex.")

    parser.set_defaults(covalent=False)

    parsed_args = parser.parse_args()

    return parsed_args

def run_covalent_replacement(complex_pdb, fragment_pdb, resname,connectivity1,
                             connectivity2, output,
                             SCH_PATH='/opt/schrodinger/suites2020-4'):
    """
    Fragment replacement protocol for covalently bound ligands. Given a complex
    and a fragment, the new fragment is merged with the scaffold by breaking and
    adding a new bond in the specified connectivy atoms of both scaffold and
    fragment.

    Parameters
    ----------
    complex_pdb : str
        The path to the complex PDB.
    fragment_pdb : str
        The path to the fragment PDB.
    connectivity1 : str
        Connection from the ligand to the fragment.
    connectivity2 : str
        onnection from the fragment to the ligand.
    output : str
        Path to the output folder.
    SCH_PATH : str
        Schrodinger’s installation path.
    """

    OUTPUT_FOLDER = os.path.join(output, 'out_rep')
    BOND_ATOMS = [connectivity2.split('-'), connectivity1.split('-')]

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    extract_residue(input_pdb=complex_pdb,
                    resname=resname,
                    output_pdb=os.path.join(OUTPUT_FOLDER, 'RES.pdb'))

    run_preprocess(SCH_PATH=SCH_PATH,
                   folder=OUTPUT_FOLDER,
                   pdb_in='RES.pdb',
                   pdb_out='RES_p.pdb')

    PrepareFragment(initial_complex=complex_pdb,
                    fragment=fragment_pdb,
                    bond_atoms=BOND_ATOMS,
                    out_folder=OUTPUT_FOLDER,
                    resname=resname)

    Replacer(ligand_pdb=os.path.join(OUTPUT_FOLDER, 'RES_p.pdb'),
             fragment_pdb=os.path.join(OUTPUT_FOLDER, 'frag_prepared.pdb'),
             ref_fragment_pdb=fragment_pdb,
             bond_atoms=BOND_ATOMS,
             out_folder=OUTPUT_FOLDER)

    perform_residue_substitution(complex_pdb=complex_pdb,
                                 new_residue_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'merged.pdb'),
                                 new_complex_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'complex_merged.pdb'),
                                 resname=resname)
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
    >>> python replacement/main.py TestingFiles/covalent_ligand/covalent_covid_scaffold.pdb TestingFiles/covalent_ligand/S1p_6367.pdb -c1 C13-N7 -c2 C1-H4
    """
    RESNAME = 'GRW'
    if args.covalent:
        run_covalent_replacement(complex_pdb = args.complex_pdb,
                                 fragment_pdb = args.fragment_pdb,
                                 resname = RESNAME,
                                 connectivity1 = args.connectivity1,
                                 connectivity2 = args.connectivity2,
                                 output = args.output)
    else:
        print('Non covalent version not fully implemented yet.')

if __name__ == '__main__':
    args = parse_args()
    main(args)
