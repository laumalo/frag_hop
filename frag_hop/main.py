"""
This module is designed to run FragHop through the command-line.
"""

__author__ = "Laura Malo"
__license__ = "GPL"
__maintainer__ = "Laura Malo"
__email__ = "laura.maloroset@bsc.es"


# General imports
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
    parser = ap.ArgumentParser(description=
                               "FRAGHOP: PERFORM FRAGMENT REPLACEMENT IN A " +
                               "PROTEIN-LIGAND COMPLEX")
    parser.add_argument("complex_pdb", type=str,
                        help="Path to the protein-ligand complex.")

    parser.add_argument("fragment_pdb", type=str,
                        help="Path to the hit fragment.")

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

def run_covalent_replacement(complex_pdb, fragment_pdb, resname,
                             connectivity1, connectivity2, output,
                             SCH_PATH='/opt/schrodinger/suites2020-4'):
    """
    Fragment replacement protocol for protein-ligand complexes covalently bound.

    Given a complex and a hit fragment, the hit fragment is replaced in the
    ligand by breaking and adding a new bond in the specified connectivy points.

    Parameters
    ----------
    complex_pdb : str
        The path to the protein-ligand complex.
    fragment_pdb : str
        The path to the hit fragment.
    connectivity1 : str
        Connection from the ligand to the fragment.
    connectivity2 : str
        Connection from the fragment to the ligand.
    output : str
        Path to the output folder.
    SCH_PATH : str
        Schrodinger’s installation path.
    """

    OUTPUT_FOLDER = os.path.join(output, 'out_rep')
    BOND_ATOMS = [connectivity2.split('-'), connectivity1.split('-')]

    from frag_hop.utils.covalent import get_resnum_from_resname
    RESNUM = get_resnum_from_resname(complex_pdb, resname)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    from frag_hop.utils.covalent import extract_residue
    extract_residue(input_pdb=complex_pdb,
                    resname=resname,
                    output_pdb=os.path.join(OUTPUT_FOLDER, 'RES.pdb'))

    from frag_hop.utils.tools import SchrodingerTools
    schrodinger_tools = SchrodingerTools(SCH_PATH=SCH_PATH)
    schrodinger_tools.run_preprocess(folder=OUTPUT_FOLDER,
                                     pdb_in='RES.pdb',
                                     pdb_out='RES_p.pdb')

    from frag_hop.replacement import Fragment
    fragment = Fragment(initial_complex=complex_pdb,
                        fragment=fragment_pdb,
                        bond_atoms=BOND_ATOMS,
                        resname=resname)
    fragment.to_file(OUTPUT_FOLDER)

    from frag_hop.replacement import Replacer
    Replacer(ligand_pdb=os.path.join(OUTPUT_FOLDER, 'RES_p.pdb'),
             fragment_pdb=os.path.join(OUTPUT_FOLDER, 'frag_prepared.pdb'),
             ref_fragment_pdb=fragment_pdb,
             resname=resname,
             resnum=RESNUM,
             bond_atoms=BOND_ATOMS,
             out_folder=OUTPUT_FOLDER)

    from frag_hop.utils.covalent import perform_residue_substitution
    perform_residue_substitution(complex_pdb=complex_pdb,
                                 new_residue_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'merged.pdb'),
                                 new_complex_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'complex_merged.pdb'),
                                 resname=resname)

def run_replacement(complex_pdb, fragment_pdb, chain_id,
                    connectivity1, connectivity2, output,
                    SCH_PATH='/opt/schrodinger/suites2020-4'):
    """
    Fragment replacement protocol for protein-ligand complexes with non-covalent
    interactions.

    Given a complex and a hit fragment, the hit fragment is replaced in the
    ligand by breaking and adding a new bond in the specified connectivy points.

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

    from frag_hop.utils.noncovalent import get_resname_resnum
    RESNAME, RESNUM = get_resname_resnum(complex_pdb, chain_id)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    from frag_hop.utils.noncovalent import extract_chain
    extract_chain(chain_id=chain_id,
                  input_file=complex_pdb,
                  output_file=os.path.join(OUTPUT_FOLDER, 'LIG_original.pdb'))

    from frag_hop.utils.tools import SchrodingerTools
    schrodinger_tools = SchrodingerTools(SCH_PATH=SCH_PATH)
    schrodinger_tools.run_preprocess(folder=OUTPUT_FOLDER,
                                     pdb_in='LIG_original.pdb',
                                     pdb_out='LIG_p.pdb')

    from frag_hop.replacement import Fragment
    fragment = Fragment(initial_complex=complex_pdb,
                        fragment=fragment_pdb,
                        bond_atoms=BOND_ATOMS,
                        resname=RESNAME)
    fragment.to_file(OUTPUT_FOLDER)

    from frag_hop.replacement import Replacer
    Replacer(ligand_pdb=os.path.join(OUTPUT_FOLDER, 'LIG_p.pdb'),
             fragment_pdb=os.path.join(OUTPUT_FOLDER, 'frag_prepared.pdb'),
             ref_fragment_pdb=fragment_pdb,
             resname=RESNAME,
             resnum=RESNUM,
             bond_atoms=BOND_ATOMS,
             out_folder=OUTPUT_FOLDER)

    # Replace the ligand chain in the output PDB to create the complex
    from frag_hop.utils.noncovalent import perform_chain_substitution

    perform_chain_substitution(complex_pdb=complex_pdb,
                               new_chain_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'LIG.pdb'),
                               new_complex_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'complex_merged.pdb'),
                               chain_id = 'L')

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

    >>> python replacement/main.py TestingFiles/covalent_ligand/covalent_covid_scaffold.pdb TestingFiles/covalent_ligand/S1p_6367.pdb -c1 C13-N7 -c2 C1-H4

    """
    CHAIN_ID = 'L'
    RESNAME_COV = 'GRW' # Default parameters for covalent ligands

    if args.covalent:
        run_covalent_replacement(complex_pdb = args.complex_pdb,
                                 fragment_pdb = args.fragment_pdb,
                                 resname = RESNAME_COV,
                                 connectivity1 = args.connectivity1,
                                 connectivity2 = args.connectivity2,
                                 output = args.output)
    else:
        run_replacement(complex_pdb = args.complex_pdb,
                        fragment_pdb = args.fragment_pdb,
                        chain_id = CHAIN_ID,
                        connectivity1 = args.connectivity1,
                        connectivity2 = args.connectivity2,
                        output = args.output)

if __name__ == '__main__':
    args = parse_args()
    main(args)
