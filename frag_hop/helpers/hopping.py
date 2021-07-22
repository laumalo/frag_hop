"""
This module containts all the workflows to run the fragments or scaffolds
hopping.
"""

import os
import logging

def run_core_replacement(complex_pdb, scaffold_pdb, conf_file,
                         output, chain_id='L'):
    """
    Scaffold replacement protocol for protein-ligand complexes non-covalently
    bound.

    Parameters
    ----------
    complex_pdb : str
        The relative path to the protein-ligand complex PDB file.
    scaffold_pdb : str
        The relative path to the hit scaffold PDB file.
    conf_file : str
        The relative path to the configuration file containing the
        connectivities.
    output : str
        Path to the output folder.
    chain_id : str
        Chain ID of the ligand. Default: L
    """

    # Paths
    PATH_TO_COMPLEX = complex_pdb
    PATH_TO_LIGAND = os.path.join(output, 'LIG.pdb')
    PATH_TO_SCAFFOLD = scaffold_pdb

    OUTPUT_FOLDER = output
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    # Connectivities scaffold-ligand
    from frag_hop.utils import parse_conf_file
    dict_connectivities = parse_conf_file(file_path=conf_file)
    BONDS_SCAFFOLD = dict_connectivities.get(scaffold_pdb.replace(os.getcwd(),
                                                                  '')[1:])
    BONDS_LIGAND = dict_connectivities.get(complex_pdb.replace(os.getcwd(),
                                                                  '')[1:])

    # Extract the target ligand from the complex
    from frag_hop.utils.noncovalent import extract_chain
    logging.info(' - Extracting ligand from chain %s',chain_id)
    extract_chain(chain_id=chain_id,
                  input_file=PATH_TO_COMPLEX,
                  output_file=PATH_TO_LIGAND)

    # Load the target ligand
    logging.info(' - Preparing scaffold.')
    from frag_hop.replacement.scaffold import Scaffold
    target_scaffold = Scaffold(path_ligand=PATH_TO_LIGAND,
                               bonds=BONDS_LIGAND)
    target_scaffold.to_file(OUTPUT_FOLDER, file_name='target_scaffold.pdb')

    # Load the hit scaffold and prepare it for replacement
    hit_scaffold = Scaffold(path_scaffold=PATH_TO_SCAFFOLD,
                            bonds=BONDS_SCAFFOLD)
    hit_scaffold.prepare(target=target_scaffold)
    hit_scaffold.to_file(OUTPUT_FOLDER, file_name='hit_scaffold.pdb')

    # Perform the replacement
    logging.info(' - Replacing selected scaffold.')
    from frag_hop.replacement.replacer import ScaffoldReplacer
    replacer = ScaffoldReplacer(hit_scaffold=hit_scaffold.scaffold,
                                target_ligand=target_scaffold.ligand_prepared)
    replacer.to_file(OUTPUT_FOLDER)

    # Prepare the output complex
    from frag_hop.utils.noncovalent import perform_chain_substitution
    logging.info(' - Creating output protein-ligand complex PDB file.')
    perform_chain_substitution(complex_pdb=PATH_TO_COMPLEX,
                               new_chain_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'new_molecule.pdb'),
                               new_complex_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'new_complex.pdb'),
                               chain_id=chain_id)
    logging.info(' - Output files where created sucessfully under the path:')
    logging.info('    - {}'.format(os.path.join(os.getcwd(),OUTPUT_FOLDER)))


def run_frag_replacement(complex_pdb, fragment_pdb, conf_file,
                         output, chain_id = 'L',
                         SCH_PATH='/opt/schrodinger/suites2020-4'):
    """
    Fragment replacement protocol for protein-ligand complexes with non-covalent
    interactions.

    Given a complex and a hit fragment, the hit fragment is replaced in the
    ligand by breaking and adding a new bond in the specified connectivy points.

    Parameters
    ----------
    complex_pdb : str
        The path to the protein-ligand complex PDB file.
    fragment_pdb : str
        The path to the hit fragment PDB file.
    chain_id : str
        Chain Id where the ligand is.
    connectivity1 : str
        Connection from the ligand to the fragment.
    connectivity2 : str
        onnection from the fragment to the ligand.
    output : str
        Path to the output folder.
    SCH_PATH : str
        Schrodinger’s installation path.
    """

    # Connectivities fragment-ligand
    from frag_hop.utils import parse_conf_file
    dict_connectivities = parse_conf_file(file_path=conf_file)
    BONDS_FRAGMENT = dict_connectivities.get(fragment_pdb)
    BONDS_LIGAND = dict_connectivities.get(complex_pdb)

    OUTPUT_FOLDER = os.path.join(output, 'out_rep')
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    from frag_hop.utils.noncovalent import get_resname_resnum_from_chain
    RESNAME, RESNUM = get_resname_resnum_from_chain(complex_pdb, chain_id)

    logging.info(' - Extracting ligand from chain %s',chain_id)
    from frag_hop.utils.noncovalent import extract_chain
    extract_chain(chain_id=chain_id,
                  input_file=complex_pdb,
                  output_file=os.path.join(OUTPUT_FOLDER, 'LIG_original.pdb'))

    # TODO: Flag to see if schodinger is available
    logging.info(' - Preprocessing ligand with Schrodinger Protein ' +
                 'Preparation Wizzard.')
    from frag_hop.utils.tools import SchrodingerTools
    schrodinger_tools = SchrodingerTools(SCH_PATH=SCH_PATH)
    schrodinger_tools.run_preprocess(folder=OUTPUT_FOLDER,
                                     pdb_in='LIG_original.pdb',
                                     pdb_out='LIG_p.pdb')

    logging.info(' - Preparing fragment.')
    from frag_hop.replacement.fragment import Fragment
    target_fragment = Fragment(path_ligand = os.path.join(OUTPUT_FOLDER,
                                                          'LIG_p.pdb'),
                               bonds = BONDS_LIGAND)
    target_fragment.to_file(OUTPUT_FOLDER)

    # Load the hit scaffold and prepare it for replacement
    hit_fragment = Fragment(path_fragment=fragment_pdb,
                            bonds=BONDS_FRAGMENT)
    hit_fragment.prepare(target=target_fragment)
    hit_fragment.to_file(OUTPUT_FOLDER, file_name='hit_fragment.pdb')

    logging.info(' - Replacing selected fragment.')
    from frag_hop.replacement.replacer import FragmentReplacer
    replacer = FragmentReplacer(hit_fragment=hit_fragment.fragment,
                                target=target_fragment)
    replacer.to_file(OUTPUT_FOLDER, chain_id='L')

    # Replace the ligand chain in the output PDB to create the complex
    logging.info(' - Creating output protein-ligand complex PDB file.')
    from frag_hop.utils.noncovalent import perform_chain_substitution

    perform_chain_substitution(complex_pdb=complex_pdb,
                               new_chain_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'LIG.pdb'),
                               new_complex_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'new_complex.pdb'),
                               chain_id = 'L')
    logging.info(' - Output files where created sucessfully under the path:')
    logging.info('    - {}'.format(os.path.join(os.getcwd(),OUTPUT_FOLDER)))


def run_covalent_frag_replacement(complex_pdb, fragment_pdb, conf_file,
                                  output, resname, chain_id = 'A',
                                  SCH_PATH='/opt/schrodinger/suites2020-4'):
    """
    Fragment replacement protocol for covalent protein-ligand complexes.

    Given a complex and a hit fragment, the hit fragment is replaced in the
    ligand by breaking and adding a new bond in the specified connectivy points.

    Parameters
    ----------
    complex_pdb : str
        The path to the protein-ligand complex PDB file.
    fragment_pdb : str
        The path to the hit fragment PDB file.
    chain_id : str
        Chain Id where the ligand is.
    connectivity1 : str
        Connection from the ligand to the fragment.
    connectivity2 : str
        onnection from the fragment to the ligand.
    output : str
        Path to the output folder.
    SCH_PATH : str
        Schrodinger’s installation path.
    """

    # Connectivities fragment-ligand
    from frag_hop.utils import parse_conf_file
    dict_connectivities = parse_conf_file(file_path=conf_file)

    BONDS_FRAGMENT = dict_connectivities.get(fragment_pdb.replace(os.getcwd(),
                                                                  '')[1:])
    BONDS_LIGAND = dict_connectivities.get(complex_pdb.replace(os.getcwd(),
                                                               '')[1:])

    OUTPUT_FOLDER = os.path.join(output, 'out_rep')
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    from frag_hop.utils.covalent import get_resnum_from_resname
    RESNUM = get_resnum_from_resname(complex_pdb, resname)

    # Extract ligand residue
    logging.info(' - Extracting ligand residue')
    from frag_hop.utils.covalent import extract_residue
    extract_residue(input_pdb=complex_pdb,
                    resname=resname,
                    output_pdb=os.path.join(OUTPUT_FOLDER, 'RES.pdb'))

    from frag_hop.utils.tools import SchrodingerTools
    logging.info(' - Preprocessing ligand with Schrodinger Protein ' +
                 'Preparation Wizzard.')
    schrodinger_tools = SchrodingerTools(SCH_PATH=SCH_PATH)
    schrodinger_tools.run_preprocess(folder=OUTPUT_FOLDER,
                                     pdb_in='RES.pdb',
                                     pdb_out='RES_p.pdb')

    logging.info(' - Preparing fragment.')
    from frag_hop.replacement import Fragment
    target_fragment = Fragment(path_ligand = os.path.join(OUTPUT_FOLDER,
                                                          'RES_p.pdb'),
                               bonds = BONDS_LIGAND,
                               resname = resname)
    target_fragment.to_file(OUTPUT_FOLDER)

    # Load the hit scaffold and prepare it for replacement
    hit_fragment = Fragment(path_fragment=fragment_pdb,
                            bonds=BONDS_FRAGMENT)
    hit_fragment.prepare(target=target_fragment)
    hit_fragment.to_file(OUTPUT_FOLDER, file_name='hit_fragment.pdb')

    logging.info(' - Replacing selected fragment.')
    from frag_hop.replacement.replacer import FragmentReplacer
    replacer = FragmentReplacer(hit_fragment=hit_fragment.fragment,
                                target=target_fragment)
    replacer.to_file(OUTPUT_FOLDER, resname = resname, resnum = RESNUM)

    # Create output protein-ligand complex
    logging.info(' - Creating output protein-ligand complex PDB file.')
    from frag_hop.utils.covalent import perform_residue_substitution
    perform_residue_substitution(complex_pdb=complex_pdb,
                                 new_residue_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'new_molecule.pdb'),
                                 new_complex_pdb=os.path.join(OUTPUT_FOLDER,
                                                        'new_complex.pdb'),
                                 resname=resname)
    logging.info(' - Output files where created sucessfully under the path:')
    logging.info('    - {}'.format(os.path.join(os.getcwd(),OUTPUT_FOLDER)))


class HoppingSelector(object):
    """ It defines a hopping selector. """

    def __init__(self, core, covalent):
        """
        Given if it is a core fragment or not and if the ligand is covalently
        bound of not to the protein, it selects the corresponding hopping
        method.

        Parameters
        ----------
        core : bool
            It indicates if the fragment is a scaffold or not.
        covalent : bool
             It indicates if the ligand is covalently bound to the protein
              or not.
        """
        self.core = core
        self.covalent = covalent

    def run(self, complex_pdb, hit_pdb, conf_file, id_lig, output):
        """
        It runs the corresponding hopping method.

        Parameters
        ----------
        complex_pdb : str
            The path to the protein-ligand complex PDB file.
        hit_pdb : str
            The path to the hit fragment/scaffold PDB file.
        conf_file : str
            Path to the configuration file.
        id_lig : str
            Chain Id/resname where the ligand is located.
        output : str
            Path to the output folder.
        """
        # Scaffold hopping
        if self.core:
            if self.covalent:
                raise NotImplementedError
            else:
                run_core_replacement(complex_pdb = complex_pdb,
                                     scaffold_pdb = hit_pdb,
                                     conf_file = conf_file,
                                     output = output,
                                     chain_id = id_lig)
        # Terminal fragment hopping
        else:
            if self.covalent:
               run_covalent_frag_replacement(complex_pdb = complex_pdb,
                                             fragment_pdb = hit_pdb,
                                             conf_file = conf_file,
                                             output = output,
                                             resname = id_lig)
            else:
                run_frag_replacement(complex_pdb = complex_pdb,
                                     fragment_pdb = hit_pdb,
                                     conf_file = conf_file,
                                     output = output,
                                     chain_id = id_lig)
