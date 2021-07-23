"""
This module allows to prepare the system to run a PELE simulation after the
replacement protocol.
"""

import os
import logging


def generate_templates(FOLDER_TO_PREPARE):
    """
    It prepares an output folder of FragHop to run a PELE simulation by
    generating the parameters templates and copying the PELE configuration and
    run file in it.

    Parameters
    ----------
    FOLDER_TO_PREPARE : str
        Path to the folder that to prepare.
    """
    def generate_parameters(path, output_path,
                            forcefield='OPLS2005',
                            charge_method='gasteiger'):
        """
        It generates the parameters of the molecule (from the input_file)
        as DataLocal in the output folder.
        Parameters
        ----------
        smiles : str
            The smiles tag representing the molecule to minimize
        mol_id : str
            Unique id to identify the molecule to minimize
        output_path : str
            The output path where parameters will be saved
        forcefield : str
            The Open Force Field force field to generate the parameters
            with
        charge_method : str
            The charge method to calculate the partial charges with
        """

        # Create representation of a particular molecule
        from peleffy.topology import Molecule
        molecule = Molecule(path)

        from peleffy.forcefield import ForceFieldSelector
        selector = ForceFieldSelector()
        ff = selector.get_by_name(forcefield)
        parameters = ff.parameterize(molecule, charge_method=charge_method)

        from peleffy.topology import Topology
        topology = Topology(molecule, parameters)

        # Paths you may need
        from peleffy.utils import OutputPathHandler
        output_handler = OutputPathHandler(molecule, ff,
                                           output_path=output_path,
                                           as_datalocal=True)
        rotamer_library_path = output_handler.get_rotamer_library_path()
        impact_template_path = output_handler.get_impact_template_path()

        # Generate its rotamer library
        from peleffy.topology import RotamerLibrary
        rotamer_library = RotamerLibrary(molecule)
        rotamer_library.to_file(rotamer_library_path)

        # Generate its parameters and template file
        from peleffy.template import Impact
        impact = Impact(topology)
        impact.to_file(impact_template_path)

    import shutil

    # Copy complex and residue PDB files
    logging.info('  -   Preparinf folder for a PELE simulation.')
    os.makedirs(os.path.join(FOLDER_TO_PREPARE, 'lig'), exist_ok=True)
    shutil.copy(os.path.join(FOLDER_TO_PREPARE,
                             'out_rep/new_complex.pdb'), FOLDER_TO_PREPARE)
    shutil.copy(os.path.join(FOLDER_TO_PREPARE, 'out_rep/LIG.pdb'),
                os.path.join(FOLDER_TO_PREPARE, 'lig/LIG.pdb'))

    # Templatizate ligand using Peleffy
    generate_parameters(path=os.path.join(FOLDER_TO_PREPARE, 'lig/LIG.pdb'),
                        output_path=FOLDER_TO_PREPARE)

def PELErunner(control_file, PELE_exc, num_proc=16):
    """
    It runs a PELE simulation with the protocol described in the fetched control
    file.

    Parameters
    ----------
    control_file : str
        Path to the control file.
    PELE_exec : str
        Path to the PELE executable
    PELE_src : str
        Path to PELE source folder
    """
    import subprocess

    logging.info('  -   Starting PELE simulation.')
    cmd =  "mpirun -n {} {} {}".format(num_proc, PELE_exc, control_file)
    subprocess.call(cmd.split())

def FRAGrunner(control_file, frag_conf_file, PELE_exc):
    """
    It runs a FragPELE simulation with the protocol described in the fetched
    control file to grow a fragment.
    raise NotImplementedError()

