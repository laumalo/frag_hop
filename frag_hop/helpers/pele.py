"""
This module allows to prepare the system to run a PELE simulation after the
replacement protocol.
"""

import os
import logging


def prepare_folder(FOLDER_TO_PREPARE, conf_file=None, run_file=None):
    """
    It prepares an output folder of FragHop to run a PELE simulation by
    generating the parameters templates and copying the PELE configuration and
    run file in it.

    Parameters
    ----------
    FOLDER_TO_PREPARE : str
        Path to the folder that to prepare.
    conf_file : str
        Path to the configuration file to use.
    run_file : str
        Path to the PELE run file to use.
    """
    def generate_parameters(path, output_path,
                            forcefield='OPLS2005',
                            charge_method='am1bcc'):
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
    os.makedirs(os.path.join(FOLDER_TO_PREPARE, 'lig'), exist_ok=True)
    shutil.copy(os.path.join(FOLDER_TO_PREPARE,
                             'out_rep/complex_merged.pdb'), FOLDER_TO_PREPARE)
    shutil.copy(os.path.join(FOLDER_TO_PREPARE, 'out_rep/merged.pdb'),
                os.path.join(FOLDER_TO_PREPARE, 'lig/LIG.pdb'))

    # Templatizate ligand using Peleffy
    generate_parameters(path=os.path.join(FOLDER_TO_PREPARE, 'lig/LIG.pdb'),
                        output_path=FOLDER_TO_PREPARE)

    if conf_file is not None:
        shutil.copy(conf_file, os.path.join(FOLDER_TO_PREPARE, 'pele.conf'))
    else:
        logging.info('  - No configuration file was provided. It has to be put' +
                     ' under the folder path before running the PELE simulation.')

    if run_file is not None:
        shutil.copu(run_file, os.path.join(FOLDER_TO_PREPARE, 'run.sh'))

def run_PELE_Job(folder, PELE_exc=None, num_proc=None, run=False):

    os.chdir(folder)
    if not os.path.isfile(folder, 'run.sh'):
        path = os.path.abspath('data/parameters/run_example.sh')
        new_file = open('run.sh', 'w')
        with open(path, 'r') as f:
            data = f.readlines()
            for line in data:
                line = line.replace('PELE_exc', PELE_exc) \
                    if 'PELE_exc' in line else line
                line = line.replace('num_proc', str(num_proc)) \
                    if 'num_proc' in line else line
                new_file.write(line)
    if run:
        os.system('bash run.sh')

