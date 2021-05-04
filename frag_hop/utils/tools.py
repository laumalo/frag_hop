import collections
import sys
import os
import numpy as np
import time


class PDBTools():
    """
    Class that contains all methods applied to handle PDB files.
    """

    def pbd_uniqname(self, pdb_path):
        """
        It handles the atom names of a PDB file and modifies them to have unique
        atom names.
        Implementation based on: https://github.com/haddocking/pdb-tools

        Parameters
        ----------
        pdb_path : str
            Path to the PDB file to rename the atom names.
        """
        def rename_atoms(fhandle):
            """
            Renames HETATM atoms on each residue based on their element.
            """
            prev_res = None
            for line_idx, line in enumerate(fhandle):
                if line.startswith('HETATM'):
                    element = line[76:78].strip()
                    if not element:
                        emsg = \
                            'ERROR!! No element found in line {}'.format(
                                line_idx)
                        sys.stderr.write(emsg)
                        sys.exit(1)
                    resuid = line[17:27]
                    if prev_res != resuid:
                        prev_res = resuid
                        element_idx = collections.defaultdict(lambda: 1)
                    spacer = ' ' if len(element) == 1 else ''
                    name = (spacer + element +
                            str(element_idx[element])).ljust(4)
                    line = line[:12] + name + line[16:]
                    element_idx[element] += 1
                yield line

        pdbfh = open(pdb_path, 'r')
        new_pdb = rename_atoms(pdbfh)
        os.remove(pdb_path)
        outf = open(pdb_path, 'w')
        try:
            for lineno, line in enumerate(new_pdb):
                outf.write(line)
        except IOError:
            pass
        pdbfh.close()

    def rename_atoms_fragment(self, pdb_path):
        """
        It handles the atom names of a PDB file and modifies them to have unique
        atom names. It only modifies the atom names of the ligand, if the CYS is
        also in the residue, their atom names remain the same.

        Implementation based on: https://github.com/haddocking/pdb-tools

        Parameters
        ----------
        pdb_path : str
            Path to the PDB file to rename the atom names.
        """
        def rename_atoms(fhandle):
            """
            Renames HETATM or ATOM atoms on each residue based on their element.
            """
            prev_res = None
            for line_idx, line in enumerate(fhandle):

                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_name = [x for x in line.split(' ') if not x == ''][2]
                    if any(map(str.isdigit, atom_name)) and not 'HB' in line:
                        element = line[76:78].strip()
                        if not element:
                            emsg = \
                                'ERROR!! No element found in line {}'.format(
                                    line_idx)
                            sys.stderr.write(emsg)

                        resuid = line[17:27]
                        if prev_res != resuid:
                            prev_res = resuid
                            element_idx = collections.defaultdict(lambda: 1)
                        spacer = ' ' if len(element) == 1 else ''
                        name = (spacer + element +
                                str(element_idx[element])).ljust(4)
                        line = line[:12] + name + line[16:]
                        element_idx[element] += 1
                yield line

        pdbfh = open(pdb_path, 'r')
        new_pdb = rename_atoms(pdbfh)
        os.remove(pdb_path)
        outf = open(pdb_path, 'w')
        try:
            for lineno, line in enumerate(new_pdb):
                outf.write(line)
        except IOError:
            pass
        pdbfh.close()


class RDKitTools():
    """
    Class that contains all the methods using RDKit.
    """
    def get_atomid_by_atomname(self, mol, atom_name):
        """
        It returns the atom ID of a molecule given its atom name.

        Parameters
        ----------
        mol : rdkit.molecule object
            Molecule.
        atom_name : str
            Atom name.

        Returns
        -------
        atom_id : str
            Atom ID.
        """
        atom_id = [atom.GetIdx() for atom in mol.GetAtoms()
                   if atom.GetPDBResidueInfo().GetName().strip() == atom_name][0]
        return atom_id

    def get_element_by_atomname(self, mol, atom_name):
        """
        It returns the element symbol of a molecule given its atom name.

        Parameters
        ----------
        mol : rdkit.molecule object
            Molecule.
        atom_name : str
            Atom name.

        Returns
        -------
        element : str
            Element symbol.
        """
        element = [atom.GetSymbol() for atom in mol.GetAtoms()
                   if atom.GetPDBResidueInfo().GetName().strip() == atom_name][0]
        return element


def perform_residue_substitution(complex_pdb, new_residue_pdb, new_complex_pdb,
                                 resname):
    """
    Given a PDB with a complex and the PDB of the residue that has to be
    replaced. It replaces from the original complex the residue and exports
    the PDB of the new strucutre.

    Parameters
    ----------
    complex_pdb : str
        Path to the initial complex PDB.
    new_residue_pdb : str
        Path to the residue PDB you want to modidy in the complex.
    new_complex_pdb : str
        Path to the output complex PDB with the residue modified.
    resname : str
        Residue name for the modified residue.
    """

    f1 = open(complex_pdb, 'r')
    f2 = open(new_residue_pdb, 'r')
    f3 = open(new_complex_pdb, 'w')
    lines_to_delete, new_lines = ([] for x in range(2))

    # Read original complex
    lines = f1.readlines()
    lines_to_delete = [line_num for line_num, line in enumerate(lines)
                       if resname in line]

    # Load new residue atoms and correct atom numbers format
    new_lines = [line for line in f2.readlines() if 'ATOM' in line]
    new_lines = [line[:6] + str(int(line[6:11])).rjust(5) + line[11:]
                 for line in new_lines]

    # Generate output file for the PDB structure
    new_file_lines = \
        lines[:lines_to_delete[0]] + new_lines + \
        lines[lines_to_delete[-1] + 1:]
    for line in new_file_lines:
        f3.write(line)

class SchrodingerTools():
    def __init__(self, SCH_PATH):
        """
        Parameters
        ----------
        SCH_PATH : str
            Schrodinger’s installation path.
        self.sch_path = SCH_PATH
        """
        self.sch_path = SCH_PATH

    def run_preprocess(self, folder, pdb_in, pdb_out):
        """
        It preprocess a PDB file using Schrodinger Protein Preparation Wizzard.

        Parameters
        ----------

        folder : str
            Path to the folder containing the PDB file.
        pdb_in : str
            Relative path to the input PDB.
        pdb_out :str
            Relative path to the output PDB.
        """

        curr_dir = os.getcwd()
        command = '{}/utilities/prepwizard {} {}'.format(self.sch_path, pdb_in,
                                                         pdb_out)\
                    + ' -noepik -noccd -noimpref -nohtreat'

        # Run command to preprocess the PDB
        os.chdir(folder)
        os.system(command)

        # Wait until the output file is created
        condition = os.path.isfile(pdb_out)
        while not condition:
            condition = os.path.isfile(pdb_out)
            time.sleep(0.1)

        # Remove log file
        os.remove(pdb_in.replace('.pdb', '.log'))
        os.chdir(curr_dir)


