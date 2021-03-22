import collections
import sys
import os
import numpy as np


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

    def extract_chain(self, chain_id, input_file, output_file):
        """
        It extracts a chain from a PDB file and generates a new PDB file only
        containing this chain.

        Parameters
        ----------
        chain_id : str
            Chain ID.
        input_file : str
            Path to the PDB file to extract the chain from.
        output_file : str
            Path to the new PDB file to save the extracted chain.
        """
        def select_chain(fhandle, chain_set):
            """
            Filters the PDB file for specific chain identifiers.
            """

            records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
            for line in fhandle:
                if line.startswith(records):
                    if line[21] not in chain_set:
                        continue
                yield line

        pdbfh = open(input_file, 'r')
        new_pdb = select_chain(pdbfh, chain_id)
        outf = open(output_file, 'w')
        try:
            for lineno, line in enumerate(new_pdb):
                outf.write(line)
        except IOError:
            pass


class MathTools():
    """
    Class that contains all methods related to mathematical operations on points
    and vectors.
    """

    def rotation_matrix_from_vectors(self, vec1, vec2):
        """
        Finds the rotation matrix that aligns vec1 to vec2
        Parameters
        ----------
        vec1: np.array
            A 3d "source" vector
        vec2: np.array
            A 3d "destination" vector

        Returns
        -------
        rotation_matrix: np.array
            A transform matrix (3x3) which when applied to vec1,
            aligns it with vec2.
        """
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), \
               (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + \
            kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix

    def rotation_matrix_axis(self, radi, u):
        """
        Finds the general rotation matrix of a certain angle along an axis
        defined by a vector.

        Parameters
        ----------
        radi : float
            Rotation angle, in radians.
        u : np.array
            A 3d vector defining the rotation axis.

        Returns
        -------
        rotation_matrix: np.array
            A transform matrix (3x3) which when applied to a point, it rotates
            a certain angle along an axis.
        """
        ux, uy, uz = u
        rotation_matrix =  np.array([[np.cos(radi) + ux**2 * (1 - np.cos(radi)),
                          ux * uy * (1 - np.cos(radi)) - uz * np.sin(radi),
                          ux * uz * (1 - np.cos(radi)) + uy * np.sin(radi)],
                         [uy * ux * (1 - np.cos(radi)) + uz * np.sin(radi),
                             np.cos(radi) + uy**2 * (1 - np.cos(radi)),
                          uy * uz * (1 - np.cos(radi)) - ux * np.sin(radi)],
                         [uz * ux * (1 - np.cos(radi)) - uy * np.sin(radi),
                          uz * uy * (1 - np.cos(radi)) + ux * np.sin(radi),
                             np.cos(radi) + uz**2 * (1 - np.cos(radi))]],
                        dtype=np.double)
        return rotation_matrix

def distance(x, y):
    """
    It computes the distance between two points.

    Parameters
    ----------
    x : np.array
        Point
    y : np.array
        Point

    Returns
    -------
    dist : float
        Distance between the two points
    """
    import math
    return math.sqrt((x[0] - y[0]) ** 2 +
                     (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2)

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



def extract_residue(input_pdb, resname, output_pdb):
    """
    It extracts one residue of a PDB file into a separate PDB file.

    Parameters
    ----------
    input_pdb : str
        Path to the input PDB file.
    resname : str
        Name of the residue to extract.
    output_pdb : str
        Path to the output PDB file.
    """
    f1 = open(input_pdb, 'r')
    f2 = open(output_pdb, 'w')
    residue = [line for line in f1.readlines() if resname in line]

    for line in residue:
        f2.write(line)
