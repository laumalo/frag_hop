import collections
import sys
import os
import numpy as np

class PDBTools():
    """
    Class that contains all methods applieds to handle PDB files.
    """

    def pbd_uniqname(self, pdb_path):
        """
        It handles the atom names of a PDB file and modifies them to have unique
        atom names.
        Implementation based on: https://github.com/haddocking/pdb-tools
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
        sys.exit(0)

    def extract_chain(self,chain_id, input_file, output_file):
        def select_chain(fhandle, chain_set):
            """Filters the PDB file for specific chain identifiers.
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
        """ Find the rotation matrix that aligns vec1 to vec2
        :param vec1: A 3d "source" vector
        :param vec2: A 3d "destination" vector
        :return mat: A transform matrix (3x3) which when applied to vec1,
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
        ux, uy, uz = u
        return np.array([[np.cos(radi) + ux**2*(1 - np.cos(radi)),
                            ux*uy*(1-np.cos(radi))-uz*np.sin(radi),
                            ux*uz*(1-np.cos(radi)) + uy*np.sin(radi)],
                            [uy*ux*(1-np.cos(radi))+uz*np.sin(radi),
                            np.cos(radi) + uy**2*(1-np.cos(radi)),
                            uy*uz*(1-np.cos(radi))- ux*np.sin(radi)],
                            [uz*ux*(1-np.cos(radi))-uy*np.sin(radi),
                            uz*uy*(1-np.cos(radi))+ ux*np.sin(radi),
                            np.cos(radi)+uz**2*(1-np.cos(radi))]],
                            dtype=np.double)

class RDKitTools():
    def get_atomid_by_atomname(self, mol, atom_name):
        atom_id = [atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetPDBResidueInfo().GetName().strip() == atom_name][0]
        return atom_id

    def get_element_by_atomname(self, mol, atom_name):
        element = [atom.GetSymbol() for atom in mol.GetAtoms()
        if atom.GetPDBResidueInfo().GetName().strip() == atom_name][0]
        return element
