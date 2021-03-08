import collections
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
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    element = line[76:78].strip()
                    if not element:
                        print(
                            'ERROR!! No element found in line {}'.format(
                                line_idx))
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

    def extract_chain(self,chain_id, input_file, output_file):
        """
        It extracts from a PDB file a chain selected by chain_id and exports a
        PDB file with only the structure of these chain.
        Implementation based on: https://github.com/haddocking/pdb-tools
        """
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


def perform_complex_mutation(complex_pdb, new_residue_pdb,
                             resname = 'GRW',
                             output_pdb = 'out/merged_complex.pdb',
                             HIS_atoms = 10):
    """
    Given a PDB with a complex and the PDB of the residue that has to be
    replaced. It replaces from the original complex the residue and exports
    the PDB of the new strucutre.
    """

    f1 = open(complex_pdb,'r')
    f2 = open(new_residue_pdb, 'r')
    f3 = open(output_pdb, 'w')

    # Read original complex
    lines = f1.readlines()

    lines_to_delete, new_lines = ([] for x in range(2))

    # Load new residue atoms
    for line in f2.readlines():
        if 'ATOM' in line:
            new_lines.append(line)

    # Get all residue atoms to delete
    for line_num, line in enumerate(lines):
        if resname in line:
            lines_to_delete.append(line_num)

    lines_to_delete = lines_to_delete[HIS_atoms:]

    # Correct atom numbers of the residue
    for idx, line in enumerate(new_lines):
        new_lines[idx] = line[:6] + str(int(line[6:11]) + \
                                        int(HIS_atoms)).rjust(5) + line[11:]


    new_file_lines = \
        lines[:lines_to_delete[0]] + new_lines + lines[lines_to_delete[-1]+1:]

    for line in new_file_lines:
        f3.write(line)

def extract_residue(pdb_file, resname , out_file = 'out/RES.pdb'):
    """
    It extracts one residue of a PDB file into a separate PDB file.
    """
    f = open(pdb_file,'r')
    f2 = open(out_file, 'w')
    residue = [ line  for line in f.readlines() if resname in line]

    for line in residue[10:]:
        f2.write(line)
