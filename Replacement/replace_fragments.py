import mdtraj as md
from rdkit import Chem
import numpy as np
from utils import RDKitTools


class InputStructure:
    """Object to modify scaffold and fragments structures"""

    def __init__(self, input_file, bond_link=None, top_file=None):
        self.input_file = input_file
        self.top_file = top_file
        self.structure = None
        self.__load_to_mdtraj()

        if bond_link:
            self.set_bond_to_link(bond_link)
        else:
            self.bond_link = bond_link

    def __load_to_mdtraj(self):
        if self.top_file:
            self.structure = md.load(self.input_file, top_file=self.top_file)
        else:
            self.structure = md.load(self.input_file)

    def set_bond_to_link(self, new_bond):
        if type(new_bond) is list:
            if len(new_bond) == 2:
                self.bond_link = new_bond
            else:
                raise ValueError(
                    "Wrong length for {new_bond}. It must be a list of 2 elements")
        else:
            raise TypeError("Wrong type for {new_bond}. It must be a 'list'")


class PrepareSystem:
    """Class to replace fragments"""
    def __init__(self, initial_complex, fragment, bond_atoms,
                 resname = 'GRW', resnum = 145, res_chain = 'A',
                 top_complex=None, top_fragment=None, bond_type='single'):

        self.initial_complex = InputStructure(initial_complex,
                                              top_file=top_complex)
        # Residue information
        self.resname = resname
        self.resnum = resnum
        self.res_chain = res_chain

        self.fragment = InputStructure(fragment, top_file=top_fragment)

        self.set_complex_link(bond = bond_atoms[1])
        self.set_fragment_link(bond = bond_atoms[0])

        self.fragment_prepared = None
        self.prepare_fragment()

        self.ligand = None
        self.bond_type = bond_type

    def set_complex_link(self, bond):
        self.initial_complex.set_bond_to_link(bond)

    def set_fragment_link(self, bond):
        self.fragment.set_bond_to_link(bond)

    def get_atom_indices(self):
        # set atom indices
        h_idx = [atom.index for atom in self.fragment.structure.top.atoms
                 if atom.name == self.fragment.bond_link[1]][0]
        ha_idx = [atom.index for atom in self.fragment.structure.top.atoms
                  if atom.name == self.fragment.bond_link[0]][0]

        # set atom reference indices
        h_ref_idx = [atom.index for atom in self.initial_complex.structure.top.atoms
                     if atom.name == self.initial_complex.bond_link[1] and atom.residue.name == self.resname][0]
        ha_ref_idx = [atom.index for atom in self.initial_complex.structure.top.atoms
                      if atom.name == self.initial_complex.bond_link[0] and atom.residue.name == self.resname][0]

        return [h_idx, ha_idx], [h_ref_idx, ha_ref_idx]

    def superimpose_fragment_bond(self, trajA=None,
                                  trajB=None):
        """
        """
        from utils import MathTools

        trajA = self.initial_complex.structure
        trajB = self.fragment.structure

        # Get atoms of the bond indices
        atom_indices, atom_ref_indices = self.get_atom_indices()

        # Compute reference vectors
        vec_ref = trajA.xyz[0, atom_ref_indices[0], :] - \
            trajA.xyz[0, atom_ref_indices[1], :]
        vec_ref = vec_ref / (vec_ref**2).sum()**0.5
        vec = trajB.xyz[0, atom_indices[0], :] - \
            trajB.xyz[0, atom_indices[1], :]
        vec = vec / (vec**2).sum()**0.5

        # Apply rotation
        MTools = MathTools()
        rot = MTools.rotation_matrix_from_vectors(vec, vec_ref)
        for i, xyz in enumerate(trajB.xyz[0]):
            trajB.xyz[0][i] = np.dot(rot, xyz)

        # Apply translation between the heavy atoms
        distance = trajA.xyz[0, atom_ref_indices[1], :] - \
            trajB.xyz[0, atom_indices[1], :]
        for i, xyz in enumerate(trajB.xyz[0]):
            trajB.xyz[0][i] = xyz + distance

    def remove_H(self):
        """
        Remove the H in the fragment that will be the bond with the molecule
        """
        atom_indices = [
            atom.index for atom in self.fragment.structure.top.atoms]
        h_idx = [atom.index for atom in self.fragment.structure.top.atoms
                 if atom.name == self.fragment.bond_link[1]]
        new_indices = list(set(atom_indices) - set(h_idx))
        self.fragment.structure = self.fragment.structure.atom_slice(
            new_indices)

    def prepare_fragment(self):

        # Prepare fragment
        self.superimpose_fragment_bond()
        self.remove_H()
        self.fragment.structure.save_pdb('out/frag_prepared.pdb')

class Replacer:
    def __init__(self, ligand_pdb, fragment_pdb, bond_atoms,
                 bond_type='single'):

        # Set bonds to break
        self.bond_frag = bond_atoms[0]
        self.bond_lig = bond_atoms[1]
        self.bond_type = bond_type

        # Ligand
        self.ligand = self.__load_to_rdkit(ligand_pdb)
        self.ligand_prepared = None
        self.original_fragment = None
        self.__break_ligand()

        # Fragment
        self.fragment = self.__load_to_rdkit(fragment_pdb)
        self.rotated_fragment = self.fragment
        self.__get_best_dihedral_angle()

        # Merge structure
        self.merged = None
        self.__generate_merged_structure()
        self.__correct_bond_distance()
        self.__export_merged_structure()


    def __break_ligand(self):
        from rdkit.Chem.rdmolops import FastFindRings
        rdkit_tools = RDKitTools()
        lig = self.ligand
        idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                  self.bond_lig[1])
        idx2 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                  self.bond_lig[0])
        Chem.Kekulize(lig,clearAromaticFlags=True)
        em = Chem.EditableMol(lig)
        em.RemoveBond(idx1,idx2)
        nm = em.GetMol()
        nm.GetAtomWithIdx(idx1).SetNoImplicit(True)
        nm.GetAtomWithIdx(idx2).SetNoImplicit(True)
        frags = Chem.GetMolFrags(nm,asMols=True, sanitizeFrags = False)
        if frags[0].GetNumAtoms() < frags[1].GetNumAtoms():
            self.original_fragment = frags[0]
            self.ligand_prepared = frags[1]
        else:
            self.original_fragment = frags[1]
            FastFindRings(self.original_fragment)
            self.ligand_prepared = frags[0]
            FastFindRings(self.ligand_prepared)

        Chem.rdmolfiles.MolToPDBFile(self.ligand_prepared,
                                     'out/lig_prepared.pdb')
        Chem.rdmolfiles.MolToPDBFile(self.original_fragment,
                                     'out/original_fragment.pdb')

    def __load_to_rdkit(self, path):
        return Chem.rdmolfiles.MolFromPDBFile(path, removeHs=False)

    def __generate_merged_structure(self):
        self.merged = Chem.CombineMols(self.ligand_prepared, self.fragment)

    def __correct_bond_distance(self):
        from rdkit.Chem import rdMolTransforms
        from rdkit.Chem.rdmolops import FastFindRings
        from atom_constants import BONDING_DISTANCES

        rdkit_tools = RDKitTools()

        # Get the ligand ane fragment atom id and atom element for the bond
        lig_idx = rdkit_tools.get_atomid_by_atomname(self.ligand_prepared,
                                                     self.bond_lig[1])
        frag_idx = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                     self.bond_frag[0]) \
                + len(self.ligand_prepared.GetAtoms())

        lig_element = rdkit_tools.get_element_by_atomname(self.ligand_prepared,
                                                          self.bond_lig[1])

        frag_element = rdkit_tools.get_element_by_atomname(self.fragment,
                                                           self.bond_frag[0])

        # Creata and editable Mol object and correct the bond distance
        ed_merged = Chem.EditableMol(self.merged)
        if self.bond_type == 'single':
            order = Chem.rdchem.BondType.SINGLE
        ed_merged.AddBond(lig_idx, frag_idx, order=order)
        self.merged = ed_merged.GetMol()
        new_distance = BONDING_DISTANCES[lig_element.upper(),
                                         frag_element.upper(), self.bond_type]
        FastFindRings(self.merged)
        rdMolTransforms.SetBondLength(self.merged.GetConformer(), lig_idx,
                                      frag_idx, new_distance)

    def __export_merged_structure(self, resname = 'GRW', resnum = 145,
                                  chain_id = 'A', out_file='out/merged.pdb'):

        for a in self.merged.GetAtoms():
            mi = Chem.AtomPDBResidueInfo()
            mi.SetName(a.GetPDBResidueInfo().GetName())
            mi.SetIsHeteroAtom(False)
            mi.SetResidueName(resname)
            mi.SetResidueNumber(resnum)
            mi.SetChainId(chain_id)
            a.SetMonomerInfo(mi)

        Chem.rdmolfiles.MolToPDBFile(self.merged, out_file)

        from utils import PDBTools
        PDBModifier = PDBTools()
        PDBModifier.pbd_uniqname(out_file)

    def rotate_fragment(self, radi):
        """
        It rotates the fragment for the bonding vector
        """

        rdkit_tools = RDKitTools()

        def translation_vector(self):

            ix = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                     self.bond_frag[0])
            p_ref = self.rotated_fragment.GetConformer().GetAtomPosition(ix)
            p_ref_rot = np.dot(rot_mat, [p_ref.x, p_ref.y, p_ref.z])
            return p_ref - p_ref_rot

        # Computes vector to be used as axis of rotation
        idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                  self.bond_lig[1])
        idx2 = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                  self.bond_frag[0])

        p1 = self.ligand.GetConformer().GetAtomPosition(idx1)
        p2 = self.fragment.GetConformer().GetAtomPosition(idx2)

        u = np.array([p1.x, p1.y, p1.z]) - np.array([p2.x, p2.y, p2.z])
        u = u / (u**2).sum()**0.5

        from utils import MathTools
        M = MathTools()

        # Gets rotation matrix
        rot_mat = M.rotation_matrix_axis(radi, u)

        # Gets Translation vector
        translation = translation_vector(self)

        # Perform a rotation and translation
        for atom in self.rotated_fragment.GetAtoms():
            p = \
            self.rotated_fragment.GetConformer().GetAtomPosition(atom.GetIdx())
            new_p = np.dot(rot_mat, [p.x, p.y, p.z]) + translation
            self.rotated_fragment.GetConformer().SetAtomPosition(atom.GetIdx(),
                                                                 new_p)

    def compute_fragment_rmsd(self, trajA, trajB):
        """
        It computes the RMSD between two trajectories.

        # IT HAS TO BE MODIFIED FOR FRAGMENTS THAT ARE NOT REALLY SIMILAR
        """
        def traj_to_coords_list(traj):
            """
            Converts a <mdtraj.core.trajectory.Trajectory> object into a list of
            3-tuples.
            """
            coords = traj.xyz.tolist()[0]
            x = tuple([xyz[0] for xyz in coords])
            y = tuple([xyz[1] for xyz in coords])
            z = tuple([xyz[2] for xyz in coords])
            return [x, y, z]

        def squared_distance(coordsA, coordsB):
            """
            Find the squared distance between two 3-tuples
            """
            sqrdist = sum((a - b)**2 for a, b in zip(coordsA, coordsB))
            return sqrdist

        def rmsd(allcoordsA, allcoordsB):
            """
            Finds the RMSD between two lists of 3-tuples.

            Parameters
            ----------

            allcoordsA : list
                Coordenates of the first structure.
            allcoordsN : list
                Cornenates of the second structure.

            Returns
            -------

            rmsd : float
                RMSD.
            """
            import math

            deviation = sum(squared_distance(atomA, atomB) for
                            (atomA, atomB) in zip(allcoordsA, allcoordsB))
            return math.sqrt(deviation / float(len(allcoordsA)))

        # RMSD
        allcoordsA = traj_to_coords_list(trajA)
        allcoordsB = traj_to_coords_list(trajB)
        return rmsd(allcoordsA, allcoordsB)

    def compute_fragment_usr(self, rdkit_molA, rdkit_molB):
        """
        Performs a shape similarity calculation between two small molecules
        using the Ultrafast Shape Recognition (USR) descriptor with the RDKit
        implementation.

        Parameters
        ----------
        rdkit_molA : an rdkit.Chem.rdchem.Mol object
            RDKit's Molecule object of the reference fragment.
        rdkit_molB : an rdkit.Chem.rdchem.Mol object
            RDKit's Molecule object of the candidate fragment orientation.

        Returns
        -------
        s : float
            Shape similarity coefficient.
        """
        def s(USR1, USR2):
            """
            Finds the shape similarity coefficient between two USR descriptors.

            Implementation based on: Front. Chem., 25 July 2018 |
            https://doi.org/10.3389/fchem.2018.00315

            Parameters
            ----------

            USR1 : list
                USR descriptor for one conformer of a molecule.
            USR2 : list
                USR descriptor for one conformer of a molecule.

            Returns
            -------

            rmsd : float
                Shape Similarity coefficient.
            """

            deviation = sum((a - b)**2 for a, b in zip(USR1, USR2))
            return 1 / (1 + (1 / 12) * deviation)

        from rdkit.Chem import rdMolDescriptors

        # Shape similarity coefficient
        usr1 = rdMolDescriptors.GetUSR(rdkit_molA)
        usr2 = rdMolDescriptors.GetUSR(rdkit_molB)
        return s(usr1, usr2)

    def __get_best_dihedral_angle(self, angle_rotation=0.15):

        import math
        # Computes the number of rotations that will be performed
        num_rotations = int(2 * math.pi / angle_rotation)

        # Rotates the fragment and computs USR metrics
        d = {}
        for rot in range(num_rotations):
            self.rotate_fragment(radi=angle_rotation * rot)
            usr_value = self.compute_fragment_usr(
                self.original_fragment,
                self.rotated_fragment)
            d[rot] = usr_value

        # Choose the best position for the fragment
        best_rot = min(d, key=d.get)
        self.rotate_fragment(radi=best_rot)
        self.fragment = self.rotated_fragment
