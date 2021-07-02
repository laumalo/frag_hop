"""
This module containts classes and methods involved in the manipulation of
scaffolds.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from frag_hop.utils.tools import RDKitTools
from frag_hop.replacement.structure import InputStructure, _Structure
import tempfile
import os
import numpy as np
import copy


class Scaffold(_Structure):
    """Class to prepare the selected scaffold for replacement"""

    def __init__(self, path_scaffold=None, path_ligand=None,
                 bonds=None, atom_names=None):

        """
        It initializates a Scaffold object, either directly from a PDB file of
        an scaffold or with the PDB file of a ligand specifing which bonds have
        to be broken (this bonds can be specified by the atoms names of the
        scaffold of the pairs of atoms of each bond).

        Parameters
        ----------
        path_scaffold : str
            Path to a PDB file containing a scaffold.
        path_ligand : str
            Path to a PDB file containin a ligand.
        bonds : list
            List of bonds representing the connectivity points of the scaffold.
        atom_names : list
            List of the atom names comforming the scaffold.
        """
        super().__init__(core=True)

        # Scaffold
        self.path_scaffold = path_scaffold
        self.scaffold = None

        # Ligand
        self.ligand = None
        self.ligand_prepared = None

        if path_scaffold:
            self.initialize_from_pdb(path_scaffold, bonds)

        if path_ligand:
            if bonds is None and atom_names is None:
                raise AttributeError('To select a scaffold from a ligand, ' +
                                     'atom names or bonds have to be selected.')
            else:
                self.initialize_from_ligand(path_ligand, bonds, atom_names)

    def initialize_from_pdb(self, path, bonds):
        """
        It initialices an scaffold as a InputStructure object.

        Parameters
        ----------
        path : str
            Path to a PDB file containing a scaffold.
        bonds : list
            List of bonds representing the connectivity points of the scaffold.
        """
        self.scaffold = InputStructure(path, bonds_link=bonds)

    def initialize_from_ligand(self, path, bonds, atom_names):
        """
        It initialices an scaffold as a InputStructure object from an input
        ligand.

        It extracts the scaffold from the ligand and creates two structures, one
        containing the scaffold and other containing the remaning fragments of
        the ligand.

        Parameters
        ----------
        path : str
            Path to a PDB file containing a ligand.
        bonds : list
            List of bonds representing the connectivity points of the scaffold.
        atom_names : list
            List of the atom names comforming the scaffold.
        """

        def get_scaffold_from_atom_names(molecule, atom_names):
            """
            It gets an scaffold structure from a list of atom names.

            Parameters
            ----------
            molecule : an rdkit.Chem.rdchem.Mol object
                Ligand.
            atom_names : list
                List of the atom names comforming the scaffold.
            """
            raise NotImplementedError

        def get_scaffold_from_bonds(molecule, bonds):
            """
            It gets an scaffold structure from a list of the bonds that connect
            the scaffold with the rest of the molecule.

            Parameters
            ----------
            molecule : an rdkit.Chem.rdchem.Mol object
                Ligand.
            atom_names : list
                List of the atom names comforming the scaffold.
            """
            rdkit_tools = RDKitTools()

            # Checks that there are at least to bonds specified
            if len(bonds) < 2:
                raise ValueError(
                    'To define a scaffold you need at least two bonds.')
            else:
                len_bonds = [len(bond) for bond in bonds]
                if not len_bonds == [2, ] * len(len_bonds):
                    raise ValueError(
                        'Bonds have to be specified by an atom pair.')

            # Breaks the molecule at the selected bonds
            Chem.Kekulize(molecule, clearAromaticFlags=True)
            em = Chem.EditableMol(molecule)
            for bond in bonds:
                idx1 = rdkit_tools.get_atomid_by_atomname(molecule, bond[0])
                idx2 = rdkit_tools.get_atomid_by_atomname(molecule, bond[1])
                em.RemoveBond(idx1, idx2)
            nm = em.GetMol()
            frags = Chem.GetMolFrags(nm, asMols=True, sanitizeFrags=False)

            # Extracts the scaffold out of the ligand
            atom_scaffold = [bond[0].strip() for bond in bonds]
            other_frags = []
            for idx, frag in enumerate(frags):
                atom_names = [atom.GetPDBResidueInfo().GetName().strip()
                              for atom in frag.GetAtoms()]
                is_scaffold = all(elem in atom_names for elem in atom_scaffold)
                if is_scaffold:
                    with tempfile.NamedTemporaryFile(suffix='.pdb') as tmp:
                        Chem.rdmolfiles.MolToPDBFile(frag, tmp.name)
                        self.scaffold = InputStructure(tmp.name,
                                                       bonds_link=bonds)
                else:
                    other_frags.append(frag)

            # Merges the remaining fragments of the ligand
            molecule_frags = other_frags[0]
            for frag in other_frags[1:]:
                molecule_frags = Chem.CombineMols(molecule_frags, frag)

            with tempfile.NamedTemporaryFile(suffix='.pdb') as tmp:
                Chem.rdmolfiles.MolToPDBFile(molecule_frags, tmp.name)
                self.ligand_prepared = InputStructure(tmp.name,bonds_link=bonds)

        self.ligand = InputStructure(path, bonds_link=bonds)

        # Initializates the scaffold by the selected bonds
        if not bonds is None:
            get_scaffold_from_bonds(self.ligand.rdkit_mol, bonds)

        # Initializates the scaffold by the selected atom names
        else:
            if not atom_names is None:
                get_scaffold_from_atom_names(self.ligand.rdkit_mol, atom_names)

    def prepare(self, target):
        """
        It prepares the scaffold structure according to a target molecule in
        other to be replaced for the original scaffold of this target ligand.

        In other to prepare the scaffold, the scaffold bonds are superimpose to
        the connectivity points in the ligand and the hydrogens are removed.

        Parameters
        ----------
        target : a Scaffold object
            Target molecule where the scaffold will be replaced.
        """

        reference_bonds = [el[0] for el in target.ligand.bonds_link]
        scaffold_bonds = [el[0] for el in self.scaffold.bonds_link]

        # Superimposition of the hit scaffold to the target scaffold
        self.scaffold = self.superimpose_fragment_bond(self.scaffold,
                                                       target.ligand,
                                                       scaffold_bonds,
                                                       reference_bonds)

        # Update RDKit molecule with the obtained position
        with tempfile.NamedTemporaryFile(suffix='.pdb') as tmp:
            self.scaffold.structure.save_pdb(tmp.name)
            self.scaffold.rdkit_mol = \
                Chem.rdmolfiles.MolFromPDBFile(tmp.name, removeHs=False)

        ref = InputStructure(self.path_scaffold)
        self.scaffold.rdkit_mol = AllChem.AssignBondOrdersFromTemplate(
            ref.rdkit_mol,
            self.scaffold.rdkit_mol)

        # Rotate to find best structural fit
        self.get_best_rotation_pose(target_scaffold = target.scaffold)

        # Remove hydrogens of the bonds
        self.remove_hydrogens(molecule = self.scaffold)


    def get_best_rotation_pose(self, target_scaffold, rotation_angle=0.15):
        """
        It rotates the new candidate scaffold and computes different metrics
        between each rotation and the original scaffold to obtain the best pose.

        Parameters
        ----------
        target_scaffold : a Scaffold.scaffold object
            Original scaffold of the ligand.
        rotation_angle : float
            Angle of rotation at each iteration. Default: 0.15
        """

        import math

        def distance_plane_fit(molecule1, molecule2):
            """
            It gets the distance betwen the centroid of each molecule. This
            metrics is better used to compare ring fragments where the distance
            between the two geometrical centroids gives an aproximate value for
            the distance of the two planes.

            Parameters
            ----------
            molecule1 : an rdkit.Chem.rdchem.Mol object
                Query molecule.
            molecule2 : an rdkit.Chem.rdchem.Mol object
                Reference molecule.

            Returns
            -------
            distance : float
                Distance between the two centroids.
            """

            from frag_hop.utils.geometry import distance

            pos1 = molecule1.GetConformer().GetPositions()
            pos2 = molecule2.GetConformer().GetPositions()
            centroid1 = np.mean(pos1[:,-3:], axis=0)
            centroid2 = np.mean(pos2[:,-3:], axis=0)
            return distance(centroid1,centroid2)

        def rotate_fragment(radi, molecule, bonds):
            """
            It rotates a molecule along its attachment vector.
            radi : float
                Angle (in radians) to rotate the fragment.
            molecule : an InputStructure object
                Molecule to rotate.
            bonds : list[str]
                List of bonds that define the attachment vector.
            """

            def translation_distance():
                """
                It computes the translation vector between the reference
                fragment position and the rotated fragment.
                """
                ix = rdkit_tools.get_atomid_by_atomname(molecule.rdkit_mol,
                                                        bonds[0])
                p_ref = molecule.rdkit_mol.GetConformer().GetAtomPosition(ix)
                p_ref_rot = np.dot(rot_mat, [p_ref.x, p_ref.y, p_ref.z])
                return p_ref - p_ref_rot

            # Computes vector to be used as axis of rotation
            rdkit_tools = RDKitTools()
            idx1 = rdkit_tools.get_atomid_by_atomname(molecule.rdkit_mol,
                                                      bonds[0])
            idx2 = rdkit_tools.get_atomid_by_atomname(molecule.rdkit_mol,
                                                      bonds[1])

            p1 = molecule.rdkit_mol.GetConformer().GetAtomPosition(idx1)
            p2 = molecule.rdkit_mol.GetConformer().GetAtomPosition(idx2)

            u = np.array([p1.x, p1.y, p1.z]) - np.array([p2.x, p2.y, p2.z])
            u = u / (u**2).sum()**0.5

            # Gets rotation matrix
            from frag_hop.utils.geometry import rotation_matrix_axis
            rot_mat = rotation_matrix_axis(radi, u)

            # Gets Translation vector
            translation = translation_distance()

            # Perform a rotation and translation
            rotated_molecule = copy.deepcopy(molecule)
            for atom in rotated_molecule.rdkit_mol.GetAtoms():
                p = rotated_molecule.rdkit_mol.GetConformer().GetAtomPosition(
                    atom.GetIdx())
                new_p = np.dot(rot_mat, [p.x, p.y, p.z]) + translation
                rotated_molecule.rdkit_mol.GetConformer().SetAtomPosition(
                    atom.GetIdx(), new_p)
            return rotated_molecule

        # Molecules to compare
        molecule1 = self.scaffold
        molecule2 = target_scaffold
        bonds = [el[0] for el in self.scaffold.bonds_link]

        # Computes the number of rotations that will be performed
        num_rotations = int(2 * math.pi / rotation_angle)

        # Rotates the fragment and computes best plane superposition
        d = {}
        for rot in range(num_rotations):
            mol = rotate_fragment(radi=rotation_angle * rot,
                                  molecule = molecule1, bonds = bonds)
            dis = distance_plane_fit(mol.rdkit_mol, molecule2.rdkit_mol)
            d[rot] = dis

        # Get best pose
        best_rot = min(d, key=d.get)
        best_pose = rotate_fragment(radi=rotation_angle * best_rot,
                                    molecule = molecule1, bonds = bonds)
        self.scaffold = best_pose



    def to_file(self, path, file_name='scaffold.pdb'):
        """
        It exports to a PDB structure the scaffold and in case of having obtain
        the scaffold out of a ligand, it also generated the structure of the
        remaining fragments of the ligand.

        Parameters
        ----------
        path : str
            Output path.
        file_name : str
            File name. Default: scafold.pdb
        """
        output_path = os.path.join(path, file_name)

        Chem.rdmolfiles.MolToPDBFile(self.scaffold.rdkit_mol,
                                     output_path)
        if self.ligand_prepared is not None:
            Chem.rdmolfiles.MolToPDBFile(self.ligand_prepared.rdkit_mol,
                                         os.path.join(path, 'lig_prepared.pdb'))


