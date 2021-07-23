"""
This module containts classes to handle the input file and basic structures for
fragments and scaffolds.
"""

from rdkit import Chem
import mdtraj as md
import numpy as np


class InputStructure:
    """Object to initialize scaffold and fragment structures"""

    def __init__(self, input_file, bonds_link=None, top_file=None):
        """
        It initialices a InputStructure object.

        Parameters
        ----------
        input_file : str
            Path to the input file.
        bonds_link : list[str]
            List of the two atoms that conforms the bond representing the
            attachment vector.
        top_file : str
            Path to the topology file.
        """

        self.structure = self.__load_to_mdtraj(input_file, top_file)
        self.rdkit_mol = self.__load_to_rdkit(input_file)

        if bonds_link:
            self.set_bonds_to_link(bonds_link)
        else:
            self.bonds_link = bonds_link

    def __load_to_mdtraj(self, input_file, top_file):
        """
        It initializes a trajectorymd.Trajectory object from a PDB file.

        Parameters
        ----------
        input_file : str
            Path to the trajectory file.
        top_file : str
            Path to the topology file.

        Returns
        -------
        structure : a trajectorymd.Trajectory object
            The resulting trajectory, as an md.Trajectory object.
        """
        if top_file:
            return md.load(input_file, top=top_file)
        else:
            return md.load(input_file)

    def __load_to_rdkit(self, input_file):
        """
        It initializes an RDKit's Molecule object from a PDB file.

        Parameters
        ----------
        input_file : str
            Path to the PDB file.

        Returns
        -------
        rdkit_mol : an rdkit.Chem.rdchem.Mol object
            The resulting molecule, as an rdkit.Chem.rdchem.Mol object.
        """
        return Chem.rdmolfiles.MolFromPDBFile(input_file, removeHs=False)

    def set_bonds_to_link(self, new_bonds):
        """
        It sets the bond to link the scaffold and fragment.

        Parameters
        ----------
        new_bond : list[str]
            List of atom names of the atoms involved in the bond.
        """
        def check_new_bonds(new_bonds):
            check = True
            for new_bond in new_bonds:
                if type(new_bond) is list:
                    if not len(new_bond) == 2:
                        check = False
            return check
        if check_new_bonds(new_bonds):
            self.bonds_link = new_bonds
        else:
            raise ValueError(
                "Wrong length for {new_bond}. It must be a list" +
                    " of 2 elements")


class _Structure(object):
    """Class to modify and prepare for replacement scaffold and fragment
    structures"""

    def __init__(self, core=False, terminal=False):
        """
        It initializates a _Structure object

        Parameters
        ----------
        core : bool
            True if the structure is an scaffold.
        terminal : bool
            True if the structure is a terminal fragment.
        """
        self.core = core
        self.terminal = terminal

    def _get_atom_indices(self, mol, mol_ref, bond, bond_ref):
        """
        Gets the atoms indixes of the atoms in the ligand and fragment that will
        be needed for the new connectivity.

        Returns
        -------
        atom_indices : tuple
            Tuple of the pairs of atom indices.
        """
        def get_idx(atoms, atom_name):
            return [atom.index for atom in atoms if atom.name == atom_name][0]

        def get_idx_ref(atoms, atom_name):
            return [atom.index for atom in atoms if atom.name == atom_name][0]

        # Get atom indices
        h_idx = get_idx(atoms=mol.structure.top.atoms,
                        atom_name=bond[1])
        ha_idx = get_idx(atoms=mol.structure.top.atoms,
                         atom_name=bond[0])

        # Get atom reference indices
        h_ref_idx = get_idx_ref(atoms=mol_ref.structure.top.atoms,
                                atom_name=bond_ref[1])

        ha_ref_idx = get_idx_ref(atoms=mol_ref.structure.top.atoms,
                                 atom_name=bond_ref[0])

        return [h_idx, ha_idx], [h_ref_idx, ha_ref_idx]

    def superimpose_fragment_bond(self, molecule, molecule_ref, bond, bond_ref):
        """
        Given the attachment vector of the fragment and the attachment vector
        of the scaffold given by the pair of atoms in that bond, it performs
        the rotation and translation needed to the fragment to align the
        fragment attachment vector to the scaffold attachment vector.

        Parameters
        ----------
        molecule : an InputStructure object
            Molecule that will be superimpose.
        molecule_ref : an InputStructure object
            Reference molecule.
        bond : list[str]
            Selected bond in the molecule.
        bond_ref : list[str]
            Selected bond in the reference molecule.
        """
        from frag_hop.utils.geometry import rotation_matrix_from_vectors

        trajA = molecule_ref.structure
        trajB = molecule.structure

        bond, bond_ref = \
            (bond[0], bond_ref[0]) if self.terminal else (bond, bond_ref)

        # Get atoms of the bond indices
        atom_indices, atom_ref_indices = self._get_atom_indices(
            mol=molecule,
            mol_ref=molecule_ref,
            bond=bond,
            bond_ref=bond_ref)

        atom_indices = atom_indices[::-1] if self.terminal else atom_indices

        # Compute reference vectors
        vec_ref = trajA.xyz[0, atom_ref_indices[0], :] - \
            trajA.xyz[0, atom_ref_indices[1], :]
        vec_ref = vec_ref / (vec_ref**2).sum()**0.5
        vec = trajB.xyz[0, atom_indices[0], :] - \
            trajB.xyz[0, atom_indices[1], :]
        vec = vec / (vec**2).sum()**0.5

        # Apply rotation
        rot = rotation_matrix_from_vectors(vec, vec_ref)
        for i, xyz in enumerate(trajB.xyz[0]):
            trajB.xyz[0][i] = np.dot(rot, xyz)

        # Apply translation between the heavy atoms
        translation_distance = trajA.xyz[0, atom_ref_indices[1], :] - \
            trajB.xyz[0, atom_indices[1], :]
        for i, xyz in enumerate(trajB.xyz[0]):
            trajB.xyz[0][i] = xyz + translation_distance
        return molecule

    def remove_hydrogens(self, molecule):
        """
        It removes all the hydrogen atoms in the bonds to link list of a
        molecule.

        Parameters
        ----------
        molecule : an InputStructure object
            Molecule to delete hydrogen atoms from.
        """
        from rdkit.Chem.rdmolops import FastFindRings

        # Remove hydrogens
        atoms_list = molecule.bonds_link
        Hs = [el for el in [j for i in atoms_list for j in i] if 'H' in el]

        for H in Hs:
            idx = [atom.GetIdx() for atom in molecule.rdkit_mol.GetAtoms()
                   if atom.GetPDBResidueInfo().GetName().strip() == H][0]

            ed_scaffold = Chem.EditableMol(molecule.rdkit_mol)
            ed_scaffold.RemoveAtom(idx)
            molecule.rdkit_mol = ed_scaffold.GetMol()

        # Update the RDKit representation of the molecule
        FastFindRings(molecule.rdkit_mol)
        molecule.rdkit_mol.UpdatePropertyCache(strict=False)
