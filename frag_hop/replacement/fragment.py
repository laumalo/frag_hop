"""
This module containts classes and methods involved in the manipulation of
fragments.
"""

# General imports
import mdtraj as md
import numpy as np
import os
import copy


class InputStructure:
    """Object to modify scaffold and fragments structures"""

    def __init__(self, input_file, bond_link=None, top_file=None):
        """
        It initialices a InputStructure object.

        Parameters
        ----------
        input_file : str
            Path to the input file.
        bond_link : list[str]
            List of the two atoms that conforms the bond representing the
            attachment vector.
        top_file : str
            Path to the topology file.
        """

        self.structure = self.__load_to_mdtraj(input_file, top_file)

        if bond_link:
            self.set_bond_to_link(bond_link)
        else:
            self.bond_link = bond_link

    def __load_to_mdtraj(self, input_file, top_file):
        """
        It loads a PDB file into a mdtraj.structure object.

        Parameters
        ----------
        input_file : str
            Path to the trajectory file.
        top_file : str
            Path to the topology file.
        """
        if top_file:
            return md.load(input_file, top=top_file)
        else:
            return md.load(input_file)

    def set_bond_to_link(self, new_bond):
        """
        It sets the bond to link the scaffold and fragment.

        Parameters
        ----------
        new_bond : list[str]
            List of atom names of the atoms involved in the bond.
        """
        if type(new_bond) is list:
            if len(new_bond) == 2:
                self.bond_link = new_bond
            else:
                raise ValueError(
                    "Wrong length for {new_bond}. It must be a list" +
                    " of 2 elements")
        else:
            raise TypeError("Wrong type for {new_bond}. It must be a 'list'")


class Fragment:
    """Class to prepare the selected fragment for replacement"""

    def __init__(self, initial_complex, fragment, bond_atoms,
                 resname, top_complex=None, top_fragment=None):
        """
        It initialices a Fragment object and generates the prepared strcuture
        for the fragment according to the target complex.

        Parameters
        ----------
        initial_complex : str
            Path to the protein-ligand complex in which a fragment will be
            replaced.
        fragment : str
            Path to the hit fragment.
        bond_atoms : list[list[str]]
            List of atoms that represent the attachment vectors of the fragment
            and scaffold.
        resname : str
            Residue name from the complex where the replacement will be
            performed.
        """

        # Complex structure
        self.resname = resname
        self.initial_complex = InputStructure(initial_complex,
                                              top_file=top_complex)
        self.set_complex_link(bond=bond_atoms[1])

        # Initial fragment
        self.initial_fragment = InputStructure(fragment, top_file=top_fragment)
        self.set_fragment_link(bond=bond_atoms[0])

        # Prepared fragment
        self.prepared_fragment = copy.deepcopy(self.initial_fragment)
        self.prepare_fragment()

    def set_complex_link(self, bond):
        self.initial_complex.set_bond_to_link(bond)

    def set_fragment_link(self, bond):
        self.initial_fragment.set_bond_to_link(bond)

    def _get_atom_indices(self):
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

        def get_idx_ref(atoms, atom_name, resname=self.resname):
            return [atom.index for atom in atoms if atom.name == atom_name
                    and atom.residue.name == resname][0]

        # set atom indices
        h_idx = get_idx(atoms=self.initial_fragment.structure.top.atoms,
                        atom_name=self.initial_fragment.bond_link[1])
        ha_idx = get_idx(atoms=self.initial_fragment.structure.top.atoms,
                         atom_name=self.initial_fragment.bond_link[0])

        # set atom reference indices
        h_ref_idx = get_idx_ref(atoms=self.initial_complex.structure.top.atoms,
                                atom_name=self.initial_complex.bond_link[1])

        ha_ref_idx = get_idx_ref(atoms=self.initial_complex.structure.top.atoms,
                                 atom_name=self.initial_complex.bond_link[0])

        return [h_idx, ha_idx], [h_ref_idx, ha_ref_idx]

    def _superimpose_fragment_bond(self):
        """
        Given the attachment vector of the fragment and the attachment vector
        of the scaffold given by the pair of atoms in that bond, it performs
        the rotation and translation needed to the fragment to align the
        fragment attachment vector to the scaffold attachment vector.
        """
        from frag_hop.utils.geometry import rotation_matrix_from_vectors

        trajA = self.initial_complex.structure
        trajB = self.prepared_fragment.structure

        # Get atoms of the bond indices
        atom_indices, atom_ref_indices = self._get_atom_indices()

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

    def prepare_fragment(self):
        """
        It prepares the input fragment for later fragment replacement techniques
        to generate new molecules.
        """
        # Prepare fragment
        self._superimpose_fragment_bond()

    def to_file(self, output_path, file_name='frag_prepared.pdb'):
        """
        Exports the PDB structure of the prepared fragment.

        Parameters
        ----------
        output_path : str
            Output path for the PDB structure.
        file_name : str
            Output PDB file name. Default: frag_prepared.pdb
        """
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        self.prepared_fragment.structure.save_pdb(os.path.join(output_path,
                                                               file_name))
