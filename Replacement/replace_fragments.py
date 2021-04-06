import mdtraj as md
from rdkit import Chem
import numpy as np
import os
from utils import RDKitTools

class InputStructure:
    """
    Object to modify scaffold and fragments structures
    """

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


class PrepareFragment:
    """Class to prepare the fragment that will be replaced"""

    def __init__(self, initial_complex, fragment, bond_atoms, out_folder,
                 resname, top_complex=None, top_fragment=None,
                 bond_type='single'):

        self.out = out_folder
        self.resname = resname
        self.initial_complex = InputStructure(initial_complex,
                                              top_file=top_complex)

        self.set_complex_link(bond=bond_atoms[1])

        self.fragment = InputStructure(fragment, top_file=top_fragment)
        self.set_fragment_link(bond=bond_atoms[0])

        self.fragment_prepared = None
        self.prepare_fragment()

    def set_complex_link(self, bond):
        self.initial_complex.set_bond_to_link(bond)

    def set_fragment_link(self, bond):
        self.fragment.set_bond_to_link(bond)

    def __get_atom_indices(self):
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
        h_idx = get_idx(atoms=self.fragment.structure.top.atoms,
                        atom_name=self.fragment.bond_link[1])
        ha_idx = get_idx(atoms=self.fragment.structure.top.atoms,
                         atom_name=self.fragment.bond_link[0])

        # set atom reference indices
        h_ref_idx = get_idx_ref(atoms=self.initial_complex.structure.top.atoms,
                                atom_name=self.initial_complex.bond_link[1])

        ha_ref_idx = get_idx_ref(atoms=self.initial_complex.structure.top.atoms,
                                 atom_name=self.initial_complex.bond_link[0])

        return [h_idx, ha_idx], [h_ref_idx, ha_ref_idx]

    def superimpose_fragment_bond(self):
        """
        Given the attachment vector of the fragment and the attachment vector of
        the scaffold given by the pair of atoms in that bond, it performs the
        rotation and translation needed to the fragment to align the fragment
        attachment vector to the scaffold attachment vector.
        """
        from utils import MathTools

        trajA = self.initial_complex.structure
        trajB = self.fragment.structure

        # Get atoms of the bond indices
        atom_indices, atom_ref_indices = self.__get_atom_indices()

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
        translation_distance = trajA.xyz[0, atom_ref_indices[1], :] - \
            trajB.xyz[0, atom_indices[1], :]
        for i, xyz in enumerate(trajB.xyz[0]):
            trajB.xyz[0][i] = xyz + translation_distance


    def prepare_fragment(self):
        """
        It prepares the input fragment for later fragment replacement techniques
        to generate new molecules.

        It also generates a PDB file with the prepared fragment structure.
        """
        # Prepare fragment
        self.superimpose_fragment_bond()

        # Export PDB structure
        output_path = os.path.join(self.out, 'frag_prepared.pdb')
        self.fragment.structure.save_pdb(output_path)


class Replacer:
    """Class to replace fragments"""

    def __init__(self, ligand_pdb, fragment_pdb, ref_fragment_pdb, bond_atoms,
                 out_folder, bond_type='single'):

        self.out = out_folder

        # Set bonds to break
        self.bond_frag, self.bond_lig = bond_atoms
        self.bond_type = bond_type

        # Ligand
        self.ligand = self.__load_to_rdkit(ligand_pdb)
        self.ligand_prepared = None
        self.original_fragment = None
        self.break_ligand()

        # Fragment
        self.fragment_pdb = fragment_pdb
        self.ref_fragment_pdb = ref_fragment_pdb
        self.fragment = self.__load_to_rdkit(fragment_pdb)
        self.correct_fragment()

        self.rotated_fragment = None
        self.get_best_dihedral_angle()

        # Merge structure
        self.merged = None
        self.get_merged_structure()

    def __load_to_rdkit(self, path):
        """
        It laods a PDB file into an rdkit.Chem.rdchem.Mol object.

        Parameters
        ----------
        path : str
            Path to the PDB file.

        Returns
        -------
        mol : rdkit.Chem.rdchem.Mol object
            Molecule.
        """
        return Chem.rdmolfiles.MolFromPDBFile(path, removeHs=False)

    def break_ligand(self):
        """
        It breaks the ligand at a given bond and generated two structures
        containing the original fragment and the ligand prepared for the
        fragment replacement technique.

        These two structures are also exported as PDB files.
        """
        from rdkit.Chem.rdmolops import FastFindRings
        rdkit_tools = RDKitTools()
        lig = self.ligand
        idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                  self.bond_lig[1])
        idx2 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                  self.bond_lig[0])

        Chem.Kekulize(lig, clearAromaticFlags=True)
        em = Chem.EditableMol(lig)
        em.RemoveBond(idx1, idx2)
        nm = em.GetMol()
        nm.GetAtomWithIdx(idx1).SetNoImplicit(True)
        nm.GetAtomWithIdx(idx2).SetNoImplicit(True)
        frags = Chem.GetMolFrags(nm, asMols=True, sanitizeFrags=False)
        if frags[0].GetNumAtoms() < frags[1].GetNumAtoms():
            self.original_fragment = frags[0]
            self.ligand_prepared = frags[1]
        else:
            self.original_fragment = frags[1]
            FastFindRings(self.original_fragment)
            self.ligand_prepared = frags[0]
            FastFindRings(self.ligand_prepared)

        # Export PDB structures
        ligand_output = os.path.join(self.out, 'lig_prepared.pdb')
        fragment_output = os.path.join(self.out, 'original_fragment.pdb')

        Chem.rdmolfiles.MolToPDBFile(self.ligand_prepared,
                                     ligand_output)
        Chem.rdmolfiles.MolToPDBFile(self.original_fragment,
                                     fragment_output)

    def correct_fragment(self):
        """
        It assigns the bond orders of the prepared fragment from the reference
        fragment (from the library) and removes the hydrogen where the bond will
        take place.
        """

        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Assign bond orders
        ref = self.__load_to_rdkit(self.ref_fragment_pdb)
        frag_bonds = AllChem.AssignBondOrdersFromTemplate(ref, self.fragment)

        # Remove hydrogen
        rdkit_tools = RDKitTools()
        frag_idx = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                      self.bond_frag[1])
        ed_frag = Chem.EditableMol(frag_bonds)
        ed_frag.RemoveAtom(frag_idx)
        Chem.rdmolfiles.MolToPDBFile(ed_frag.GetMol(), self.fragment_pdb)

    def __generate_merged_structure(self):
        """
        It merges two rdkit.Chem.rdchem.Mol objects into a single one.
        """
        self.merged = Chem.CombineMols(self.ligand_prepared, self.fragment)

    def __correct_bond_distance(self):
        """
        It creates a bond between the ligand and the fragment and corrects the
        bond distance of the new bond based on the Dictionary of bonding
        distances in Anstrongs in atom_constants.py
        """
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

        # Create and editable Mol object and correct the bond distance
        ed_merged = Chem.EditableMol(self.merged)
        if self.bond_type == 'single':
            order = Chem.rdchem.BondType.SINGLE
        elif self.bond_type == 'double':
            order = Chem.rdchem.BondType.DOUBLE
        elif self.bond_type == 'triple':
            order = Chem.rdchem.BondType.TRIPLE
        else:
            raise ValueError(
                "Wrong bond type for {self.bond_type}.")

        ed_merged.AddBond(lig_idx, frag_idx, order=order)
        self.merged = ed_merged.GetMol()
        new_distance = BONDING_DISTANCES[lig_element.upper(),
                                         frag_element.upper(), self.bond_type]
        FastFindRings(self.merged)
        rdMolTransforms.SetBondLength(self.merged.GetConformer(), lig_idx,
                                      frag_idx, new_distance)

    def __export_merged_structure(self, resname='GRW', resnum=145,
                                  chain_id='A'):
        """
        Exports the PDB files for the merged structure in the necessary formats
        for later complex modification of a residue and for templetizate the
        residue for PELE simulations.

        Parameters
        ----------
        rename : str
            Residue name
        resnum : int
            Residue number
        chain_id : str
            Chain Id.
        """
        out_file = os.path.join(self.out, 'merged.pdb')
        out_ligand = os.path.join(self.out, resname + '.pdb')

        for heteroatom, path in zip((False, True), (out_file, out_ligand)):
            for a in self.merged.GetAtoms():
                mi = Chem.AtomPDBResidueInfo()
                mi.SetName(a.GetPDBResidueInfo().GetName())
                mi.SetIsHeteroAtom(heteroatom)
                mi.SetResidueName(resname)
                mi.SetResidueNumber(resnum)
                mi.SetChainId(chain_id)
                a.SetMonomerInfo(mi)

            Chem.rdmolfiles.MolToPDBFile(self.merged, path)

            from utils import PDBTools
            PDBModifier = PDBTools()
            PDBModifier.rename_atoms_fragment(path)

    def get_merged_structure(self):
        """
        It generates a merged structure by creating a bond between the prepared
        ligand and fragment. It also exports this structure in a PDB format.
        """
        self.__generate_merged_structure()
        self.__correct_bond_distance()
        self.__export_merged_structure()

    def get_best_dihedral_angle(self, rotation_angle=0.15):
        """
        It rotates the new candidate fragment and computs the USR metrics
        between each rotated fragment and the original fragment to obtain the
        best dihedral angle for the new fragment.

        Parameters
        ----------
        rotation_angle : float
            Angle (in radians) to rotate at each iteration.
        """
        import math

        def rotate_fragment(radi):
            """
            It rotates the object self.rotated_fragment along its attachment
            vector.
            radi : float
                Angle (in radians) to rotate the fragment.
            """

            rdkit_tools = RDKitTools()

            def translation_vector(self):
                """
                It computes the translation vector between the reference
                fragment position and the rotated fragment.
                """
                ix = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                        self.bond_frag[0])
                p_ref = self.fragment.GetConformer().GetAtomPosition(ix)
                p_ref_rot = np.dot(rot_mat, [p_ref.x, p_ref.y, p_ref.z])
                return p_ref - p_ref_rot

            # Resets the rotated fragment
            self.rotated_fragment = self.__load_to_rdkit(self.fragment_pdb)

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
                p = self.rotated_fragment.GetConformer().GetAtomPosition(
                    atom.GetIdx())
                new_p = np.dot(rot_mat, [p.x, p.y, p.z]) + translation
                self.rotated_fragment.GetConformer().SetAtomPosition(
                    atom.GetIdx(), new_p)

        def check_atoms_distances():
            """
            Checks that the position of the fragment in the merged structure
            will not have any of the atoms of the fragment within bond distance
            to the scaffold atoms, leading to the formation of unwanted bonds in
            the new molecule of overlapings between fragment-scaffold.

            Returns
            -------
            keep_position : bool
                True if the position could be keeped for the merged structure.
            """

            from atom_constants import BONDING_DISTANCES
            from utils import distance

            keep_position = True
            for atom in self.rotated_fragment.GetAtoms():
                p1 = self.rotated_fragment.GetConformer().GetAtomPosition(
                    atom.GetIdx())
                for a in self.ligand_prepared.GetAtoms():
                    p2 = self.ligand_prepared.GetConformer().GetAtomPosition(
                        a.GetIdx())
                    #Distance between the atoms
                    distance_atoms = distance(np.array([p1.x, p1.y, p1.z]),
                                              np.array([p2.x, p2.y, p2.z]))
                    # Search for the specific bond distance of the pair of atoms
                    try:
                        cut_off = BONDING_DISTANCES[atom.GetSymbol().upper(),
                                        a.GetSymbol().upper(), 'single']
                    # Default value if not found
                    except KeyError:
                        cut_off = 1.8
                    # Check atoms distances
                    scaffold_atom = a.GetPDBResidueInfo().GetName().strip()
                    frag_atom = atom.GetPDBResidueInfo().GetName().strip()
                    if distance_atoms < cut_off \
                            and not frag_atom == self.bond_frag[0].strip() \
                            and not scaffold_atom == self.bond_lig[1].strip():
                        # Discard this position
                        keep_position = False
                        break
            return keep_position

        def compute_fragment_usr(rdkit_molA, rdkit_molB):
            """
            Performs a shape similarity calculation between two small molecules
            using the Ultrafast Shape Recognition (USR) descriptor with the
            RDKit implementation.

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
                Finds the shape similarity coefficient between two USR
                descriptors.

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

        # Computes the number of rotations that will be performed
        num_rotations = int(2 * math.pi / rotation_angle)

        # Rotates the fragment and computes USR metrics
        d = {}
        for rot in range(num_rotations):
            rotate_fragment(radi=rotation_angle * rot)
            # Checks overlaping fragment-scaffold
            if check_atoms_distances():
                usr_value = compute_fragment_usr(
                    self.original_fragment,
                    self.rotated_fragment)
                d[rot] = usr_value

        # Choose the best position for the fragment
        try:
            best_rot = min(d, key=d.get)
            rotate_fragment(radi=rotation_angle * best_rot)
            self.fragment = self.rotated_fragment
        except ValueError:
            print('Warning: Skipping fragment, there is no possible position ' +
                  'without overlapping. Default initial position returned. ')



