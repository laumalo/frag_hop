"""
This module containts classes that handle the fragment replacement in a ligand.
"""

#General imports
from rdkit import Chem
import numpy as np
import os

#Local imports
from frag_hop.utils.geometry import rotation_matrix_axis
from frag_hop.utils.tools import RDKitTools


class FragmentReplacer(object):
    """Class to replace fragments"""

    def __init__(self, hit_fragment, target):
        """
        It initialices a Replacer object.

        Parameters
        ----------
        ligand_pdb : str
            Path to the ligand involved in the fragment replacement.
        fragment_pdb : str
            Path to the prepared fragment to be repalced for.
        ref_fragment_pdb : str
            Path to the fragment in the library.
        bond_atoms : list[list[str]]
            List of atoms that represent the attachment vectors of the fragment
            and scaffold.
        resname : str
            Residue name of the ligand.
        resnum : str
            Residu sequence number of the ligand
        out_folder : str
            Path to the output folder. Default: out
        bond_type : str
            Type of bond of the connection scaffold-fragment. Dafault: single
        """

        # Fragment
        self.fragment = hit_fragment.rdkit_mol
        self.bond_fragment = hit_fragment.bonds_link[0]

        # Ligand
        self.ligand = target.ligand_prepared.rdkit_mol
        self.bond_ligand = target.ligand_prepared.bonds_link[0]

        # Perform replacement
        self.original_fragment = target.fragment.rdkit_mol
        self.best_fragment_position = self.get_best_dihedral_angle(
                                        ref_fragment = self.original_fragment,
                                        fragment = self.fragment)

        self.new_molecule = None
        self.replacement()

    def replacement(self):
        self.merge_structures()


    def merge_structures(self, bond_type = 'single'):
        """
        It creates a bond between the ligand and the fragment and corrects the
        bond distance of the new bond based on the Dictionary of bonding
        distances in Anstrongs in atom_constants.py
        """
        from rdkit.Chem import rdMolTransforms
        from rdkit.Chem.rdmolops import FastFindRings
        from frag_hop.data.parameters.atom_constants import BONDING_DISTANCES

        merged_structure= Chem.CombineMols(self.ligand,
                                           self.best_fragment_position)

        # Get the ligand and fragment atom id and atom element for the bond
        rdkit_tools = RDKitTools()
        atom1 = self.bond_ligand[0]
        atom2 = self.bond_fragment[0]
        idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand, atom1)
        idx2 = rdkit_tools.get_atomid_by_atomname(self.fragment, atom2) \
            + len(self.ligand.GetAtoms())

        el1 = rdkit_tools.get_element_by_atomname(self.ligand, atom1).upper()
        el2 = rdkit_tools.get_element_by_atomname(self.fragment, atom2).upper()

        # Create and editable Mol object and correct the bond distance
        ed_merged = Chem.EditableMol(merged_structure)
        bond_order = rdkit_tools.get_bond_type(bond_type)
        ed_merged.AddBond(idx1, idx2, order=bond_order)
        merged_structure = ed_merged.GetMol()
        new_distance = BONDING_DISTANCES[el1,el2, bond_type]
        FastFindRings(merged_structure)
        rdMolTransforms.SetBondLength(
                merged_structure.GetConformer(), idx1, idx2, new_distance)

        self.new_molecule = merged_structure

    def to_file(self, out_folder, resname = 'LIG', resnum =  '1', chain_id='A'):
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
        from rdkit.Chem.rdmolops import FastFindRings

        out_file = os.path.join(out_folder, 'new_molecule.pdb')
        out_ligand = os.path.join(out_folder, resname + '.pdb')

        rdkit_tools = RDKitTools()
        for heteroatom, path in zip((True, False), (out_ligand,out_file)):
            if heteroatom is False:
                try:
                    em = Chem.EditableMol(self.new_molecule)
                    hn_idx = rdkit_tools.get_atomid_by_atomname(self.new_molecule,
                                                                'HN')
                    em.RemoveAtom(hn_idx)
                    self.new_molecule = em.GetMol()
                    em = Chem.EditableMol(self.new_molecule)
                    hxt_idx = rdkit_tools.get_atomid_by_atomname(self.new_molecule,
                                                                 'HXT')
                    em.RemoveAtom(hxt_idx)
                    self.new_molecule = em.GetMol()
                except Exception:
                    pass
            for a in self.new_molecule.GetAtoms():
                mi = Chem.AtomPDBResidueInfo()
                mi.SetName(a.GetPDBResidueInfo().GetName())
                mi.SetIsHeteroAtom(heteroatom)
                mi.SetResidueName(resname)
                mi.SetResidueNumber(int(resnum))
                mi.SetChainId(chain_id)
                a.SetMonomerInfo(mi)

            FastFindRings(self.new_molecule)
            self.new_molecule.UpdatePropertyCache(strict=False)
            Chem.rdmolfiles.MolToPDBFile(self.new_molecule, path)

            from utils import PDBTools
            PDBModifier = PDBTools()
            PDBModifier.rename_atoms_fragment(path)


    def get_best_dihedral_angle(self, ref_fragment, fragment,
                                rotation_angle=0.15):
        """
        It rotates the new candidate fragment and computs the USR metrics
        between each rotated fragment and the original fragment to obtain the
        best dihedral angle for the new fragment.

        Parameters
        ----------
        rotation_angle : float
            Angle (in radians) to rotate at each iteration. Dafault: 0.15 rad
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
                                                        self.bond_fragment[0])
                p_ref = self.fragment.GetConformer().GetAtomPosition(ix)
                p_ref_rot = np.dot(rot_mat, [p_ref.x, p_ref.y, p_ref.z])
                return p_ref - p_ref_rot

            # Resets the rotated fragment
            rotated_fragment = self.fragment

            # Computes vector to be used as axis of rotation
            new_bond = [self.bond_ligand[0], self.bond_fragment[0]]

            idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                      new_bond[0])
            idx2 = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                      new_bond[1])

            p1 = self.ligand.GetConformer().GetAtomPosition(idx1)
            p2 = self.fragment.GetConformer().GetAtomPosition(idx2)

            u = np.array([p1.x, p1.y, p1.z]) - np.array([p2.x, p2.y, p2.z])
            u = u / (u**2).sum()**0.5

            # Gets rotation matrix
            rot_mat = rotation_matrix_axis(radi, u)

            # Gets Translation vector
            translation = translation_vector(self)

            # Perform a rotation and translation
            for atom in rotated_fragment.GetAtoms():
                p = rotated_fragment.GetConformer().GetAtomPosition(
                    atom.GetIdx())
                new_p = np.dot(rot_mat, [p.x, p.y, p.z]) + translation
                rotated_fragment.GetConformer().SetAtomPosition(atom.GetIdx(),
                                                                new_p)
            return rotated_fragment

        def check_atoms_distances(rotated_fragment = None):
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

            from frag_hop.data.parameters.atom_constants import BONDING_DISTANCES
            from frag_hop.utils.geometry import distance

            keep_position = True

            molecule1 = rotated_fragment
            atom_bonds1 = self.bond_fragment
            atom_bonds2 =self.bond_ligand
            molecule2 = self.ligand

            for atom in molecule1.GetAtoms():
                p1 = molecule1.GetConformer().GetAtomPosition(atom.GetIdx())
                for a in molecule2.GetAtoms():
                    p2 = molecule2.GetConformer().GetAtomPosition(a.GetIdx())
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
                    mol2_atom = a.GetPDBResidueInfo().GetName().strip()
                    mol1_atom = atom.GetPDBResidueInfo().GetName().strip()
                    if distance_atoms < cut_off \
                            and not mol1_atom in atom_bonds1 \
                            and not mol2_atom in atom_bonds2:
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
            rotated_fragment = rotate_fragment(radi=rotation_angle * rot)
            # Checks overlaping fragment-scaffold
            if check_atoms_distances(rotated_fragment):
                usr_value = compute_fragment_usr(
                    self.original_fragment,
                    rotated_fragment)
                d[rot] = usr_value

        # Choose the best position for the fragment
        try:
            best_rot = min(d, key=d.get)
            best_fragment = rotate_fragment(radi=rotation_angle * best_rot)

        except ValueError:
            raise Exception('Warning: Skipping fragment, there is no  ' +
                            'possible position without overlapping.')
        return best_fragment

class ScaffoldReplacer(object):
    """Class to replace scaffolds"""

    def __init__(self, hit_scaffold, target_ligand):
        """
        It initialices a ScaffoldReplacer object.

        Parameters
        ----------
        hit_scaffold : an InputStructure object
            New hit scaffold.
        target_ligand : an InputStructure object
            Prepared ligand where the scaffold will be replaced.
        """

        # Scaffold
        self.scaffold = hit_scaffold.rdkit_mol
        self.bonds_scaffold = hit_scaffold.bonds_link

        # Ligand
        self.ligand = target_ligand.rdkit_mol
        self.bonds_ligand = target_ligand.bonds_link

        # Perform replacement
        self.new_molecule = None
        self.replacement()

    def replacement(self):
        """
        It performs the scaffold replacement to generate a new molecule.
        """

        # Merge structures
        self.merge_structures()

        # Check for internal clashes
        valid_molecule = self.check_atom_distances()
        if not valid_molecule:
            raise Exception(
                'The selected scaffold produces internal clashes when replaced '
                + 'in the selected ligand structure. A valid molecule could '
                + 'not be generated.')


    def merge_structures(self, bond_type = 'single'):
        """
        It merges the scaffold structure and the ligand structure in order to
        create the new molecule.
        """
        from frag_hop.data.parameters.atom_constants import BONDING_DISTANCES
        from rdkit.Chem import rdMolTransforms
        from rdkit.Chem.rdmolops import FastFindRings

        merged_structure = Chem.CombineMols(self.ligand, self.scaffold)

        # Create new bonds
        rdkit_tools = RDKitTools()
        for atom1, atom2 in zip([el[1] for el in self.bonds_ligand],
                                [el[0] for el in self.bonds_scaffold]):
            # Indexes
            idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand, atom1)
            idx2 = rdkit_tools.get_atomid_by_atomname(self.scaffold, atom2) \
                + len(self.ligand.GetAtoms())

            # Elements
            el1 = rdkit_tools.get_element_by_atomname(self.ligand, atom1)
            el2 = rdkit_tools.get_element_by_atomname(self.scaffold, atom2)

            # Add bond
            ed_merged = Chem.EditableMol(merged_structure)
            bond_order = rdkit_tools.get_bond_type(bond_type)
            ed_merged.AddBond(idx1, idx2, order=bond_order)
            merged_structure = ed_merged.GetMol()

            # Prepare molecule and correct bond distance
            Chem.SanitizeMol(merged_structure)
            new_distance = BONDING_DISTANCES[el1, el2, bond_type]
            FastFindRings(merged_structure)
            rdMolTransforms.SetBondLength(
                merged_structure.GetConformer(), idx1, idx2, new_distance)

        self.new_molecule = merged_structure

    def check_atom_distances(self):
        """
        Checks that the position of the atoms in the merged structure
        will not have any of the atoms of the fragment within bond distance
        to the scaffold atoms, leading to the formation of unwanted bonds in
        the new molecule of overlapings between fragment-scaffold.

        Returns
        -------
        keep_position : bool
            True if the position could be keeped for the merged structure.
        """

        from frag_hop.data.parameters.atom_constants import BONDING_DISTANCES
        from frag_hop.utils.geometry import distance

        keep_position = True

        molecule1 = self.scaffold
        atom_bonds1 = sum(self.bonds_scaffold, [])
        atom_bonds2 = sum(self.bonds_ligand, [])
        molecule2 = self.ligand

        for atom in molecule1.GetAtoms():
            p1 = molecule1.GetConformer().GetAtomPosition(atom.GetIdx())
            for a in molecule2.GetAtoms():
                p2 = molecule2.GetConformer().GetAtomPosition(a.GetIdx())
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
                mol2_atom = a.GetPDBResidueInfo().GetName().strip()
                mol1_atom = atom.GetPDBResidueInfo().GetName().strip()
                if distance_atoms < cut_off \
                        and not mol1_atom in atom_bonds1 \
                        and not mol2_atom in atom_bonds2:
                    # Discard this position
                    keep_position = False
                    break
        return keep_position

    def to_file(self, path, file_name='new_molecule.pdb', resname = 'LIG',
                resnum = 1, chain_id = 'L', heteroatom = True):
        """
        Exports PDB
        """
        output_path = os.path.join(path, file_name)
        from rdkit.Chem.rdmolops import FastFindRings

        for a in self.new_molecule.GetAtoms():
            mi = Chem.AtomPDBResidueInfo()
            mi.SetName(a.GetPDBResidueInfo().GetName())
            mi.SetIsHeteroAtom(heteroatom)
            mi.SetResidueName(resname)
            mi.SetResidueNumber(int(resnum))
            mi.SetChainId(chain_id)
            a.SetMonomerInfo(mi)

        FastFindRings(self.new_molecule)
        self.new_molecule.UpdatePropertyCache(strict=False)
        Chem.rdmolfiles.MolToPDBFile(self.new_molecule, output_path)

        from utils import PDBTools
        PDBModifier = PDBTools()
        PDBModifier.rename_atoms_fragment(output_path)

