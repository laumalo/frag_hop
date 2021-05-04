"""
This module containts classes that handle the fragment replacement in a complex.
"""

#General imports
from rdkit import Chem
import numpy as np
import os

#Local imports
from frag_hop.utils.geometry import rotation_matrix_axis
from frag_hop.utils.tools import RDKitTools

class Replacer:
    """Class to replace fragments"""

    def __init__(self, ligand_pdb, fragment_pdb, ref_fragment_pdb,
                 bond_atoms, resname, resnum,
                 out_folder='out', bond_type='single'):
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

        self.out = out_folder

        # Set bonds to break
        self.bond_frag, self.bond_lig = bond_atoms
        self.bond_type = bond_type

        # Ligand
        self.resname = resname
        self.resnum = resnum
        self.ligand = self.__load_to_rdkit(ligand_pdb)
        self.ligand_prepared = None
        self.original_fragment = None
        self.break_ligand()

        # Fragment
        self.fragment_path = fragment_pdb
        self.ref_fragment_path = ref_fragment_pdb
        self.fragment = self.__load_to_rdkit(fragment_pdb)
        self.correct_fragment()

        self.rotated_fragment = None
        self.get_best_dihedral_angle()

        # Merge structure
        self.merged = None
        self.get_merged_structure()

    def __load_to_rdkit(self, path):
        """
        It loads a PDB file into an rdkit.Chem.rdchem.Mol object.

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

        idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                  self.bond_lig[1])
        idx2 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                  self.bond_lig[0])
        lig = self.ligand
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
        Chem.rdmolfiles.MolToPDBFile(self.ligand_prepared,
                                     os.path.join(self.out,'lig_prepared.pdb'))
        Chem.rdmolfiles.MolToPDBFile(self.original_fragment,
                                     os.path.join(self.out,'original_frag.pdb'))

    def correct_fragment(self):
        """
        It assigns the bond orders of the prepared fragment from the reference
        fragment (from the library) and removes the hydrogen where the bond will
        take place.
        """

        from rdkit.Chem import AllChem
        from rdkit.Chem.rdmolops import FastFindRings

        # Assign bond orders
        ref = self.__load_to_rdkit(self.ref_fragment_path)
        Chem.rdmolfiles.MolToPDBFile(ref,'ref.pdb')
        Chem.rdmolfiles.MolToPDBFile(self.fragment,'frag.pdb')
        frag_bonds = AllChem.AssignBondOrdersFromTemplate(ref, self.fragment)

        # Remove hydrogen
        rdkit_tools = RDKitTools()
        frag_idx = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                      self.bond_frag[1])
        ed_frag = Chem.EditableMol(frag_bonds)
        ed_frag.RemoveAtom(frag_idx)

        # Save fragment
        mol = ed_frag.GetMol()
        FastFindRings(mol)
        mol.UpdatePropertyCache(strict=False)
        Chem.rdmolfiles.MolToPDBFile(mol, self.fragment_path)


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
        from frag_hop.data.parameters.atom_constants import BONDING_DISTANCES

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
        from rdkit.Chem.rdmolops import FastFindRings

        out_file = os.path.join(self.out, 'merged.pdb')
        out_ligand = os.path.join(self.out, resname + '.pdb')

        rdkit_tools = RDKitTools()

        for heteroatom, path in zip((True, False), (out_ligand,out_file)):
            if heteroatom is False:
                try:
                    em = Chem.EditableMol(self.merged)
                    hn_idx = rdkit_tools.get_atomid_by_atomname(self.merged, 'HN')
                    em.RemoveAtom(hn_idx)
                    self.merged = em.GetMol()
                    em = Chem.EditableMol(self.merged)
                    hxt_idx = rdkit_tools.get_atomid_by_atomname(self.merged, 'HXT')
                    em.RemoveAtom(hxt_idx)
                    self.merged = em.GetMol()
                except Exception:
                    pass
            for a in self.merged.GetAtoms():
                mi = Chem.AtomPDBResidueInfo()
                mi.SetName(a.GetPDBResidueInfo().GetName())
                mi.SetIsHeteroAtom(heteroatom)
                mi.SetResidueName(resname)
                mi.SetResidueNumber(resnum)
                mi.SetChainId(chain_id)
                a.SetMonomerInfo(mi)

            FastFindRings(self.merged)
            self.merged.UpdatePropertyCache(strict=False)
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
        self.__export_merged_structure(resname = self.resname,
                                       resnum=self.resnum)

    def get_best_dihedral_angle(self, rotation_angle=0.15):
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
                                                        self.bond_frag[0])
                p_ref = self.fragment.GetConformer().GetAtomPosition(ix)
                p_ref_rot = np.dot(rot_mat, [p_ref.x, p_ref.y, p_ref.z])
                return p_ref - p_ref_rot

            # Resets the rotated fragment
            self.rotated_fragment = self.__load_to_rdkit(self.fragment_path)

            # Computes vector to be used as axis of rotation
            idx1 = rdkit_tools.get_atomid_by_atomname(self.ligand,
                                                      self.bond_lig[1])
            idx2 = rdkit_tools.get_atomid_by_atomname(self.fragment,
                                                      self.bond_frag[0])

            p1 = self.ligand.GetConformer().GetAtomPosition(idx1)
            p2 = self.fragment.GetConformer().GetAtomPosition(idx2)

            u = np.array([p1.x, p1.y, p1.z]) - np.array([p2.x, p2.y, p2.z])
            u = u / (u**2).sum()**0.5

            # Gets rotation matrix
            rot_mat = rotation_matrix_axis(radi, u)

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

            from frag_hop.data.parameters.atom_constants import BONDING_DISTANCES
            from frag_hop.utils.geometry import distance

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
