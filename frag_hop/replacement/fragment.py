"""
This module containts classes and methods involved in the manipulation of
fragments.
"""

# General imports
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem

from frag_hop.replacement.structure import InputStructure, _Structure
from frag_hop.utils.tools import RDKitTools

class Fragment(_Structure):
    """Class to prepare the selected fragment for replacement"""

    def __init__(self, path_ligand = None, path_fragment = None,
                 bonds = None, atom_names = None,
                 top_complex=None, top_fragment=None, resname = 'LIG'):
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
        super().__init__(terminal=True)

        # Scaffold
        self.path_fragment = path_fragment
        self.fragment = None

        # Ligand
        self.ligand = None
        self.ligand_prepared = None
        self.resname = resname

        if path_fragment:
            self.initialize_from_pdb(path_fragment, bonds)

        if path_ligand:
            if bonds is None and atom_names is None:
                raise AttributeError('To select a fragment from a ligand, ' +
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
        self.fragment = InputStructure(path, bonds_link=bonds)

    def initialize_from_ligand(self, path, bonds, atom_names):
        """
        It initialices a fragment as a InputStructure object from an input
        ligand.

        It extracts the fragment from the ligand and creates two structures, one
        containing the fragment and other containing the remaning scaffold of
        the ligand.

        Parameters
        ----------
        path : str
            Path to a PDB file containing a ligand.
        bonds : list
            List of bonds representing the connectivity points of the fragment.
        atom_names : list
            List of the atom names comforming the fragment.
        """

        def get_fragment_from_atom_names(molecule, atom_names):
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

        def get_fragment_from_bonds(molecule, bonds):
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
            bond = bonds[0] if len(bonds) == 1 else bonds

            # Breaks the molecule at the selected bonds
            idx1 = rdkit_tools.get_atomid_by_atomname(molecule, bond[1])
            idx2 = rdkit_tools.get_atomid_by_atomname(molecule, bond[0])

            Chem.Kekulize(molecule, clearAromaticFlags=True)
            em = Chem.EditableMol(molecule)
            em.RemoveBond(idx1, idx2)
            nm = em.GetMol()
            nm.GetAtomWithIdx(idx1).SetNoImplicit(True)
            nm.GetAtomWithIdx(idx2).SetNoImplicit(True)

            frags = Chem.GetMolFrags(nm, asMols=True, sanitizeFrags=False)

            # Extracts the fragment out of the ligand
            atom_fragment = bond[1].strip()
            for frag in frags:
                atom_names = [atom.GetPDBResidueInfo().GetName().strip() for
                              atom in frag.GetAtoms()]
                is_fragment = atom_fragment in atom_names
                if is_fragment:
                    with tempfile.NamedTemporaryFile(suffix='.pdb') as tmp:
                        Chem.rdmolfiles.MolToPDBFile(frag, tmp.name)
                        self.fragment = InputStructure(tmp.name,
                                                       bonds_link=bonds)
                else:
                    with tempfile.NamedTemporaryFile(suffix='.pdb') as tmp:
                        Chem.rdmolfiles.MolToPDBFile(frag, tmp.name)
                        self.ligand_prepared = InputStructure(tmp.name,
                                                              bonds_link=bonds)

        self.ligand = InputStructure(path, bonds_link=bonds)

        # Initializates the scaffold by the selected bonds
        if not bonds is None:
            get_fragment_from_bonds(self.ligand.rdkit_mol, bonds)

        # Initializates the scaffold by the selected atom names
        else:
            if not atom_names is None:
                get_fragment_from_atom_names(self.ligand.rdkit_mol, atom_names)

    def prepare(self, target):
        """
        It prepares the input fragment for later fragment replacement techniques
        to generate new molecules.
        """
        # Prepare fragment
        self.superimpose_fragment_bond(self.fragment, target.ligand,
                                       self.fragment.bonds_link,
                                       target.ligand.bonds_link)

        # Update RDKit molecule with the obtaind position
        with tempfile.NamedTemporaryFile(suffix='.pdb') as tmp:
            self.fragment.structure.save_pdb(tmp.name)
            self.fragment.rdkit_mol = \
                Chem.rdmolfiles.MolFromPDBFile(tmp.name, removeHs=False)
        ref = InputStructure(self.path_fragment)
        self.fragment.rdkit_mol = AllChem.AssignBondOrdersFromTemplate(
            ref.rdkit_mol,
            self.fragment.rdkit_mol)


        # Remove hydrogens of the scaffold
        self.remove_hydrogens(molecule = self.fragment)

    def to_file(self, path, file_name='fragment.pdb'):
        """
        Exports the PDB structure of the prepared fragment.

        Parameters
        ----------
        output_path : str
            Output path for the PDB structure.
        file_name : str
            Output PDB file name. Default: frag_prepared.pdb
        """
        output_path = os.path.join(path, file_name)
        Chem.rdmolfiles.MolToPDBFile(self.fragment.rdkit_mol,
                                     output_path)
        if self.ligand_prepared is not None:
            Chem.rdmolfiles.MolToPDBFile(self.ligand_prepared.rdkit_mol,
                                         os.path.join(path, 'lig_prepared.pdb'))
