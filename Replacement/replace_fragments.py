import mdtraj as md
from lib_prep.FragmentTools import tree_detector
from rdkit import Chem

class InputStructure:
    """Object to modify scaffold and fragments structures"""
    
    def __init__(self, input_file, bond_link=None, top_file=None, chain="L", resnum=None):
        self.input_file = input_file
        self.top_file = top_file
        self.structure = None
        self.chain = chain
        self.resnum = resnum
        self.__load_to_mdtraj()
        self.rdkit_representation = None
        self.__load_to_rdkit()

        if bond_link:
            self.set_bond_to_link(bond_link)
        else:
            self.bond_link = bond_link

    # Add any attibute or method that you consider

    def __load_to_mdtraj(self):
        if self.top_file:
            self.structure = md.load(self.input_file, top_file=self.top_file)
        else:
            self.structure = md.load(self.input_file)

    def __load_to_rdkit(self):
        # Depending on the metric we might need the rdkit representation
        self.rdkit_representation = Chem.rdmolfiles.MolFromPDBFile(
                                            self.input_file, removeHs=False)

    def set_bond_to_link(self, new_bond):
        if type(new_bond) is list:
            if len(new_bond) == 2:
                self.bond_link = new_bond
            else:
                raise ValueError(
                    "Wrong length for {new_bond}. It must be a list of 2 elements")
        else:
            raise TypeError(f"Wrong type for {new_bond}. It must be a 'list'")


    def get_atoms_to_delete(self):
        if self.bond_link:
            atoms_to_delete = tree_detector.main(self.input_file, self.bond_link,
                                                 chain=self.chain, resnum=self.resnum)
        else:
            raise ValueError("Non-defined bond to link. Set it before deleting atoms!")
        return atoms_to_delete
    
class Replacer:
    """Class to replace fragments"""
    def __init__(self, initial_complex, fragment, ligand_resnum, 
                 top_complex=None, top_fragment=None, chain_complex="A", 
                 chain_fragment="L", bond_type='single'):
        self.initial_complex = InputStructure(initial_complex, top_file=top_complex,
                                              chain=chain_complex, resnum=ligand_resnum)
        self.fragment = InputStructure(fragment, top_file=top_fragment,
                                       chain=chain_fragment)
        self.combination = None
        self.bond_type = bond_type

    # Add any attibute or method that you consider

    # Example of methods
    def set_complex_link(self, bond):
        self.initial_complex.set_bond_to_link(bond)

    def set_fragment_link(self, bond):
        self.fragment.set_bond_to_link(bond)

    def change_fragment(self, fragment_file, fragment_top=None):
        self.fragment = InputStructure(fragment_file, top_file=fragment_top)

    def get_atoms_to_delete(self):
        atoms_core = self.initial_complex.get_atoms_to_delete()
        atoms_fragment = self.fragment.get_atoms_to_delete()
        d = { 'core' : atoms_core,
              'fragment' : atoms_fragment }
        return d

    def superimpose_fragment_bond(self):
        # To FILL
        pass

    def rotate_fragment(self):
        # To FILL
        pass

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

            deviation = sum( (a-b)**2 for a, b in zip(USR1, USR2) )
            return 1/(1 + (1/12)*deviation)

        from rdkit.Chem import rdMolDescriptors

        # Shape similarity coefficient
        usr1 = rdMolDescriptors.GetUSR(rdkit_molA)
        usr2 = rdMolDescriptors.GetUSR(rdkit_molB)
        return s(usr1,usr2)
