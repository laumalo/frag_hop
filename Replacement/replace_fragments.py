import mdtraj as md
from lib_prep.FragmentTools import tree_detector


class InputStructure:
    """Object to modify scaffold and fragments structures"""
    
    def __init__(self, input_file, bond_link=None, top_file=None, chain="L", resnum=None):
        self.input_file = input_file
        self.top_file = top_file
        self.structure = None
        self.chain = chain
        self.resnum = resnum
        self.__load_to_mdtraj()
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

    def set_bond_to_link(self, new_bond):
        if type(new_bond) is list:
            if len(new_bond) == 2:
                self.bond_link = new_bond
            else:
                raise ValueError(f"Wrong length for {new_bond}. It must be a list of 2 elements")
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
        #To FILL
        pass

    def compute_fragment_rmsd(self):
        #To FILL
        pass
