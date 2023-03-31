"""
This module has different functions to help prepare the system to run simulations for the combination AI-FragHop.
"""

import os
import logging


class AIHelpers(object): 

    @staticmethod
    def parse_set(path): 
        """
        It parses a set of smiles coming from an output of the generative models into a dictionary. 

        Parameters
        ----------
        path : str
            Path to the .txt/.csv of the set. 

        Returns
        -------
        d : dictionary
            Dictionary with the set of SMILES.
        """
        with open(path, 'r') as f:
            lines = f.readlines()
            if 'smiles' in lines[0]:
                d = {fields[0]:fields[1] for fields in [line.split(',') \
                    for line in lines[1:]]}
            else: 
                d = {i:line.split()[0] for i,line in enumerate(lines)}
        return d 

    @staticmethod
    def fragment_from_smiles(smiles, reference_ligand, output_PDB):
        """
        If a smiles with a different fragment is feched, it extracts the fragment and prepares a PDB to use FragHop with it.

        Parameters
        ----------
        smiles : str
            Smiles of the new ligand. 

        reference_ligand: str
            Path to the reference ligand PDB file.

        Returns
        -------
        fragment : a frag_hop.replacement.Fragment object
            Fragment extracted from SMILES.
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdFMCS
        
        def get_selected_bond(fragment):
            """
            It returns the selected bond that connects the new fragment with the scaffold. 

            Parameters
            ----------
            fragment : a rdkit.Molecule object.
                Selected fragment.
            """
            
            # Get dummy atom
            dummy_idx = next((a.GetIdx() for a in fragment.GetAtoms() if a.GetAtomicNum() == 0),None)
            
            # Get list of bonds
            bonds = [[bond.GetBeginAtom().GetIdx(),bond.GetEndAtom().GetIdx()] for bond in fragment.GetBonds()]

            #Get selected bond
            selected_bond = next((bond for bond in bonds if dummy_idx in bond), None)

            # Atom Idx to atom names
            PDBBlock = Chem.MolToPDBBlock(fragment)
            PDBinfo = [line.split() for line in PDBBlock.split('\n') if line.startswith('HETATM')]
            bond_atom_names = [line[2].replace('*', 'H') for line in PDBinfo 
                               if int(line[1]) - 1 in selected_bond]
            return '{}-{}'.format(bond_atom_names[0], bond_atom_names[1])

        # Reference ligand
        smi_original = Chem.MolToSmiles(Chem.MolFromPDBFile(reference_ligand))
        m_target = Chem.MolFromSmiles(smi_original)

        # New ligand
        m = Chem.MolFromSmiles(smiles)
        mcs = Chem.rdFMCS.FindMCS([m,m_target],
                                completeRingsOnly = True)
        
        match = m.GetSubstructMatch(Chem.MolFromSmarts(mcs.smartsString))
        frag = Chem.ReplaceCore(m, Chem.MolFromSmarts(mcs.smartsString))
        
        f = open('fragments.conf', 'a+') 
        # Check if the molecule is valid
        if not (len(match) == m_target.GetNumAtoms()): 
            raise AttributeError('The input structure differs from the reference structure with more than one fragment.', 
                + ' FragHop can not handle this cases for now.')
        
        # Generate the PDB of the fragment and its connectivity to the scaffold
        else:
            # Add Hydrogens 
            selected_bond = get_selected_bond(frag)
            
            # Replace the dummy atom with a hydrogen
            frag = AllChem.ReplaceSubstructs(frag,Chem.MolFromSmarts('[#0]'),
                Chem.MolFromSmarts('[#1]'))[0]

            # Correct and prepare fragment 
            frag.UpdatePropertyCache()
            frag = Chem.rdmolops.AddHs(frag)
            AllChem.EmbedMolecule(frag)
            
            # Generate PDB and save the connectivity
            Chem.MolToPDBFile(frag, output_PDB)

            # Load the Fragment 
            from frag_hop.replacement import Fragment
            fragment = Fragment(path_fragment = output_PDB, bonds = selected_bond)
                
        return fragment



