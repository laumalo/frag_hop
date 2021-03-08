from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import os
import glob
import pandas as pd


class Ranking:
    """
    Class that contains all the tools to rank the fragments according to
    similarity using the Tanimoto coefficient.
    """

    def __init__(self, original_fragment_path, library_path):
        """
        It initializes a Ranking object through the path of the reference
        fragment and the full library of fragments.

        Parameters
        ----------
        original_fragment_path : str
            Path to the original fragment of the ligand
        library_path : str
            Path to the library containing the replacement fragment candidates.

        Example
        -------

        Generates a CSV file with the fragment ranked acoording to similarity.

        >>> ranking_tool = Ranking(path/to/fragment.pdb, path/to/library)
        >>> ranking_tool.rank_fragments(to_file = True)

        """
        self.original_fragment = Chem.rdmolfiles.MolFromPDBFile(
            original_fragment_path, removeHs=False)
        self.library_path = library_path
        self.frag_names = None
        self.frags = None

    def _load_all_candidates(self):
        """
        It loads all the fragment candidates to be ranked.
        """

        path = os.path.join(self.library_path, '*pdb')
        fragment_files = glob.glob(path)

        self.frag_names = fragment_files
        self.frags = [Chem.rdmolfiles.MolFromPDBFile(frag, removeHs=False)
                      for frag in fragment_files]

    def _compute_Tanimoto_coeff(self, mol1, mol2):
        """
        It computes the Tanimoto coefficient between two molecules.
        """
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=1024)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=1024)
        tm_coeff = DataStructs.FingerprintSimilarity(fp1, fp2)
        return tm_coeff

    def _diff_Tanimoto(self, main_mol, mols):
        """
        It computes the Tanimoto coefficient from all the fragments in the
        library against the reference fragment.
        """
        tm_coeffs = [self._compute_Tanimoto_coeff(main_mol, mol)
                     for mol in mols]
        return tm_coeffs

    def rank_fragments(self, to_file=False, output_file='out.csv'):
        """
        It ranks all the fragments of the library according to their similarity
        to the reference fragment.
        """

        self._load_all_candidates()
        tm_coeffs = self._diff_Tanimoto(self.original_fragment, self.frags)
        sorted_frags = sorted(list(zip(self.frag_names, tm_coeffs)),
                              key=lambda x: x[1], reverse=True)
        if to_file:
            df = pd.DataFrame(sorted_frags, columns=['Path', 'Tm coeff'])
            df.to_csv(output_file)
        return sorted_frags
