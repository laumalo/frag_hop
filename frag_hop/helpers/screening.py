"""
This module containts all the workflows to run the different library screening
protocols.
"""
import os
import logging
import numpy as np
import time

from rdkit import Chem


class Screening():
    """Class to perform the different methods for screening a library."""

    def __init__(self, complex_pdb, lib_path, conf_file,
                 core, covalent, output='out'):
        """
        It initialices a Screening object

        Parameters
        ----------
        complex_odb : str
            Path to the target complex PDB file.
        lib_path : str
            Path to the fragments library.
        conf_file : str
            Path to the configuration file for connectivities of the fragments.
        """
        self.complex_pdb = complex_pdb
        self.lib_path = lib_path
        self.conf_file = conf_file
        self.output = output

    def run(self, type_screening='basic', CONTROL_FILE=None, PELE_exc=None):
        """
        It runs the selected screening method.

        Parameters
        ----------
        type_screening : str
            Screening method.
        CONTROL_FILE : str
            Path to a PELE control file.
        PELE_exc : str
            Path to a PELE executable version.
        """
        if CONTROL_FILE is None or PELE_exc is None:
            logging.warning('You have not selected a control file or a PELE' +
                            ' executable version, only the replacement of' +
                            ' fragments (witout refinement) will be performed.')
            self.basic_screening(CONTROL_FILE, PELE_exc)

        TYPES_AVAILABLE = ['basic', 'sequential', 'cluster']
        if type_screening in TYPES_AVAILABLE:
            if type_screening == 'basic':
                self.basic_screening(CONTROL_FILE, PELE_exc)
            if type_screening == 'sequential':
                self.sequential_screening(CONTROL_FILE, PELE_exc)
            if type_screening == 'cluster':
                self.cluster_screening(CONTROL_FILE, PELE_exc)
        else:
            raise AttributeError('Selected screening method {} is invalid.'
                                 .format(type_screening))

    def _align_molecules(self, mols):
        """
        It aligns all the molecules of a set based on Crippen logP atom
        contributions.

        Parameters
        ----------
        mols : list[rdkit.molecule object]
            Set of molecules.

        Returns
        -------
        aligned_mols : list[rdkit.molecule object]
            Aligned set of molecules.
        """
        from rdkit.Chem import rdMolDescriptors, rdMolAlign
        crippen_contribs = [rdMolDescriptors._CalcCrippenContribs(mol)
                            for mol in mols]
        crippen_ref_contrib, crippen_prob_contribs = \
            crippen_contribs[0], crippen_contribs[1:]
        ref_mol, prob_mols = mols[0], mols[1:]

        for i, mol in enumerate(prob_mols):
            crippenO3A = rdMolAlign.GetCrippenO3A(mol, ref_mol,
                                                  crippen_prob_contribs[i],
                                                  crippen_ref_contrib,)
            crippenO3A.Align()
        prob_mols.append(ref_mol)
        return prob_mols

    def _dist_matrix(self, mols):
        """
        It returns the distance matrix of a list of molecules using the Shape
        Tanimoto distance.

        Parameters
        ----------
        mols : list[rdkit.molecule object]
            Molecules

        Returns
        -------
        dists : np.array
            Distance matrix
        """
        from rdkit.Chem import rdShapeHelpers
        mols_aligned = self._align_molecules(mols)
        dists = [[rdShapeHelpers.ShapeTanimotoDist(ref_mol, mol)
                  for mol in mols_aligned] for ref_mol in mols_aligned]
        return np.array(dists)

    def _dist_vector(self, ref_mol, mols):
        """
        It returns the distance vector of a list of molecules using the Shape
        Tanimoto distance to a reference molecule.

        Parameters
        ----------
        ref_mol : a rdkit.molecule object
            Reference molecule
        mols : list[rdkit.molecule object]
            Molecules

        Returns
        -------
        sims : np.array
            Distance vector
        """
        from rdkit.Chem import rdShapeHelpers
        mols_aligned = self._align_molecules(mols)
        sims = [rdShapeHelpers.ShapeTanimotoDist(ref_mol, mol)
                for mol in mols_aligned]
        return np.array(sims)

    def basic_screening(self, CONTROL_FILE, PELE_exc):
        """
        This protocol is the most basic way of screening the library, each
        fragment is replaced using a hopping tecnique to original complex and a
        short PELE simulation is runned in order refine and rescore the new
        molecule.

        Parameters
        ----------
        CONTROL_FILE : str
            Path to a PELE control file.
        PELE_exc : str
            Path to a PELE executable version.
        """

        # Perform replacement
        fragments = [f for f in os.listdir(
            self.lib_path) if f.endswith('.pdb')]
        for i, fragment in enumerate(fragments):
            folder_name = os.path.join(
                self.output, fragment.replace('.pdb', ''))
            logging.info(' -   Replacing fragment {}:{}'.format(i, fragment))
            try:
                from frag_hop.helpers.hopping import HoppingSelector
                hopping_method = HoppingSelector(core=self.core,
                                                 covalent=self.covalent)
                hopping_method.run(complex_pdb=self.complex_pdb,
                                   hit_pdb=os.path.join(self.lib_path,
                                                        fragment),
                                   conf_file=self.conf_file,
                                   id_lig='GRW' if self.covalent else 'L',
                                   output=folder_name)
            except:
                logging.warning(
                    '     -   Skipping fragment {}'.format(fragment))

        if CONTROL_FILE is not None and PELE_exc is not None:
            # Prepare and run short PELE exploration
            # TODO
            # This loop could be done in parallel?? + merge with the other loop
            for fragment in enumerate(fragments):
                folder_name = os.path.join(
                    self.output, fragment.replace('.pdb', ''))

                from frag_hop.helpers.pele import generate_templates
                generate_templates(FOLDER_TO_PREPARE=folder_name)

                from frag_hop.helpers.pele import PELErunner
                PELErunner(control_file=CONTROL_FILE, PELE_exc=PELE_exc)


    def sequential_screening(self, CONTROL_FILE, PELE_exc):
        """
        This protocol for screening a fragment library consists on ranking the
        fragments accoding to shape similarity. So when a fragment is replaced,
        the cavity was adapted to the most similiar fragmnt thus the difference
        is the minimum possible.

        Parameters
        ----------
        CONTROL_FILE : str
            Path to a PELE control file.
        PELE_exc : str
            Path to a PELE executable version.
        """

        def get_next_fragment(ref_mol, mols):
            """
            Given a reference molecule and a list of molecules, it returns the
            most similiar molecule of the list to the refenrece molecule based
            on the Shape Tanimoto distance.

            It also removes this molecule from the list and adds its label to
            the ranking list.

            Parameters
            ----------
            ref_mol : a rdkit.molecule object
                Reference molecule
            mols : list[rdkit.molecule object]
                List of molecules
            """

            # Compute the distance from the ref mol to all the other molecules
            distances = self._dist_vector(ref_mol, mols)
            d_scores = dict(zip(distances, fragments))

            # Get the most similiar molecule from the list
            top_fragment = d_scores.get(min(distances))
            idx = fragments.index(top_fragment)

            # Remove this molecule and append to the ranking
            fragments.pop(idx), mols.pop(idx)
            ranking.append(top_fragment)
            return top_fragment

        fragments = [f for f in os.listdir(
            self.lib_path) if f.endswith('.pdb')]
        mols = [Chem.MolFromPDBFile(os.path.join(self.lib_path, f))
                for f in fragments]
        d_frags = dict(zip(fragments, mols))
        ref_fragment = '/Users/laura/LAURA/BSC/SARA/FragsLib/S1/S1_100.pdb'

        # Get the most similiar fragment from the library to the target fragment
        ranking = []
        ref_mol = Chem.MolFromPDBFile(ref_fragment)
        top_fragment = get_next_fragment(ref_mol, mols)

        # Order the other fragments of the library based on similarity
        while fragments:
            ref_mol = Chem.MolFromPDBFile(os.path.join(self.lib_path,
                                                       top_fragment))
            top_fragment = get_next_fragment(ref_mol, mols)

        # Sort the fragments according to this ranking
        index_map = {v: i for i, v in enumerate(ranking)}
        d = dict(sorted(d_frags.items(), key=lambda pair: index_map[pair[0]]))

        target_complex_path = self.complex_pdb
        for k, v in d.items():

            hit_fragment_path = os.path.join(self.lib_path, k)
            folder_name = os.path.join(
                    self.output, k.replace('.pdb', ''))
            # Run hopping
            try:

                from frag_hop.helpers.hopping import HoppingSelector
                hopping_method = HoppingSelector(core=self.core,
                                                 covalent=self.covalent)
                hopping_method.run(complex_pdb=target_complex_path,
                                   hit_pdb=hit_fragment_path,
                                   conf_file=self.conf_file,
                                   id_lig='GRW' if self.covalent else 'L',
                                   output=folder_name)

                # Prepare folders
                from frag_hop.helpers.pele import generate_templates
                generate_templates(FOLDER_TO_PREPARE=folder_name)

                from frag_hop.helpers.pele import PELErunner
                PELErunner(control_file=CONTROL_FILE, PELE_exc=PELE_exc)

                # Get top structure for next initial complex path
                target_complex_path = os.path.join(folder_name,
                                                   'top_structure.pdb')

                # Wait until the PELE simulation has finisehd
                condition = os.path.isfile(target_complex_path)
                while not condition:
                    condition = os.path.isfile(target_complex_path)
                    time.sleep(0.1)

            except:
                logging.warning(
                    '     -   Skipping fragment {}'.format(k))


    def cluster_screening(self, CONTROL_FILE, PELE_exc, n_clusters = 10):
        """
        This protocol for screening clusters the library and strarts by growing
        one fragment of each cluster using FragPELE, then this grown structure
        is used as initial structure for the fragment hopping within each
        cluster.

        Parameters
        ----------
        CONTROL_FILE : str
            Path to a PELE control file.
        PELE_exc : str
            Path to a PELE executable version.
        n_clusters : int
            Number of clusters. Default : 10
        """

        # Generate distance matrix
        fragments = [f for f in os.listdir(
            self.lib_path) if f.endswith('.pdb')]
        mols = [Chem.MolFromPDBFile(os.path.join(self.lib_path, f))
                for f in fragments]
        dists = self._dist_matrix(mols)

        # Clustering algorithm
        # TODO: find a better clustering method
        from sklearn.preprocessing import normalize
        normed_dists = normalize(dists, axis=1, norm='l1')

        from sklearn.cluster import AgglomerativeClustering
        clustering = AgglomerativeClustering(n_clusters=20,
                                             affinity='precomputed',
                                             linkage='average')
        clustering.fit(normed_dists)

        # Dicitonary with the population of each cluster
        from collections import defaultdict
        d = defaultdict(list)
        for l, f in zip(clustering.labels_, fragments):
            d.setdefault(l, []).append(f)

        # Get one representative fragment of each cluster and grow the fragment
        # with FragPELE

        # Replace the other fragments to their representative of the fragment
        # using FragHop techniques
