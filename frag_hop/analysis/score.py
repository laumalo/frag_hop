"""
This module containts useful tools to analyze the results obtained for the
different hopping+PELE exploration options.
"""

import os
import re
import shutil
import logging


def score_simulation(folder, score_file=None):
    """
    It scores the results of a PELE simulation and extracts the best structure
    in a PDB file (top_structure.pdb) in the simulation folder.

    If a path for a scoring file is fetched, it will record the path to the best
    structure and its scoring.

    Parameters
    ----------
    folder : str
        Path to the folder to score the simulation.
    score_file : str
        Path to a file where the scoring will be recorded. Default: None
    """
    def get_best_structure(folder):
        """
        Given a folder that must contain inside a out_pele folder with the
        results of a PELE simulation, it gets the best structure and its score.

        Parameters
        ----------
        folder : str
            Path to the folder to analyze.

        Returns
        -------
        score : float
            Score of the best structure.
        struc_path : str
            Path to the best structure.
        """

        # Load files
        files = os.listdir(os.path.join(folder, 'out_pele'))
        reports_files = [i for i in files if i.startswith('report_')]

        # Metric to select best structure and best molecule
        idx_structure, idx_molecule = 0, 1

        d = {}
        for report in reports_files:
            report_path = os.path.join(folder, 'out_pele', report)
            with open(report_path, 'r') as f:
                d2 = {}
                for line in f.readlines()[1:]:
                    params = [x for x in line.split(' ') if not x == '']
                    d2[params[2]] = tuple((float(params[3]), float(params[4])))
                best_step = min(
                    d2.items(), key=lambda x: x[1][idx_structure])[0]
                best_energy = d2.get(best_step)
            d[report] = tuple((int(best_step), best_energy))
        best_traj = min(d.items(), key=lambda x: x[1][1][idx_molecule])[0]
        step, score = d.get(best_traj)

        traj = 'trajectory_{}.pdb'.format(best_traj.replace('report_', ''))
        traj_path = os.path.join(folder, 'out_pele', traj)
        model = step + 1
        trajectory_selected = re.search('MODEL\s+%d(.*?)ENDMDL' % model,
                                        open(traj_path, 'r').read(),
                                        re.DOTALL)
        traj = []
        struc_path = os.path.join(folder, 'top_structure.pdb')
        with open(struc_path, 'w') as f:
            traj.append("MODEL     %d" % model)
            traj.append(trajectory_selected.group(1))
            traj.append("ENDMDL\n")
            f.write("\n".join(traj))
        return score, struc_path

    # Gets score and best trajectory
    score, struc_path = get_best_structure(folder, extract=True)

    # Writes out the obtained information in a file
    if not score_file is None:
        with open(score_file, 'a+') as f:
            f.write('{} \t {} \n'.format(struc_path, score))


def score_library(multiple_folders_dir):
    """
    It scores the results of a full library.

    Parameters
    ----------
    multiple_folders_dir : str
        Path to the directory with all the simulations results.
    """
    MULTIPLE_FOLDERS = os.listdir(multiple_folders_dir)
    SCORE_FILE = os.path.join(multiple_folders_dir, 'results_hopping.txt')
    for FOLDER in MULTIPLE_FOLDERS:
        folder_path = os.path.join(multiple_folders_dir, FOLDER)
        try:
            score_simulation(folder_path, score_file=SCORE_FILE)
        except Exception:
            logging.warning('   - Skipping {}'.format(FOLDER))


def rank_scoring(score_file, top=10, extract=False):
    """
    It orders all the simulations in a score file and extract the top structures
    and their scorings.

    Parameters
    ----------
    score_file : str
        Path to the file with the scoring of the new ligand-protein complexes.
    top : int
        Number of the top structures you want to obtain. Default: 10.
    extract : bool
        Either you want to extract the top structures into a folder or not.
        Default: False
    """
    with open(score_file, 'r') as f:
        d = {line.split()[0]: float(line.split()[1]) for line in f.readlines()}
        sorted_d = sorted(d.items(), key=lambda kv: kv[1])
    with open(score_file.replace('.txt', '_top{}.txt'.format(top)), 'w') as f:
        for struc_path, score in sorted_d[:top]:
            f.write('{} \t {} \n'.format(struc_path, score))

    # Extract the top scored poses to a new folder
    if extract:
        folder_top = os.path.join(os.path.dirname(os.path.abspath(score_file)),
                                  'top_structures')
        os.makedirs(folder_top, exist_ok=True)
        for i, (struc_path, score) in enumerate(sorted_d[:top]):
            shutil.copy(struc_path, os.path.join(folder_top,
                                                 'top_{}.pdb'.format(i + 1)))
