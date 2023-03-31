"""
This module contains all the funcions needed to handle protein-ligand complexes
covalently bound.
"""


def extract_residue(input_pdb, resname, output_pdb):
    """
    It extracts one residue of a PDB file into a separate PDB file.

    Parameters
    ----------
    input_pdb : str
        Path to the input PDB file.
    resname : str
        Name of the residue to extract.
    output_pdb : str
        Path to the output PDB file.
    """
    f1 = open(input_pdb, 'r')
    f2 = open(output_pdb, 'w')
    residue = [line for line in f1.readlines() if resname in line]

    for line in residue:
        f2.write(line)


def perform_residue_substitution(complex_pdb, new_residue_pdb, new_complex_pdb,
                                 resname):
    """
    Given a PDB with a complex and the PDB of the residue that has to be
    replaced. It replaces from the original complex the residue and exports
    the PDB of the new strucutre.

    Parameters
    ----------
    complex_pdb : str
        Path to the initial complex PDB.
    new_residue_pdb : str
        Path to the residue PDB you want to modidy in the complex.
    new_complex_pdb : str
        Path to the output complex PDB with the residue modified.
    resname : str
        Residue name for the modified residue.
    """

    f1 = open(complex_pdb, 'r')
    f2 = open(new_residue_pdb, 'r')
    f3 = open(new_complex_pdb, 'w')
    lines_to_delete, new_lines = ([] for x in range(2))

    # Read original complex
    lines = f1.readlines()
    lines_to_delete = [line_num for line_num, line in enumerate(lines)
                       if resname in line]

    # Load new residue atoms and correct atom numbers format
    new_lines = [line for line in f2.readlines() if 'ATOM' in line]
    new_lines = [line[:6] + str(int(line[6:11])).rjust(5) + line[11:]
                 for line in new_lines]

    # Generate output file for the PDB structure
    new_file_lines = \
        lines[:lines_to_delete[0]] + new_lines + \
        lines[lines_to_delete[-1] + 1:]

    atom_ids = [int(line[6:11].strip()) for line in new_file_lines
                if line.startswith('ATOM') or line.startswith('TER')
                and line[6:11].strip() != '']

    # Check if the PDB atom ids are unique, if not correct them
    if len(atom_ids) == len(set(atom_ids)):
        for line in new_file_lines:
            f3.write(line)
    else:
        max_id = max(list(set(atom_ids)))
        corrected_lines = []
        for line in new_file_lines:
            if resname in line[17:20]:
                new_atom_id = str(int(line[6:11].strip()) + max_id).rjust(5)
                corrected_lines.append(
                    line[0:6] + new_atom_id + line[11:80] + '\n')
            else:
                corrected_lines.append(line)
        for line in corrected_lines:
            f3.write(line)


def get_resnum_from_resname(path, resname):
    """
    It gets the residue sequence number from the residue name in a PDb file.

    Parameters
    ----------
    path : str
        The path to the PDB file.
    resname : str
        Residue name

    Returns
    -------
    resnum :
        Residue sequence number
    """
    with open(path, 'r') as f:
        lines = f.readlines()
        resnum = [line[22:26] for line in lines if line[17:20] == resname][0]
        return resnum.strip()
