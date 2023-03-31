"""
This module contains all the funcions needed to handle protein-ligand complexes
with non-covalent interactions.
"""
def extract_chain(chain_id, input_file, output_file):
    """
    It extracts a chain from a PDB file and generates a new PDB file only
    containing this chain.

    Parameters
    ----------
    chain_id : str
        Chain ID.
    input_file : str
        Path to the PDB file to extract the chain from.
    output_file : str
        Path to the new PDB file to save the extracted chain.
    """
    def select_chain(fhandle, chain_set):
        """
        Filters the PDB file for specific chain identifiers.
        """
        records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
        for line in fhandle:
            if line.startswith(records):
                if line[21] not in chain_set:
                    continue
            yield line

    pdbfh = open(input_file, 'r')
    new_pdb = select_chain(pdbfh, chain_id)
    outf = open(output_file, 'w')
    for line in new_pdb:
        outf.write(line)

def perform_chain_substitution(complex_pdb, new_chain_pdb, new_complex_pdb,
                               chain_id):
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
    f2 = open(new_chain_pdb, 'r')
    f3 = open(new_complex_pdb, 'w')
    lines_to_delete, new_lines = ([] for x in range(2))

    # Read original complex
    lines = f1.readlines()
    atoms_to_delete = [line[6:11] for line in lines if line[21:22]==chain_id]

    protein_pdb = [line for line in lines if line.startswith('ATOM')]
    protein_connect_pdb = [line for line in lines if line.startswith('CONECT')
                           and not line[6:11] in atoms_to_delete]




    # Load new residue atoms and correct atom numbers format
    new_lines = [line for line in f2.readlines() if 'HETATM' in line
    or 'CONECT' in line]
    new_lines = [line[:6] + str(int(line[6:11])).rjust(5) + line[11:]
                 for line in new_lines]
    # Generate output file for the PDB structure
    new_file_lines = \
        protein_pdb + new_lines + \
        protein_connect_pdb

    for line in new_file_lines:
        f3.write(line)

def get_resname_resnum_from_chain(path, chain_id):
    """
    It gets the residue name and the residue sequence number of the selected
    chain. COnsider that this is design for ligands where the chain contains a
    single non-standard residue.

    Parameters
    ----------
    path : str
        The path to the protein-ligand complex PDB file.
    chain_id : str
        Chain identifier.

    Returns
    -------
    resname : str
        Residue name.
    resnum : str
        Residue sequence number.
    """
    with open(path, 'r') as f:
        lines = f.readlines()
        resname = [line[17:20] for line in lines if line[21:22]==chain_id][0]
        resnum = [line[22:26] for line in lines if line[21:22]==chain_id][0]
        return resname, resnum.strip()
