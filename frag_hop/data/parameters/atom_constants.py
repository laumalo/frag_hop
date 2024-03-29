"""
Dictionary of bonding distances in Anstrongs

Info extracted from:
https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Chemical_Bonds/Bond_Lengths_and_Energies
"""
BONDING_DISTANCES = {
    # Atom1, Atom2, BondType: Distance(A)
    # Bonding with itself
    ("H", "H", "single"): 0.74,
    ("C", "C", "single"): 1.54,
    ("N", "N", "single"): 1.45,
    ("N", "S", "single"): 1.76,
    ("N", "O", "single"): 1.40,
    ("O", "O", "single"): 1.48,
    ("F", "F", "single"): 1.42,
    ("CL", "CL", "single"): 1.99,
    ("BR", "BR", "single"): 2.28,
    ("I", "I", "single"): 2.67,
    # Bonding with C
    ("C", "N", "single"): 1.47,
    ("C", "O", "single"): 1.43,
    ("C", "S", "single"): 1.82,
    ("C", "F", "single"): 1.35,
    ("C", "CL", "single"): 1.77,
    ("C", "BR", "single"): 1.94,
    ("C", "I", "single"): 2.18,
    ("N", "C", "single"): 1.47,
    ("O", "C", "single"): 1.43,
    ("S", "C", "single"): 1.82,
    ("F", "C", "single"): 1.35,
    ("CL", "C", "single"): 1.77,
    ("BR", "C", "single"): 1.94,
    ("I", "C", "single"): 2.18,
    # Other single bonds
    ("N", "S", "single"): 1.77,
    ("S", "N", "single"): 1.77,
    # Bonding with H
    ("H", "C", "single"): 1.08,
    ("H", "N", "single"): 1.01,
    ("H", "O", "single"): 0.98,
    ("H", "F", "single"): 0.99,
    ("H", "CL", "single"): 1.38,
    ("H", "BR", "single"): 1.53,
    ("H", "S", "single"): 1.34,
    ("H", "I", "single"): 1.73,
    ("C", "H", "single"): 1.08,
    ("N", "H", "single"): 1.01,
    ("O", "H", "single"): 0.98,
    ("F", "H", "single"): 0.99,
    ("CL", "H", "single"): 1.38,
    ("BR", "H", "single"): 1.53,
    ("I", "H", "single"): 1.73,
    ("S", "H", "single"): 1.34,
    # Double bonds
    ("C", "C", "double"): 1.34,
    ("O", "O", "double"): 1.21,
    ("C", "O", "double"): 1.23,
    ("O", "C", "double"): 1.23,
    ("S", "O", "double"): 1.48,
    ("O", "S", "double"): 1.48,
    ("C", "N", "double"): 1.28,
    ("N", "C", "double"): 1.28,
    ("S", "C", "double"): 1.66,
    ("C", "S", "double"): 1.66,
    # Triple bonds
    ("N", "N", "triple"): 1.1,
    ("C", "N", "triple"): 1.14,
    ("N", "C", "triple"): 1.14
}
