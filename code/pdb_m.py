import numpy as np
import matplotlib.pyplot as plt
import nglview as nv
import pandas as pd
import math
import colorsys
import matplotlib as mpl
import argparse
import random
import seaborn as sns
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from scipy.optimize import minimize


def pdb_2_df(pdb_file):
    # converting the pdb file to a dataframe

    params = [[] for _ in range(12)]

    with open(pdb_file) as f:
        for line in f:
            cols = line.split()
            if len(cols) == 12:
                for i in range(12):
                    if cols[0] == 'ATOM':
                        params[i].append(cols[i])
    new_column_names = ["Type", "Atom Number", "Atom Type", "Residue Name",
                        "Chain ID", "Residue Number", "x", "y", "z", "Occupancy",
                        "Mass", "Chem Atom"]
    # print(params)
    df = pd.DataFrame(np.fliplr(params))
    df = df.T
    df.columns = new_column_names
    return df  # think whether you need the sequence of dataframe ....

df = pdb_2_df('tau.pdb')
df['Residue Number'] = df['Residue Number'].astype(int)  # Convert the column to integer type
df.sort_values(by='Residue Number')

import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Custom dictionary for converting three-letter codes to one-letter codes
three_to_one = {
    'ALA': 'A',
    'CYS': 'C',
    'ASP': 'D',
    'GLU': 'E',
    'PHE': 'F',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LYS': 'K',
    'LEU': 'L',
    'MET': 'M',
    'ASN': 'N',
    'PRO': 'P',
    'GLN': 'Q',
    'ARG': 'R',
    'SER': 'S',
    'THR': 'T',
    'VAL': 'V',
    'TRP': 'W',
    'TYR': 'Y'
}

def calculate_center_of_mass(atoms):
    total_mass = 0.0
    center_of_mass = [0.0, 0.0, 0.0]

    for _, atom in atoms.iterrows():
        mass = float(atom['Mass'])
        coord = atom[['x', 'y', 'z']]
        total_mass += mass
        center_of_mass[0] += mass * float(coord[0])
        center_of_mass[1] += mass * float(coord[1])
        center_of_mass[2] += mass * float(coord[2])

    center_of_mass[0] /= total_mass
    center_of_mass[1] /= total_mass
    center_of_mass[2] /= total_mass

    return center_of_mass

def calculate_center_of_mass_for_residue_numbers(df):
    result = pd.DataFrame(columns=['Residue Number', 'Residue Name', 'Center of Mass','Mass'])
    grouped = df.groupby('Residue Number')

    for residue_number, residue_atoms in grouped:
        center_of_mass = calculate_center_of_mass(residue_atoms)
        residue_name = three_to_one[residue_atoms['Residue Name'].iloc[0]]
        mass = float(residue_atoms['Mass'].iloc[0])
        result = result.append({'Residue Number': residue_number, 'Residue Name': residue_name, 'Center of Mass': center_of_mass,'Mass': mass}, ignore_index=True)
    result['Residue Number'] = result['Residue Number'].astype(int)
    result = result.sort_values(by = 'Residue Number')
    #df['Residue Number'] = df['Residue Number'].astype(int)  # Convert the column to integer type
    #df.sort_values(by='Residue Number')
    return result

# Assuming you have a DataFrame named 'df' representing the PDB file
com_df = calculate_center_of_mass_for_residue_numbers(df)
calculate_center_of_mass_for_residue_numbers(df)
com_df[['x', 'y', 'z']] = pd.DataFrame(com_df['Center of Mass'].tolist(), index=com_df.index)


def com(amino_acids_df):
    total_mass = com_df['Mass'].sum()
    center_of_mass_x = np.sum(amino_acids_df['x'] * amino_acids_df['Mass']) / total_mass
    center_of_mass_y = np.sum(amino_acids_df['y'] * amino_acids_df['Mass']) / total_mass
    center_of_mass_z = np.sum(amino_acids_df['z'] * amino_acids_df['Mass']) / total_mass

    center_of_mass = (center_of_mass_x, center_of_mass_y, center_of_mass_z)
    return center_of_mass


def find_v(amino_acids_df):
    total_mass = com_df['Mass'].sum()
    amino_acids_df[['x', 'y', 'z']] = pd.DataFrame(amino_acids_df['Center of Mass'].tolist(),
                                                   index=amino_acids_df.index)

    center_of_mass_x, center_of_mass_y, center_of_mass_z = com(amino_acids_df)

    # Calculate vectors from the center of mass to each amino acid

    amino_acids_df['vector_x'] = amino_acids_df['x'] - center_of_mass_x
    amino_acids_df['vector_y'] = amino_acids_df['y'] - center_of_mass_y
    amino_acids_df['vector_z'] = amino_acids_df['z'] - center_of_mass_z

    return amino_acids_df


find_v(com_df)


def pair_connectivity(matrix):
    dir = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    pair_lst = []
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] == 0:
                continue
            for k in range(4):
                if 0 <= i + dir[k][0] < len(matrix) and 0 <= j + dir[k][1] < len(matrix[i]):
                    if matrix[i][j] + 4 <= matrix[i + dir[k][0]][j + dir[k][1]] or matrix[i][j] - 4 >= \
                            matrix[i + dir[k][0]][j + dir[k][1]]:
                        pair_lst.append((min(matrix[i][j], matrix[i + dir[k][0]][j + dir[k][1]]),
                                         max(matrix[i][j], matrix[i + dir[k][0]][j + dir[k][1]])))

    p_lst = []
    p_lst.append(pair_lst[0])
    for row in pair_lst:
        if row not in p_lst:
            p_lst.append(row)

    return p_lst

# Constants
NUM_PARTICLES = len()
DIMENSIONS = 3
CLOSE_RADIUS = 0.4

# Define the connectivity of the polymer as a list of tuples (particle1, particle2)
connectivity = pair_connectivity(mat)
#define connectivity

# Define the initial positions of the particles as a 1D array [x1, y1, z1, x2, y2, z2, ...]
initial_positions = np.random.rand(NUM_PARTICLES * DIMENSIONS)
#define this

def distance(p1, p2):
    # Calculate the Euclidean distance between two particles given their positions
    return np.linalg.norm(p1 - p2)

def energy_function(positions):
    # Calculate the total energy of the system based on the positions of the particles

    # Reshape the 1D array of positions into a 2D array of shape (NUM_PARTICLES, DIMENSIONS)
    positions = positions.reshape(NUM_PARTICLES, DIMENSIONS)

    total_energy = 0

    for particle1, particle2 in connectivity:
        # Calculate the distance between the connected particles
        dist = distance(positions[particle1], positions[particle2])
        # Add a penalty to the energy if the distance is not satisfied
        total_energy += (dist - CLOSE_RADIUS)**2

    return total_energy

# Constraints for the optimization (keep the particles connected)
constraints = {'type': 'eq', 'fun': lambda positions: np.array([distance(positions[particle1], positions[particle2]) - CLOSE_RADIUS for particle1, particle2 in connectivity])}

# Run the optimization
result = minimize(energy_function, initial_positions, constraints=constraints)

# Get the optimized positions
optimized_positions = result.x.reshape(NUM_PARTICLES, DIMENSIONS)
initial_positions =  initial_positions.reshape(NUM_PARTICLES, DIMENSIONS)