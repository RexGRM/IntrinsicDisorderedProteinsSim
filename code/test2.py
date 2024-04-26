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
import json
from scipy.optimize import minimize
#pdb_file path#
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


'''task number 2'''

def df_2_order(df):
    df['Atom Number'] = df['Atom Number'].astype(int)  # Convert the column to integer type
    return df.sort_values(by='Atom Number')


def seq_2_fasta(seq):
    '''
    input: list of sequence of the amino acid
    output: FASTA file

    Description:converts the sequence to a Fasta file

    '''
    # fasta......
    with open('pdb_file' + ".fasta", "w") as f:
        f.write(">" + 'pdb_file' + "\n")
        f.write(seq + "\n")


def df_to_sequence(df):
    '''
    Input: df, type=dataframe
    output: sequence, type=str

    Description:
    convert the aminoacids in PDB dataframe from 3-letter code to 1-letter code


    '''
    new_df = df.drop_duplicates(subset='Residue Number')
    residues = new_df['Residue Name'].tolist()

    one_letter_code = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y"
    }

    sequence = "".join([one_letter_code[residue] if residue in one_letter_code else "X" for residue in residues])

    return sequence

#ph#
def seq_to_charge(seq, ph):
    """

    input:the sequence of amino acid
    output:an array of the charges of the sequence
    description: record the charges of the sequence and stores it in an array

    """
    charge_seq = []

    for i in range(len(seq)):
        charge_seq.append(get_charge(seq[i], ph))
    return charge_seq


'''task number 4'''

#ph#
def get_charge(s, ph):
    '''

    input:the amino acid residue;the ph of environment
    output: the charge of the amino acid
    description: calculate the charge of the corresponding amino acid
    '''
    # pka = pd.read_csv("pka_pi_aa.csv")
    # convert all the amino acids to charges

    # pka_dict = df.set_index('Abr').to_dict()['pI']
    # this might need testing

    pka_dict = {
        "A": 6.01,
        "R": 10.76,
        "N": 5.41,
        "D": 2.77,
        "C": 5.07,
        "Q": 5.65,
        "E": 3.22,
        "G": 6.06,
        "H": 7.64,
        "I": 6.05,
        "L": 6.01,
        "K": 9.74,
        "M": 5.74,
        "F": 5.48,
        "P": 6.30,
        "S": 5.68,
        "T": 5.60,
        "W": 5.89,
        "Y": 5.66,
        "V": 5.96
    }

    if pka_dict[s] > ph:
        return 1 * ((10 ** (ph - pka_dict[s])) / ((10 ** (ph - pka_dict[s])) + 1))
    else:
        return -1 * ((10 ** (ph - pka_dict[s])) / ((10 ** (ph - pka_dict[s])) + 1))

# you need to output all the amino acids using their PIs instead of just the charges

def seq_disorder(seq):
    disorder_dict = {'G': 0.166, 'A': 0.06, 'V': -0.121, 'L': -0.326, 'I': -0.486, 'F': -0.697, 'Y': -0.51, 'W': -0.884,
                     'M': -0.397, 'P': 0.987,
                     'C': 0.02, 'S': 0.341, 'T': 0.059, 'N': 0.007, 'Q': 0.318, 'D': 0.192, 'E': 0.736, 'K': 0.586,
                     'R': 0.180, 'H': 0.303}
    disorder_seq = []
    for i in range(len(seq)):
        disorder = disorder_dict[seq[i]]
        disorder_seq.append(disorder)

    return disorder_seq

def seq_hydrophobicity(seq):
    hydrophobic_dict = {'G': 1.14, 'A': 0.33, 'V': -0.53, 'L': -0.69, 'I': -0.81, 'F': -0.58, 'Y': 0.23, 'W': -0.24,
                        'M': -0.44, 'P': -0.31,
                        'C': 0.22, 'S': 0.33, 'T': 0.11, 'N': 0.19, 'Q': 0.43, 'D': 2.41, 'E': 1.61, 'K': 1.81,
                        'R': 1.0, 'H': 1.37}
    hydro_seq = []
    for i in range(len(seq)):
        hydro = hydrophobic_dict[seq[i]]
        hydro_seq.append(hydro)

    return hydro_seq


def initialize_polymer(n, m, charges, hydro, poly_length):
    matrix = np.zeros((n, m))
    hydro_matrix = np.zeros((n, m))
    polymers = []
    reference_matrix = np.zeros((n, m))

    polymer_positions = []
    for i in range(poly_length):
        if i == 0:
            x = np.random.randint(0, n)
            y = np.random.randint(0, m)
            while matrix[x][y] != 0:
                x = np.random.randint(0, n)
                y = np.random.randint(0, m)
            matrix[x][y] = charges[i]
            reference_matrix[x][y] = i + 1
            hydro_matrix[x][y] = hydro[i]

            polymer_positions.append((x, y))
        else:
            x, y = polymer_positions[-1]
            trapped = True
            direction_lst = [[1, 0], [0, 1], [1, 1], [0, -1], [-1, 0], [-1, -1], [1, -1], [-1, 1]]
            for p in range(8):
                nx = int(x + direction_lst[p][0])
                ny = int(y + direction_lst[p][1])

                if 0 <= nx < n and 0 <= ny < m and matrix[nx][ny] == 0 and hydro_matrix[nx][ny] == 0:
                    trapped = False
            if trapped == True:
                return [], [], []
                # return []

            dx, dy = np.random.randint(-1, 2, 2)
            x_new, y_new = x + dx, y + dy
            while x_new < 0 or x_new >= n or y_new < 0 or y_new >= m or matrix[x_new][y_new] != 0:
                dx, dy = np.random.randint(-1, 2, 2)
                x_new, y_new = x + dx, y + dy
            matrix[x_new][y_new] = charges[i]
            reference_matrix[x_new][y_new] = i + 1
            hydro_matrix[x_new][y_new] = hydro[i]

            polymer_positions.append((x_new, y_new))

    return matrix, hydro_matrix, reference_matrix
    # return [matrix, hydro_matrix, reference_matrix]
# , polymers,hydro_matrix

'''task number 5'''


def initial_energy(matrix):
    up = np.roll(matrix, -1, axis=0)  # shift up
    up[-1] = 0  # set last row to zero

    bottom = np.roll(matrix, 1, axis=0)
    bottom[0] = 0

    left = np.roll(matrix, 0, axis=-1)
    left[:, 0] = 0

    right = np.roll(matrix, 0, axis=1)
    left[:, -1] = 0

    return np.exp(-np.sum(matrix * up)) + np.exp(-np.sum(matrix * bottom)) + np.exp(-np.sum(matrix * left)) + np.exp(
        -np.sum(matrix * right))


def initial_hydro(hydro_matrix):
    # this def takes in the matrix of hydrophobicity, the number of cols/rows, and an negative arbitrary number
    # returns the energy
    rows = len(hydro_matrix)
    cols = len(hydro_matrix[0])

    direction = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
    # direction vector used to visit the neighbors
    energy_sum = 0
    for i in range(cols):
        for j in range(rows):
            for k in range(4):
                # visit the neighbors
                nx = i + direction[k][0]
                ny = j + direction[k][1]
                if 0 <= nx < cols and 0 <= ny < rows:
                    if hydro_matrix[i][j] == 0:
                        energy_sum += 0.01 * hydro_matrix[ny][nx]
                        if hydro_matrix[ny][nx] < 0:
                            energy_sum -= 0.02 * hydro_matrix[ny][nx]
                    if hydro_matrix[i][j] < 0 and hydro_matrix[ny][nx] < 0:
                        energy_sum += hydro_matrix[i][j] * hydro_matrix[ny][nx]

                    # if matrix = 0 and neighbor >0:
                    # energy is positive  (neightbor * 0.01)
                    # if matrix = 0 and neighbor <0:
                    # energy is negative (-neightbor * 0.01)/
                    # if matrix <0 and neighbor <0:
                    # enegy is positive  (neightbor * matrix)
                    # if matrxi > 0 and neightbor is >0:
                    # energy is zero
                    # if matrix < 0 and neightbor is >0:

    return -1 * energy_sum

def calculate_energy(matrix):
    # Shift the matrix in four directions (up, down, left, right)
    up_shifted = np.roll(matrix, -1, axis=0)
    down_shifted = np.roll(matrix, 1, axis=0)
    left_shifted = np.roll(matrix, -1, axis=1)
    right_shifted = np.roll(matrix, 1, axis=1)

    # Calculate the energy for each case
    positive_neighbor = np.where(up_shifted > 0, up_shifted * 0.01, 0)
    negative_neighbor = np.where(down_shifted < 0, down_shifted * (-0.01), 0)
    both_negative = np.where(np.logical_and(left_shifted < 0, np.less(matrix, 0)), matrix * left_shifted, 0)

    # Calculate the total energy
    total_energy = positive_neighbor + negative_neighbor + both_negative

    return -1*total_energy.sum()
#w1 #w2
def update_energy(w1,w2,charge_matrix,hydro_matrix):
    # assign weights at 1 for now
    return w1*initial_energy(charge_matrix)+w2*calculate_energy(hydro_matrix)
#temperature
def accept_move(new_energy, energy, temperature):
    k = 1.380694*10**-23
    if new_energy<energy:
        return True
    else:
        """
        probability = math.e**-((ew_energy-energy)/
        random_integer = random.randint(0, 1)
        return random_integer<=probability
        """

        newprob_oldprob = math.e**-((new_energy-energy)/(k*temperature))
        probability = newprob_oldprob/(newprob_oldprob+1)
        random_integer = random.uniform(0, 1)
        return random_integer<=probability
        #get physical energy
        #normalize energy from 0,1
        #sigmoid function
"""IS THE ENERGY PROPERLY SET UP?"""


# temperature, iteration
def metropolis(matrix_charge, matrix_hydro, ref_matrix, energy, temperature, energy_lst, it):
    row = len(matrix_charge)
    col = len(matrix_charge[0])
    charge_copy = matrix_charge.copy()
    hydro_copy = matrix_hydro.copy()

    x, y = np.random.randint(0, row - 1), np.random.randint(0, col - 1)
    while (charge_copy[x][y] == 0):
        x, y = np.random.randint(0, row - 1), np.random.randint(0, col - 1)
    directions = random.choice([(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, -1), (-1, 1), (1, 1)])
    new_x, new_y = x + directions[0], y + directions[1]
    if 0 <= new_x < row and 0 <= new_y < col and charge_copy[new_x][new_y] == 0:
        original_charge = charge_copy[x][y]
        original_hydro = hydro_copy[x][y]
        charge_copy[x][y] = charge_copy[new_x][new_y]
        hydro_copy[x][y] = hydro_copy[new_x][new_y]
        charge_copy[new_x][new_y] = original_charge
        hydro_copy[new_x][new_y] = original_hydro

        new_energy = update_energy(1, 1, charge_copy, hydro_copy)
        # check connectivity

        connectivity = False
        direction_lst = [[1, 0], [0, 1], [1, 1], [0, -1], [-1, 0], [-1, -1], [1, -1], [-1, 1]]
        cnt = 0
        for i in range(8):
            nx = new_x + direction_lst[i][0]
            ny = new_y + direction_lst[i][1]
            if 0 <= nx < row and 0 <= ny < col and (
                    ref_matrix[nx][ny] == ref_matrix[x][y] + 1 or ref_matrix[nx][ny] == ref_matrix[x][y] - 1):
                cnt += 1
        if cnt == 2:
            connectivity = True

        if accept_move(new_energy, energy, temperature) and connectivity == True:
            matrix_charge = charge_copy
            matrix_hydro = hydro_copy
            energy = new_energy
            energy_lst.append(energy)
            it += 1
            # change the ref_matrix

            ref_temp = ref_matrix[new_x][new_y]
            ref_matrix[new_x][new_y] = ref_matrix[x][y]
            ref_matrix[x][y] = ref_temp

    return matrix_charge, matrix_hydro, ref_matrix, energy_lst, energy, it
    # return [matrix_charge, matrix_hydro,ref_matrix, energy_lst, energy, it]



def plot_graph(data):
    x = list(range(len(data[0])))  # Create x-axis values
    plt.plot(x, data[0], color='blue', label='Line 1')
    plt.plot(x, data[1], color='red', label='Line 2')
    plt.plot(x, data[2], color='yellow', label='Line 3')# Plot the data
    plt.plot(x, data[3], color='green', label='Line 4')
    plt.plot(x, data[4], color='black', label='Line 5')
    plt.xlabel("Number of Simulations")
    plt.ylabel("Energy")
    plt.title("Energy through iteration")
    plt.savefig('foo.png',dpi = 400)
    return plt.show()


def test2(ff2,p_pdb_file,p_ph,p_w1,p_w2,p_iteration, p_temp):

    df = df_2_order(pdb_2_df(p_pdb_file))
    protein_seq = df_to_sequence(df)

    charge_sequence = seq_to_charge(protein_seq, p_ph)
    # print(charge_sequence)
    dis_sequence = seq_disorder(protein_seq)
    hydro_sequence = seq_hydrophobicity(protein_seq)
    poly_length_lst = len(protein_seq)
    print("start")
    charge_matrix, hydro_matrix, ref_matrix = initialize_polymer(int(len(protein_seq) / 4), int(len(protein_seq) / 4),
                                                             charge_sequence, hydro_sequence, poly_length_lst)
    cnt = 0
    while len(charge_matrix) == 0:
        cnt+=1
        charge_matrix, hydro_matrix, ref_matrix = initialize_polymer(int(len(protein_seq) / 4),
                                                                     int(len(protein_seq) / 4 ), charge_sequence,
                                                                     hydro_sequence, poly_length_lst)
        #print(cnt)
    # print(charge_matrix)
    # print(hydro_matrix)
    # print(ref_matrix)
    #print("done")
    #charge_hm = sns.heatmap(np.array(charge_matrix),cbar=False)
    ref_hm = sns.heatmap(np.array(ref_matrix),cbar=False)
    #hydro_hm = sns.heatmap(np.array(hydro_matrix),cbar=False)
    #figure1 = charge_hm.get_figure()
    #figure1.savefig(ff1+'.png', dpi=400)
    #ff1+='i'
    #plt.show()
    figure2 = ref_hm.get_figure()
    figure2.savefig(ff2+'.png', dpi=400)
    ff2+='i'
    plt.show()
    #figure3 = hydro_hm.get_figure()
    #figure3.savefig(ff3+'.png', dpi=400)
    #ff3+='i'
    #plt.show()
    energy = update_energy(p_w1, p_w2, charge_matrix, hydro_matrix)
    energy_lst = []
    energy_lst.append(energy)
    it = 0
    charge_alignment = []
    hydro_alighment = []
    min_energy = 99999
    while it < p_iteration:
        charge_matrix, hydro_matrix, ref_matrix, energy_lst, energy, it = metropolis(charge_matrix, hydro_matrix,
                                                                                     ref_matrix, energy, p_temp,
                                                                                     energy_lst, it)
        if energy < min_energy:
            min_energy = energy
            charge_alignment = charge_matrix
            hydro_alignment = hydro_matrix
    # print(charge_alignment)
    # print(hydro_alignment)
    print("ok")
    return energy_lst, ref_matrix

def parse_arguments():
    """
    Parse command-line arguments using the argparse module.

    Returns:
        argparse.Namespace: object containing the parsed command-line arguments
    """

    '''task number 6'''
    # think about the parser variable you want to input later
    # filename for pdb
    # output directory
    # ph
    parser = argparse.ArgumentParser(description='cgRNA STS design scanner from introns')

    parser.add_argument('-type', '--input_type', type=str, default="PDB", help='input type FASTA or PDB')
    parser.add_argument('-fasta', '--fasta_file', type=str, help='FASTA sequence')
    parser.add_argument('-pdb', '--pdb_file', type=str, help='path to the input pdb file')

    parser.add_argument('-rt', '--trial', type=int, default=5, help='number of different random versions of the protein generated(default 5)')
    parser.add_argument('-d', '--directory', type=str, default='~', help='output_directory')
    parser.add_argument('-t', '--temperature', type=float, default=0, help='Temperature in K')
    parser.add_argument('-i', '--iteration', type=int, default=100, help='Simulation iterations defaulted 100 (recommand 500-5000)')
    parser.add_argument('-w1', '--weight1', type=float, default=1, help='weight of charge used in energy calculation')
    parser.add_argument('-w2', '--weight2', type=float, default=1,
                        help='weight of hydrophobicity used in energy calculation')
    parser.add_argument('-ph', '--ph', type=float, default=7.4, help='ph level of the simultation environment')

    parser.add_argument('-n', '--number_protein', type=int, default=1, help='number of protein in the simulation')

    parser.add_argument('-lower', '--lower', type=int, default=16, help='minimum window stride size for visualization')
    parser.add_argument('-upper', '--upper', type=int, default=103, help='maximum window stride size for visualization')
    parser.add_argument('-range', type=str,
                        help='The range of the number of sequences in FASTA file to scan (e.g., "0-100")')
    parser.add_argument('-duplex', '--duplex_length', type=int, default=103, help='length of the nucleotides in duplex')
    parser.add_argument('-u', '--unstructure_length', type=int, default=0,
                        help='length of the unstructured nucleotides')
    return parser.parse_args()

args = parse_arguments()


p_input_type = args.input_type
print(p_input_type)
p_pdb_file = args.pdb_file
print(p_pdb_file)
p_temperature = args.temperature
print(p_temperature)
p_weight1 = args.weight1
print(p_weight1)
p_weight2 = args.weight2
print(p_weight2)
p_ph = args.ph
print(p_ph)
p_iteration = args.iteration
print(p_iteration)
p_trial = args.trial
print(p_trial)
p_dir_out = args.directory
print(p_dir_out)

energy_matrix = []
ref_mat = []
file_name_fig = "cf"
f1 = 'ci'
f2 = 'ri'
f3 = 'hi'
for i in range(5):

    a , ref = test2(f2,p_pdb_file,p_ph,p_weight1,p_weight2,p_iteration,p_temperature)
    f1+='i'
    f2+='i'
    f3+='i'
    energy_matrix.append(a)
    ref_mat.append(ref)
    ref_hm = sns.heatmap(np.array(ref),cbar=False)
    figure4 = ref_hm.get_figure()
    file_name_fig+='i'
    figure4.savefig(file_name_fig + '.png', dpi=400)
    plt.show()

plot_graph(energy_matrix)
print(energy_matrix)
ref = ref.tolist()
json_string = json.dumps(ref, indent=None)
    # Save the JSON string to the text file

with open('ref_matrix.txt', 'w') as file:
    file.write(json_string)




# Convert the list to a JSON-formatted string

#print(energy_matrix)
#print(ref)


print("p1 done")




