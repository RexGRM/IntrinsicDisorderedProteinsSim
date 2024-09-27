import argparse


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

    parser.add_argument('-type', '--input_type', type=str, default="FASTA", help='input type FASTA or PDB')
    parser.add_argument('-fasta', '--fasta_file', type=str, help='FASTA sequence')
    parser.add_argument('-pdb', '--pdb_file', type=str, help='path to the input pdb file')



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
temp = args.pdb1_file
ph_level = args.pdb2_file
print(temp)
print(ph_level)



