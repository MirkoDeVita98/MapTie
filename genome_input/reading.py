import os
import sys; sys.path.append('..')

import numpy as np
import pandas as pd
import ntpath
import constants
from utils.parser import load_config

def _to_numeric(s):
    s = s.translate(s.maketrans('$ATCGN#', '0123456'))
    s = s.replace('', ' ')[1: -1]
    return np.fromstring(s, dtype=np.uint8, sep=' ')


    

# TODO: pickle
def read_reference_genome(genome_file_path):
    data = load_config()

    genome_id = None
    genome = None
    genome_complement = None

    with open(data['input']['reference_abs_path']) as f:
        lines = f.readlines()
        genome_id = lines[0][1:-1]
        lines = list(map(lambda line: line[:-1], lines))
        genome = ''.join(lines[1:])

    return dict(id=genome_id, genome=genome)


def read_queries(input_file):

    ids = []
    sequences = []
    qualities = []

    with open(input_file) as f:
        for n, line in enumerate(f):
            if n % 4 == 0:  # read the identifier
                read_id, direction = line.split('/')
                ids.append(read_id[1:])
            elif n % 4 == 1:  # read the sequence
                # sequences.append(_to_numeric(line[:-1]))
                sequences.append(line[:-1])
            elif n % 4 == 3:  # read the qualities
                # qualities.append([ord(c) - 33 for c in line[:-1]])
                pass

    return dict(ids=ids, sequences=sequences, qualities=qualities)


def get_lines(file1,file2):
    data = load_config()
    combined_reads_file_name = f'{path_leaf(file1)[:-4]}_combined.csv' 
    
    combined_reads_file_path = os.path.join(
        data['output']['cache_dir'], combined_reads_file_name)

    if not os.path.isfile(combined_reads_file_path):
        preprocess_reads(file1, file2)
    
    return pd.read_csv(combined_reads_file_path)['0']

    
def preprocess_reads(file1, file2):
    forward_reads_file_name = file1
    reverse_reads_file_name = file2
    cache_dir = load_config()['output']['cache_dir']


    with open(forward_reads_file_name, 'r') as file_forward, \
        open(reverse_reads_file_name, 'r') as file_reverse:

        lines_forward = file_forward.readlines()
        lines_reverse = file_reverse.readlines()

        # print('\n'.join(lines_forward[:25]))
        # print('-----------------------------------------------------------------------------------------------------------------')

        # del(lines_forward[2::4])
        # del(lines_forward[2::4])
        # del(lines_reverse[2::4])
        # del(lines_reverse[2::4])

        # print('\n'.join(lines_forward[:25]))
        # print(lines_reverse)

        lines = [','.join(
                [lines_forward[n].split('/')[0][1:],
                 lines_forward[n + 1][:-1],
                 lines_reverse[n + 1][:-1]]
            ) for n in range(0, len(lines_forward), 4)]

        df = pd.DataFrame(lines)

        parts = forward_reads_file_name.split('.')
        df.to_csv(
            os.path.join(
                cache_dir,f'{path_leaf(file1)[:-4]}_combined.csv'))


def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)