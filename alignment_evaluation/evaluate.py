import os
import sys; sys.path.append('..')

from alignment_output import output
from constants import c_to_score, penalty
import re

def evaluate_output(filename):
    # read from the baseline
    with open('/Users/sat/Desktop/ETH_PROJECTS/comp_bio/team02/data_small/output_tiny_30xCov.sam') as baseline:
        baseline_arguments = []
        counter = 0
        for l in baseline:
            read_arguments = l.strip().split("\t")
            if read_arguments[0][0] != '@':
                baseline_arguments.append(read_arguments)
                counter += 1
            if counter >= 400000:
                break
        

    # read from our output
    with open(filename) as f:
        output_arguments = []
        counter = 0
        for l in f:
            read_arguments = l.strip().split("\t")
            if read_arguments[0][0] != '@':
                output_arguments.append(read_arguments)
                counter += 1
            if counter >= 400000:
                break
    
    # compare the CIGAR string (avg. edit distance as metric)
    # s_o = 0
    # s_b = 0
    # i = 0
    # for b, o in zip(baseline_arguments, output_arguments):
    #     if b[5] != o[5]:
    #         print("\n\n")
    #         print(f'Correct: id->{b[0]}, pos-> {b[3]}, cigar->{b[5]}')
    #         print(f'Ours: id->{o[0]}, pos-> {o[3]}, cigar->{o[5]}')
    tmp = list(zip(baseline_arguments, output_arguments))
    for ii in range(0, len(tmp), 2):
        c1 = tmp[ii][0][5]
        c2 = tmp[ii + 1][0][5]
        o1 = tmp[ii][1][5]
        o2 = tmp[ii + 1][1][5]

        if not(o1 == c1 and o2 == c1 or o2 == c1 and o1 == c2):
            print("\n\n")
            print(f'Correct: id->{c1} or {c2}')
            print(f'Correct: id->{o1} or {o2}')

# Levenshtein Distance Implementation

def edit_distance(string1, string2):
    matrix = [[0 for j in range(len(string1)+1)] for j in range(len(string2)+1)]
    for i in range(len(string1)+1):
        matrix[0][i] = i
        matrix[i][0] = i
    for i in range(1, len(string2)+1):
        for j in range(1, len(string1)+1):
            possible_values = [matrix[i-1][j]+1, matrix[i][j-1]+1]
            if string1[j-1] == string2[i-1]:
                possible_values.append(matrix[i-1][j-1])
            else:
                possible_values.append(matrix[i-1][j-1]+1)
            matrix[i][j] = min(possible_values)

    return matrix[len(string2)][len(string1)]


evaluate_output("/Users/sat/Desktop/ETH_PROJECTS/comp_bio/team02/results/results_output_tiny_30xCov.sam")
