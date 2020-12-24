from numba import njit
from enum import Enum

class Alphabet:
    S = -1
    A = 0
    T = 1
    C = 2
    G = 3
    N = 4
    H = 5



@njit
def score(arg):
    # if arg == MATCH:
    #     return 0
    # if arg == MISMATCH:
    #     return -1
    # if arg == GAP:
    #     return -2
    return -arg

MATCH = 0
MISMATCH = 1
GAP = 2
SCORE = [0,-1,-2]
    

c_to_score = {
    'I': score(GAP),
    'D': score(GAP),
    'X': score(MISMATCH),
    'N': score(MISMATCH),
    '=': score(MATCH),
    'S': score(MISMATCH),
}


penalty = {
    'I': -6,
    'D': -6,
    'X': 0,
    'N': 0,
    '=': 0,
    'S': 0,
}



complement_mapper = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'N': 'N',
}


numeric_mapper = {
    'A': 0,
    'C': 1,
    'G': 2,
    'N': 3,
    'T': 4,
}


# seed_lengths = [48]
seed_lengths = [48,32, 16]
# seed_lengths = [16]


UP = 0
DOWN = 1



FORWARD = 0
REVERSE = 1


looking_for = "22-11221624-"