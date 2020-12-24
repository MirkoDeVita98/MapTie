
from random import randint
from random import seed

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def get_test_lines(genome_string,n=10,length=1000,r_seed=42, initial_offset=16050000):
    lines = []
    genome_length = len(genome_string['genome'])
    seed(r_seed)
    for i in range(1,n+1):
        index = randint(0, genome_length-length-initial_offset) + initial_offset
        id = "22K-" + str(i)
        read_a = genome_string['genome'][index  : index + length]
        read_b = rc(read_a)
        # read_a = list(read_a)
        # read_b = list(read_b)
        # read_a[621] = 'A'
        # read_b[261] = 'A'
        # read_a[1621] = 'A'
        # read_b[2261] = 'A'
        # read_a[6231] = 'A'
        # read_b[2361] = 'A'
        # read_a[6521] = 'A'
        # read_b[2661] = 'A'
        # read_a[6721] = 'A'
        # read_b[2631] = 'A'
        # read_a = ''.join(read_a)
        # read_b = ''.join(read_b)
        if 'N' in read_a or 'N' in read_b:
            continue
        lines.append(id + "," + read_a + "," + read_b)
    return lines


def rc(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))

