import constants
from constants import looking_for
from genome_input import reading
from collections import defaultdict
import numpy as np
import os
import time
from utils.parser import load_config
from numba import njit
from numba.typed import Dict
from numba import types


@njit
def search(genome: str, suffix_array: np.ndarray, sequence: str) -> tuple:
    n = len(suffix_array) - 1

    l, r = 0, n
    while l < r:
        mid = (l + r) // 2
        if sequence > genome[suffix_array[mid]: suffix_array[mid] + len(sequence)]:
            l = mid + 1
        else:
            r = mid
    s = l

    r = n
    while l < r:
        mid = (l + r) // 2
        if genome[suffix_array[mid]: suffix_array[mid] + len(sequence)] == sequence:
            l = mid + 1
        else:
            r = mid

    while s < len(suffix_array) and \
            genome[suffix_array[s]: suffix_array[s] + len(sequence)] < sequence:
        s += 1

    while r >= 0 and \
            genome[suffix_array[r]: suffix_array[r] + len(sequence)] > sequence:
        r -= 1

    return s, r


@njit
def find_all_indices(genome_string: str, suffix_array: np.ndarray, 
                     sequence: str, region_boundaries: tuple,
                     seed_length: int, read_id: str) -> tuple:

    bins = Dict.empty(
        key_type=types.int64,
        value_type=types.int64,
    )
    positions = []
    ixs = []
    matched_regions = []
    slack = min(6, len(sequence) // 25)

    for ix in range(0, len(sequence), seed_length):
        
        if len(sequence) - ix < 16:
            continue

        k_min, k_max = search(
            genome_string, suffix_array, sequence[ix: ix + seed_length])

        candidates = suffix_array[k_min: k_max + 1]

        if region_boundaries is not None:
            candidates = candidates[candidates >= region_boundaries[0]]
            candidates = candidates[candidates <= region_boundaries[1]]

        n = ix // seed_length
        for candidate in candidates:
            position = candidate // len(sequence)
            if position not in positions:
                bins[position] = 0
                positions.append(position)
            if bins[position] >= n - slack:
                ixs.append(candidate)
                matched_regions.append(n)
                if position - 1 not in positions:
                    bins[position - 1] = 0
                    positions.append(position - 1)
                if position + 1 not in positions:
                    bins[position + 1] = 0
                    positions.append(position + 1)
                bins[position - 1] += 1
                bins[position] += 1
                bins[position + 1] += 1

    ixs = np.asarray(ixs)
    matched_regions = np.asarray(matched_regions)

    no_regions = len(sequence) // seed_length + 1
    mask = np.zeros(shape=(ixs.size,))
    for ii in range(len(mask)):
        mask[ii] = bins[ixs[ii] // len(sequence)] >= no_regions - slack
    mask = np.where(mask)[0]
    ixs = ixs[mask]
    matched_regions = matched_regions[mask]

    # print(f'Gathering candidates took {time.time() - start_time} seconds.')

    return ixs, matched_regions, None


def get_indexes(whole_genome, size, seed_lengths):

    genome = whole_genome['genome']
    return {
        seed_length: get_index(genome, size, seed_length)
        for seed_length in seed_lengths
    }


def get_index(genome, size, seed_length):
    data = load_config()
    suffix_array_file = f'suffix_array_index_{size}_{seed_length}_0.npy'
    suffix_array_0_path = os.path.join(data['output']['cache_dir'], suffix_array_file)
    
    return load_index(size, seed_length) \
        if os.path.isfile(suffix_array_0_path) \
        else build_suffix_array(genome, size, seed_length)


def build_suffix_array(genome, size, seed_length):

    print(f'Building the index with seed length {seed_length}.')

    suffix_array = list(range(0, len(genome) - seed_length))
    suffix_array.sort(key=lambda ix: genome[ix: ix + seed_length])

    suffix_array = np.asarray(suffix_array)

    print('Index built.')

    save_index(suffix_array, size, seed_length)

    return suffix_array


def save_index(suffix_array, size, seed_length):
    cache_dir = load_config()['output']['cache_dir']
    os.makedirs(cache_dir, exist_ok=True)

    np.save(
        os.path.join(
            cache_dir, f'suffix_array_index_{size}_{seed_length}_0'),
        suffix_array, allow_pickle=True)

    print('Index saved.')


def load_index(size, seed_length):

    print('Loading the stored index.')

    data = load_config()
    forward_suffix_array = np.load(
        os.path.join(data['output']['cache_dir'],
                     f'suffix_array_index_{size}_{seed_length}_0.npy'),
        allow_pickle=True)

    return forward_suffix_array
