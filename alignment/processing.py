import time
import numpy as np
from numba import njit
from numba.typed import Dict
from numba.types import int32, int64
import re

from constants import REVERSE, FORWARD, complement_mapper, looking_for, \
    score, c_to_score, penalty
from alignment_output.output import TempOutput
from genome_indexing import indexing
from alignment import aligner
# from termcolor import colored


@njit
def turn(direction: int) -> int:
    return FORWARD if direction == REVERSE else REVERSE


def process_alignment_result(writer, genome: str, read_id: str,
                             read_length: int, result: tuple, reverse: bool,
                             first: bool=False, last: bool=False) -> tuple:

    if result[1] == '':
        return False, None, None

    # if read_id == looking_for:
    #     print(f'Cigar {result[1]} is for first={first} and last={last}')

    out_id = writer.append_alignment(
        read_id, result[0] + 1, result[1],
        genome[result[0]: result[0] + read_length], len(genome), read_length,
        reverse=reverse, secondary=False, first=first, last=last)

    return True, result[0], out_id  # the read has been processed


@njit
def rec(locations: np.ndarray, cluster_ixs: np.ndarray,
        valid_regions: np.ndarray, ii: int, tmp: list,
        seed_length: int) -> list:

    if ii == len(valid_regions):
        return tmp

    diff = valid_regions[ii] - valid_regions[ii - 1]

    new_tmp = []
    for sset in tmp:
        last_region_start = cluster_ixs[sset[-1]]
        for next_region_ix in np.where(locations[valid_regions[ii]])[0]:
            if seed_length * diff <= cluster_ixs[next_region_ix] - last_region_start <= (seed_length + 10) * diff:
                new_tmp.append(sset[:] + [next_region_ix])
            # else:  # This makes it very slow, but catches everything
            #     new_tmp.append([next_region_ix])

    return rec(
        locations, cluster_ixs, valid_regions, ii + 1, new_tmp, seed_length)


@njit
def rec_init(
        locations: np.ndarray, cluster_ixs: np.ndarray,
        valid_regions: np.ndarray, seed_length: int):

    new_tmp = [[location]
               for location in np.where(locations[valid_regions[0]])[0]]

    return rec(locations, cluster_ixs, valid_regions, 1, new_tmp, seed_length)


@njit
def get_candidates(
        ixs: np.ndarray, matched_regions: np.ndarray, seed_length: int,
        read_id: str) -> list:

    if len(ixs) <= 1:  # TODO
        return None

    args = np.argsort(ixs)
    ixs = ixs[args]
    matched_regions = matched_regions[args]
    diffs = np.diff(ixs)
    borders = np.where(diffs > 3 * seed_length)[0] + 1
    all_borders = np.hstack((np.asarray([0]), borders, np.asarray([ixs.size])))

    precandidates = []

    for ix in range(len(all_borders) - 1):

        cluster = np.arange(all_borders[ix], all_borders[ix + 1])

        if len(cluster) <= 1:
            continue

        cluster_ixs = ixs[cluster]
        cluster_matched_regions = matched_regions[cluster]

        unique_matched_regions = np.unique(cluster_matched_regions)

        locations = np.zeros(
            shape=(unique_matched_regions[-1] + 1, cluster_matched_regions.size))
        for region_no in unique_matched_regions:
            # locations[region_no] = np.where(
            #     cluster_matched_regions == region_no)[0]
            locations[region_no, :] = cluster_matched_regions == region_no

        tmps = rec_init(locations, cluster_ixs,
                        unique_matched_regions, seed_length)

        if len(tmps) > 0:
            precandidates.extend(
                [(cluster_ixs[np.asarray(tmp)],
                  cluster_matched_regions[np.asarray(tmp)])
                 for tmp in tmps])

    precandidates.sort(key=lambda el: -el[0].size)

    return precandidates


@njit
def process_ixs(
        ixs: np.ndarray, matched_regions: np.ndarray, seed_length: int,
        read_id: str) -> tuple:

    if len(ixs) <= 1:  # TODO
        return None, None, None, None

    args = np.argsort(ixs)
    ixs = ixs[args]
    matched_regions = matched_regions[args]
    diffs = np.diff(ixs)
    borders = np.where(diffs > 3 * seed_length)[0] + 1
    all_borders = np.hstack((np.asarray([0]), borders, np.asarray([ixs.size])))

    actual_ixs = []
    actual_regions = []
    tmp_cluster_sizes = []
    tmp_all_borders = [0]
    tmp_last_border = 0
    for ix in range(len(all_borders) - 1):

        cluster = np.arange(all_borders[ix], all_borders[ix + 1])

        if len(cluster) <= 1:
            continue

        cluster_ixs = ixs[cluster]
        cluster_matched_regions = matched_regions[cluster]

        unique_matched_regions = np.unique(cluster_matched_regions)

        min_region_in_cluster = unique_matched_regions[0]

        cluster_filter = cluster_matched_regions == min_region_in_cluster
        cluster_filter[np.argmax(
            cluster_matched_regions == min_region_in_cluster)] = False
        cluster_filter = np.where(cluster_filter)[0]
        cluster_matched_regions = np.delete(
            cluster_matched_regions, cluster_filter)
        cluster_ixs = np.delete(cluster_ixs, cluster_filter)

        for region_no in cluster_matched_regions[
                cluster_matched_regions > min_region_in_cluster]:

            target = unique_matched_regions[unique_matched_regions < region_no]
            if len(target) == 0:
                continue

            smaller_region = np.max(target)

            mask = cluster_matched_regions == smaller_region
            if np.sum(mask) == 0:
                continue

            first_smaller = np.where(mask)[0][0]

            r = np.arange(len(cluster_matched_regions))
            region_in_cluster = cluster_matched_regions == region_no

            regions_filter_template = np.asarray(region_in_cluster)

            regions_filter_1 = regions_filter_template * \
                np.asarray(r < first_smaller)

            regions_filter_2 = regions_filter_template * np.asarray(
                cluster_ixs < cluster_ixs[first_smaller] + seed_length)

            regions_filter = regions_filter_1 + regions_filter_2

            regions_filter = np.where(regions_filter)[0]
            cluster_matched_regions = np.delete(
                cluster_matched_regions, regions_filter)
            cluster_ixs = np.delete(cluster_ixs, regions_filter)

            regions_filter = cluster_matched_regions == region_no
            regions_filter[np.argmax(
                cluster_matched_regions == region_no)] = False
            regions_filter = np.where(regions_filter)[0]
            cluster_matched_regions = np.delete(
                cluster_matched_regions, regions_filter)
            cluster_ixs = np.delete(cluster_ixs, regions_filter)

        actual_ixs.extend(cluster_ixs)
        actual_regions.extend(cluster_matched_regions)
        tmp_cluster_sizes.append(len(cluster_matched_regions))
        tmp_last_border += len(cluster_matched_regions)
        tmp_all_borders.append(tmp_last_border)

    if len(actual_ixs) <= 1:  # TODO
        return None, None, None, None

    return (np.asarray(actual_ixs),
            np.asarray(actual_regions),
            np.asarray(tmp_all_borders),
            np.asarray(tmp_cluster_sizes))


@njit
def process_read_all_indices(genome: str, suffix_array, read: str,
                             seed_length: int, region_boundaries: tuple,
                             read_id: str, consider_all: bool = False) -> list:

    # start_time = time.time()
    ixs, matched_regions, _ = indexing.find_all_indices(
        genome, suffix_array, read, region_boundaries, seed_length, read_id)
    # print(f'Finding all ixs took {time.time() - start_time} seconds.')

    if consider_all:

        return get_candidates(ixs, matched_regions, seed_length, read_id)

    else:
        # start_time = time.time()
        ixs, matched_regions, all_borders, cluster_sizes = process_ixs(
            ixs, matched_regions, seed_length, read_id)
        # print(f'Processing ixs took {time.time() - start_time} seconds.')

        if ixs is None:
            return None

        args = np.argsort(cluster_sizes)
        cluster_sizes = cluster_sizes[args]
        ixs_sorted = []
        while len(cluster_sizes) > 0:
            cluster_start, args = args[-1], args[:-1]
            cluster_size, cluster_sizes = cluster_sizes[-1], cluster_sizes[:-1]
            if cluster_size > 1:
                m = (ixs[all_borders[cluster_start]:
                        all_borders[cluster_start] + cluster_size],
                    matched_regions[all_borders[cluster_start]:
                                    all_borders[cluster_start] + cluster_size])
                ixs_sorted.append(m)
        return ixs_sorted


def report_time_and_status(start_time, milestone_time, successful_queries,
                           failed_queries, period, processed_queries):
    print(f'There have been {successful_queries} successful '
          f'and {failed_queries} failed queries')
    print(f'The success rate is {successful_queries / processed_queries}.')
    print(f'The failure rate is {failed_queries / processed_queries}.')
    elapsed_time = time.time() - start_time
    print(f'{processed_queries}: {elapsed_time} seconds - '
          f'{(processed_queries) / elapsed_time} queries/second')
    elapsed_milestone_time = time.time() - milestone_time
    print(f'{period}: {elapsed_milestone_time} seconds - '
          f'{(period) / elapsed_milestone_time} queries/second')


## @njit
def rc(read: str) -> str:
    # rread = read[::-1]
    # result = ''
    # for ch in rread:
    #     result += complement_mapper(ch)
    # return result
    return ''.join(complement_mapper.get(b, b) for b in reversed(read))


def match(genome: str, suffix_arrays: tuple, seed_lengths: list,
          read_id: str, read_a: str, read_b: str, consider_all: bool,
          out_writer: TempOutput) -> bool:
    both_results = []
    forward_score = -len(genome)

    seed_ix = 0
    while seed_ix < len(seed_lengths):

        seed_length = seed_lengths[seed_ix]
        suffix_array = suffix_arrays[seed_length]

        ixs_sorted_a = process_read_all_indices(
            genome, suffix_array, read_a,
            seed_length, None, read_id=read_id, consider_all=consider_all)
        
        # if read_id == looking_for:
        #     print(f'Ixs a: {ixs_sorted_a}')

        if ixs_sorted_a is not None:

            success, b_0, b_1, a_0, a_1, best_score = align_ixs_sorted(
                ixs_sorted_a, genome, suffix_array, read_a, read_b,
                seed_lengths, seed_length, consider_all=consider_all,
                direction=FORWARD, read_id=read_id)

            if success:

                if best_score == 0:
                    process_alignment_result(
                        out_writer, genome, read_id, len(read_a),
                        (a_0, a_1), FORWARD, first=True, last=False)
                    process_alignment_result(
                        out_writer, genome, read_id, len(read_b),
                        (b_0, b_1), REVERSE, first=False, last=True)
                    return True

                both_results.append(
                    ((out_writer, genome, read_id,
                        len(read_a), (a_0, a_1), FORWARD, True, False),
                        (out_writer, genome, read_id,
                            len(read_b), (b_0, b_1), REVERSE, False, True)))
                forward_score = best_score

        read_a_rc = rc(read_a)

        ixs_sorted_a = process_read_all_indices(
            genome, suffix_array, read_a_rc,
            seed_length, None, read_id=read_id, consider_all=consider_all)

        # if read_id == looking_for:
        #     print(f'Ixs a~: {ixs_sorted_a}')

        if ixs_sorted_a is not None:

            success, b_0, b_1, a_0, a_1, best_score = align_ixs_sorted(
                ixs_sorted_a, genome, suffix_array, read_a_rc, read_b,
                seed_lengths, seed_length, consider_all=consider_all,
                direction=REVERSE, read_id=read_id)

            if best_score > forward_score:
                process_alignment_result(
                    out_writer, genome, read_id, len(read_a),
                    (a_0, a_1), REVERSE, first=True, last=False)
                process_alignment_result(
                    out_writer, genome, read_id, len(read_b),
                    (b_0, b_1), FORWARD, first=False, last=True)
                return True

        if len(both_results) > 0:
            process_alignment_result(*both_results[0][0])
            process_alignment_result(*both_results[0][1])
            return True

        seed_ix += 1

    return False


def process_lines(whole_genome: dict, suffix_arrays: dict, lines, fname: str,
                  seed_lengths:list, process_no: int, period: int):

    print(f'I have to process {len(lines)} lines.')

    genome_id = whole_genome['id']
    genome = whole_genome['genome']

    out_writer = TempOutput(genome_id, process_no, fname)

    lines = list(map(lambda line: line.split(','), lines))

    start_time = time.time()
    milestone_time = start_time

    successful, failed = 0, 0
    for n, line in enumerate(lines):

        read_id, read_a, read_b = line

        matched = match(genome, suffix_arrays, seed_lengths,
                        read_id, read_a, read_b, False, out_writer)

        if matched:
            successful += 1
        else:

            corrected = match(genome, suffix_arrays, seed_lengths,
                              read_id, read_a, read_b, True, out_writer)

            if not corrected:
                print(f'Alignment for read {read_id} failed.')
                failed += 1
            else:
                successful += 1

        if n > 0 and n % period == 0:
            report_time_and_status(
                start_time, milestone_time, successful, failed, period, n + 1)
            milestone_time = time.time()

    out_writer.close()
    

## @njit
def align_ixs_sorted(ixs_sorted_a: list, genome: str, suffix_array: np.ndarray,
                     read_a: str, read_b: str, seed_lengths: list,
                     seed_length: int, direction: int, consider_all: bool,
                     read_id: str) -> tuple:

    results = []
    while ixs_sorted_a is not None and len(ixs_sorted_a) > 0:

        tmp_a = ixs_sorted_a.pop(0)

        if check_perfect_match(tmp_a[0], tmp_a[1], read_a, seed_length):

            success, result_0, result_1, best_score = check_backward_match(
                genome, suffix_array, read_b,
                seed_lengths, tmp_a[0][0], consider_all=consider_all,
                direction=turn(direction), read_id=read_id)

            if success:
                results.append((
                    result_0, result_1,
                    tmp_a[0][0], str(len(read_a)) + '=',
                    best_score))
        else:
            cigar_a, position_a = aligner.match_missing_regions(
                genome, read_a, tmp_a[1], tmp_a[0], seed_length,
                direction=direction, read_id=read_id)

            if cigar_a != '' and position_a != -1:
                success, result_0, result_1, best_score = check_backward_match(
                    genome, suffix_array, read_b,
                    seed_lengths, tmp_a[0][0], consider_all=consider_all,
                    direction=turn(direction), read_id=read_id)
                if success:
                    results.append((
                        result_0, result_1,
                        position_a, cigar_a,
                        best_score))

    if len(results) == 0:
        return False, -1, '', -1, '', -len(genome)
    else:
        return best_result(results, len(read_a))


## @njit
def check_backward_match(genome: str, suffix_array: np.ndarray, read: str,
                         seed_lengths: list, match_pos: int, direction: int,
                         consider_all: bool, read_id: str) -> dict:

    region_boundaries = (match_pos - 8 * len(read), match_pos + 8 * len(read))

    read = read if direction == FORWARD else rc(read)

    seed_ix = 0
    while seed_ix < len(seed_lengths):

        results = []

        seed_length = seed_lengths[seed_ix]

        ixs_sorted_b = process_read_all_indices(
            genome, suffix_array, read, seed_length,
            region_boundaries, read_id=read_id, consider_all=consider_all)

        while ixs_sorted_b is not None and len(ixs_sorted_b) > 0:

            tmp_b = ixs_sorted_b.pop(0)

            if check_perfect_match(tmp_b[0], tmp_b[1], read, seed_length):
                results.append((tmp_b[0][0], str(len(read)) + '='))
            else:
                cigar_b, position_b = aligner.match_missing_regions(
                    genome, read, tmp_b[1], tmp_b[0],
                    seed_length, direction=direction, read_id=read_id)

                if cigar_b != '' and position_b != -1:
                    results.append((position_b, cigar_b))

        if len(results) > 0:
            return (True, ) + best_single_result(results, len(read))

        seed_ix += 1

    return False, -1, '', -1


@njit
def check_perfect_match(ixs: np.ndarray, matched_regions: np.ndarray,
                        read: str, seed_length: int) -> bool:
    if np.max(ixs) - np.min(ixs) >= len(read) - seed_length \
            and np.max(np.abs(np.diff(ixs))) == seed_length \
                and np.all(np.diff(matched_regions) >= 0):
        return True
    return False


## @njit
def cigar_score(cigar: str, read_length: int) -> int:
    ops = [(cigar[m.end()], cigar[m.start(): m.end()])
           for m in re.finditer(r'\d+', cigar)]
    p = -read_length \
        if np.sum([int(op[1]) for op in ops]) >= 1.1 * read_length else 0
    return np.sum(
        [c_to_score[op[0]] * int(op[1]) + penalty[op[0]] for op in ops]) + p


## @njit
def best_single_result(results: list, read_length: int) -> tuple:
    scores = np.array([cigar_score(result[1], read_length) for result in results])
    best = np.argmax(scores)
    return results[best][0], results[best][1], scores[best]


## @njit
def best_result(results: list, read_length: int) -> tuple:
    scores = np.array(
        [cigar_score(result[3], read_length) + result[4] for result in results])
    best = np.argmax(scores)
    return True, results[best][0], results[best][1], results[best][2], results[best][3], scores[best]
