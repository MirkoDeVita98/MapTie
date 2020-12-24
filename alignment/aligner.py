import constants
from constants import score, MATCH, MISMATCH, GAP
import numpy as np
from numba import njit
import time





@njit
def build_cigar(
        cigar: str, op: str, num: int, new_op: str, matches: int = 0) -> tuple:
    if op == '':
        return cigar, num + 1, new_op
    else:
        if op == "=" and matches != 0:
            return str(num + matches) + op + cigar, 1, new_op
        elif op != "=" and matches != 0:
            return str(matches) + "=" + str(num) + op + cigar, 1, new_op
        else:
            return str(num) + op + cigar, 1, new_op


# @njit
def cigarFixing(cigar: str, new_cigar: str) -> str:

    if len(cigar) > 0 and len(new_cigar) > 0 and cigar[-1] == "=":
        i = len(cigar) - 2
        assert i >= 0
        left = ''
        while i >= 0 and cigar[i].isdigit():
            left = cigar[i] + left
            i -= 1

        j = 0
        right = ''
        while j < len(new_cigar) and new_cigar[j].isdigit():
            right += new_cigar[j]
            j += 1

        if new_cigar[j] != "=":
            return cigar + new_cigar

        return cigar[0: i + 1] + str(int(left) + int(right)) + new_cigar[j: len(new_cigar)]
    else:
        return cigar + new_cigar


@njit
def build_alignment(reference: str, read: str,
                    dp: np.ndarray, n: int, m: int, k: int,
                    matches: int, last_region: str) -> str:
    cigar = ''
    op = ''
    num = 0
    # start_time = time.time()
    while n > 0 and m > 0:
        num += 1
        if inside_band(n - 1, m, k):
            if inside_band(n, m - 1, k):
                if dp[n, m] == dp[n - 1, m - 1] + s(reference[n - 1], read[m - 1]):
                    direction = 0
                elif dp[n, m] == dp[n - 1, m] + score(GAP):
                    direction = 1
                else:
                    direction = 2

            else:
                if dp[n, m] == dp[n - 1, m - 1] + s(reference[n - 1], read[m - 1]):
                    direction = 0
                else:
                    direction = 1

        else:
            if inside_band(n, m - 1, k):
                if dp[n, m] == dp[n - 1, m - 1] + s(reference[n - 1], read[m - 1]):
                    direction = 0
                else:
                    direction = 2
            else:
                direction = 0

        if direction == 0:
            if reference[n - 1] == read[m - 1] and op != "=":
                cigar, num, op = build_cigar(cigar, op, num - 1, "=")
            elif reference[n - 1] != read[m - 1] and op != "X":
                cigar, num, op = build_cigar(cigar, op, num - 1, "X")

            n -= 1
            m -= 1
        elif direction == 1:
            if op != "D":
                if op == '' and last_region:
                    cigar, num, op = build_cigar(cigar, op, num - 1, "S")
                elif op != "S":
                    cigar, num, op = build_cigar(cigar, op, num - 1, "D")

            n -= 1
        else:
            if op != "I":
                if op == '' and last_region:
                    cigar, num, op = build_cigar(cigar, op, num - 1, "S")
                elif op != "S":
                    cigar, num, op = build_cigar(cigar, op, num - 1, "I")
            m -= 1

    cigar, _, _ = build_cigar(cigar, op, num, '', matches)

    # print(f'Backward pass took {time.time() - start_time} seconds.')

    return cigar


@njit
def full_building(reference, read, dp, n, m, matches, last_region):
    cigar = ''
    op = ''
    num = 0
    while n > 0 and m > 0:
        num += 1

        if dp[n, m] == dp[n - 1, m - 1] + s(reference[n - 1], read[m - 1]):
            direction = 0
        elif dp[n, m] == dp[n - 1, m] + score(GAP):
            direction = 1
        else:
            direction = 2

        if direction == 0:
            if reference[n - 1] == read[m - 1] and op != "=":
                cigar, num, op = build_cigar(cigar, op, num - 1, "=")
            elif reference[n - 1] != read[m - 1] and op != "X":
                cigar, num, op = build_cigar(cigar, op, num - 1, "X")

            n -= 1
            m -= 1
        elif direction == 1:
            if op != "D":
                if op == '' and last_region:
                    cigar, num, op = build_cigar(cigar, op, num - 1, "S")
                elif op != "S":
                    cigar, num, op = build_cigar(cigar, op, num - 1, "D")

            n -= 1
        else:
            if op != "I":
                if op == '' and last_region:
                    cigar, num, op = build_cigar(cigar, op, num - 1, "S")
                elif op != "S":
                    cigar, num, op = build_cigar(cigar, op, num - 1, "I")
            m -= 1

    cigar, _, _ = build_cigar(cigar, op, num, '', matches)

    return cigar


@njit
def s(x: int, y: int) -> int:
    return score(MATCH) if x == y else score(MISMATCH)


@njit
def inside_band(i: int, j: int, k: int) -> bool:
    return abs(i - j) <= k


@njit
def banded_global_alignment(reference: str, read: str, matches: int,
                            last_region: bool, read_id: str, k: int = 6) -> tuple:
    n = len(reference)
    m = len(read)

    k = min(k, m)

    if n != m:
        # if read_id == constants.looking_for:
        #     print('Sending to global alignment from banded.')
        return global_alignment(reference, read, matches, last_region)

    # if read_id == constants.looking_for:
    #     print(f'The read   = {read}')
    #     print(f'The region = {reference}')

    # assert len(reference) == len(read)

    # start_time = time.time()
    dp = np.zeros((n + 1, m + 1))
    dp[0: k + 1, 0] = np.arange(k + 1) * score(GAP)
    dp[0, 0: k + 1] = np.arange(k + 1) * score(GAP)

    # compute the score in d[n + 1, m + 1]
    for i in range(1, n + 1):

        # jrange = np.arange(max(0, i - k), min(n + 1, i + k + 1), dtype=np.int32)

        # scores = np.array([s(reference[i - 1], read[j - 1]) for j in jrange])

        # dp[i, jrange] = dp[i - 1, jrange - 1] + scores

        # dp[i, jrange] = np.maximum(
        #     dp[i, jrange], dp[i - 1, jrange] + score(GAP))

        for j in range(max(0, i - k), min(n + 1, i + k + 1)):
            dp[i, j] = dp[i - 1, j - 1] + s(reference[i - 1], read[j - 1])
            if inside_band(i - 1, j, k):
                dp[i, j] = max(dp[i, j], dp[i - 1, j] + score(GAP))
            if inside_band(i, j - 1, k):
                dp[i, j] = max(dp[i, j], dp[i, j - 1] + score(GAP))

    # print(f'Forward pass took {time.time() - start_time} seconds.')

    # if read_id == constants.looking_for:
    #     print(f'The score = {dp[n, m]}')

    return build_alignment(reference, read, dp, n, m, k, matches, last_region), dp[n, m]


@njit
def global_alignment(reference, read, matches, last_region):
    n = len(reference)
    m = len(read)

    dp = np.zeros((n + 1, m + 1))
    dp[0: n + 1, 0] = np.arange(n + 1) * score(GAP)
    dp[0, 0: m + 1] = np.arange(m + 1) * score(GAP)

    # compute the score in d[n + 1, m + 1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = dp[i - 1, j - 1] + s(reference[i - 1], read[j - 1])
            insertion = dp[i - 1, j] + score(GAP)
            deletion = dp[i, j-1] + score(GAP)

            dp[i, j] = max(match, insertion, deletion)

    # build the alignments

    return full_building(reference, read, dp, n, m, matches, last_region), dp[n, m]

# @njit
def match_missing_regions(
        genome_string: str, read: str, matched_regions: np.ndarray,
        matched_indices: np.ndarray, seed_length: int, direction: int,
        read_id: str) -> tuple:

    if np.any(np.diff(matched_regions) <= 0):
        return '', -1

    last_region = len(read) // seed_length + (len(read) % seed_length != 0) - 1

    # if read_id == constants.looking_for:
    #     print(f'Read = {read}')
    #     print(f'Seed length = {seed_length}')
    #     print(f'Matched regions = {matched_regions}')
    #     print(f'Matched indices = {matched_indices}')
    #     print(f'Genome region = {genome_string[matched_indices[0] - (matched_regions[0] - 1) * seed_length: matched_indices[0] + len(read) - (matched_regions[0] - 1) * seed_length]}')

    previous_index = -1
    num = 0
    read_pos = 0
    cigar = ''
    full_score = 0
    for index, region in enumerate(matched_regions):
        if region == 0:
            match_position = matched_indices[index]
            num += seed_length
            read_pos += seed_length
        elif region > 0 and previous_index == -1:

            end_genome = matched_indices[index]
            match_position = start_genome = matched_indices[index] - \
                region * seed_length

            start_read = read_pos
            end_read = start_read + matched_regions[index] * seed_length

            read_pos = end_read

            new_cigar, score = banded_global_alignment(
                genome_string[start_genome: end_genome],
                read[start_read: end_read], num, False, read_id=read_id)

            full_score += score
            if score < -seed_length // 2 or full_score < -len(read) / 15:
                return '', -1
            # if read_id == constants.looking_for:
            #     print(f'Obtained score = {score}, full score = {full_score}.')

            cigar = cigarFixing(cigar, new_cigar)
            if region != last_region:
                num += seed_length
                read_pos += seed_length
            else:
                if len(read) % seed_length == 0:
                    num += seed_length
                    read_pos += seed_length
                else:
                    num += len(read) % seed_length
                    read_pos += len(read) % seed_length
        elif region == matched_regions[previous_index] + 1:

            gap = matched_indices[index] - \
                matched_indices[previous_index] - seed_length
            if gap != 0:
                if num != 0:
                    cigar = cigar + str(num) + "="
                    num = 0
                cigar = cigar + str(gap) + "D"

            if gap > 5:
                return '', -1

            if region != last_region:
                num += seed_length
                read_pos += seed_length
            else:
                if len(read) % seed_length == 0:
                    num += seed_length
                    read_pos += seed_length
                else:
                    num += len(read) % seed_length
                    read_pos += len(read) % seed_length

        elif region != matched_regions[previous_index] + 1:

            end_genome = matched_indices[index]
            start_genome = matched_indices[previous_index] + seed_length

            start_read = read_pos
            end_read = start_read + \
                (matched_regions[index] -
                 matched_regions[previous_index] - 1)*seed_length

            read_pos = end_read

            new_cigar, score = banded_global_alignment(
                genome_string[start_genome: end_genome],
                read[start_read: end_read], num, False, read_id=read_id)

            full_score += score
            if score < -seed_length // 2 or full_score < -len(read) / 15:
                return '', -1
            # if read_id == constants.looking_for:
            #     print(f'Obtained score = {score}, full score = {full_score}.')

            cigar = cigarFixing(cigar, new_cigar)

            if region != last_region:
                num = seed_length
                read_pos += seed_length
            else:

                if len(read) % seed_length == 0:
                    num = seed_length
                    read_pos += seed_length
                else:
                    num = len(read) % seed_length
                    read_pos += len(read) % seed_length
        else:
            assert False

        previous_index = index

    if read_pos < len(read):

        start_read = read_pos
        end_read = len(read)

        start_genome = matched_indices[previous_index] + seed_length
        end_genome = start_genome + (end_read - start_read)

        new_cigar, score = banded_global_alignment(
            genome_string[start_genome: end_genome],
            read[start_read: end_read], num, True, read_id=read_id)

        full_score += score
        if score < -seed_length // 2 or full_score < -len(read) / 15:
            return '', -1
        # if read_id == constants.looking_for:
        #     print(f'Obtained score = {score}, full score = {full_score}.')

        cigar = cigarFixing(cigar, new_cigar)

        num = 0

    if num != 0:
        cigar = cigarFixing(cigar, str(num) + "=")

    return cigar, match_position
