import json
import sys
import itertools
import time


ROTATIONS = []

times = [0, 0]


def next_grid(dim, volume):
    # Generates all possible matrices representing a volume <volume>
    # defining distinct lattices in <dim>-dimensional space.
    # See: https://oeis.org/A001001 for dim=3:
    # These sublattices are in 1-1 correspondence with matrices
    # [a b d]
    # [0 c e]
    # [0 0 f]
    # with acf = n, b = 0..c-1, d = 0..f-1, e = 0..f-1.

    # Here we generalize for higher dimensions


    # next_diagonal() enumerates all possible diagonals by recursive
    # factorization of volume

    diagonal = [0] * dim
    def next_diagonal(n, remaining_volume):
        if n == dim - 1:
            diagonal[n] = remaining_volume
            # we need to make a copy here!
            yield list(diagonal)
            return
        i = 1
        while i * i <= remaining_volume:
            if remaining_volume % i != 0:
                i += 1
                continue
            diagonal[n] = i
            yield from next_diagonal(n+1, remaining_volume // i)
            if i * i != remaining_volume:
                diagonal[n] = remaining_volume // i
                yield from next_diagonal(n+1, i)
            i += 1

    offsets = [None] * dim
    def next_offsets(diag, n, k):
        if n == dim:
            offsets_tuple = tuple([tuple(x) for x in offsets])
            yield offsets_tuple
            return
        if k == dim:
            yield from next_offsets(diag, n + 1, n + 1)
            return
        if k == n:
            offsets[n] = [0] * dim
            offsets[n][n] = diag[n]
            yield from next_offsets(diag, n, k + 1)
            return
        for i in range(diag[k]):
            offsets[n][k] = i
            yield from next_offsets(diag, n, k + 1)

    diagonals = list(next_diagonal(0, volume))
    # Not necessary, but we sort here to check more "round" matrices first
    # and leave the "flat" ones for later. Makes visualisation in 3d better.
    diagonals = sorted(diagonals, key=max)
    for diag in diagonals:
        for m in next_offsets(diag, 0, 0):
            yield m


def norm0(block):
    mind = list(block[0])
    for b in block:
        for i in range(len(b)):
            mind[i] = min(mind[i], b[i])
    nb = []
    for b in block:
        nb.append(tuple([b[i] - mind[i] for i in range(len(b))]))
    return tuple(nb)


def rotate(rot, block):
    rotated = []
    for b in block:
        rb = []
        for i in range(len(b)):
            r = 0
            for j in range(len(b)):
                r += rot[i][j] * b[j]
            rb.append(r)
        rotated.append(tuple(rb))
    return rotated


def rotate_inv(rot, block):
    rotated = []
    for b in block:
        rb = []
        for i in range(len(b)):
            r = 0
            for j in range(len(b)):
                r += rot[j][i] * b[j]
            rb.append(r)
        rotated.append(tuple(rb))
    return rotated


def init_rotations(dim):
    def next_rotation(rot, used, n, parity):
        if n == len(rot):
            yield rot, parity % 2
            return
        p2 = parity
        for i in range(len(rot)):
            if used[i]:
                continue
            rot[n][i] = 1
            used[i] = True
            yield from next_rotation(rot, used, n + 1, p2)
            rot[n][i] = -1
            yield from next_rotation(rot, used, n + 1, p2 + 1)
            rot[n][i] = 0
            used[i] = False
            p2 += 1

    rot0 = [[0] * dim for i in range(dim)]
    used0 = [False] * dim
    for rot, p in next_rotation(rot0, used0, 0, 0):
        r = tuple([tuple(rot[i]) for i in range(dim)])
        if p == 0:
            ROTATIONS.append(r)


def norm(block):
    dim = len(block[0])
    block = norm0(block)
    best = tuple(sorted(block))
    rot0 = [[0] * dim for i in range(dim)]
    used0 = [False] * dim
    for rot in ROTATIONS:
        block = norm0(rotate(rot, block))
        sblock = tuple(sorted(block))
        if sblock < best:
            best = sblock
        block = norm0(rotate_inv(rot, block))

    return tuple(best)


def get_cover(cache, block, pos, idx, offsets):
    # block[idx] shall be positioned at pos
    cells_covered = set()
    diag = [offsets[i][i] for i in range(len(offsets))]
    delta = [pos[i] - block[idx][i] for i in range(len(block[idx]))]

    for b in block:
        p1 = [b[i] + delta[i] for i in range(len(b))]
        p2 = tuple(p1)
        if p2 in cache:
            cells_covered.add(cache[p2])
            continue

        for i in range(len(b)):
            di = p1[i] // offsets[i][i]
            for j in range(i, len(b)):
                p1[j] -= di * offsets[i][j]
        s = 0
        for i in range(len(b)):
            s *= diag[i]
            s += p1[i]
        cells_covered.add(s)
        cache[p2] = s
    return tuple(cells_covered)


def place(block, pos, idx, offsets):
    dim = len(offsets)
    delta = [pos[i] - block[idx][i] for i in range(dim)]

    translated = []
    for b in block:
        translated.append(tuple([b[i] + delta[i] for i in range(dim)]))
    return tuple(translated)


def tile_for_offsets(block, offsets):
    # offsets is expected here to by a triangular matrix
    # as generated by next_grid()

    dim = len(offsets)

    volume = 1
    for i in range(dim):
        volume *= offsets[i][i]

    # Below is the description for the 3d case, here we generalize for any dimension

    # Placing a block at some position will cover some cells in the basic
    # parallelepiped [0..D0-1][0..D1-1][0..D2-1].
    # We enumerate each x, y, z of these as x * D1 * D2 + y * D2 + z.
    # possibilities[i] will be a dictionary with keys being all
    # possible such tuples that include i (values in the dictionaries are data to
    # help reconstruct the solution at the end).
    # Then, if we try to fill cell i, we can check one by one the keys
    # in possibilities[i].
    #
    # Example:
    # If a block in some orientation and translation covers cells: (2, 5, 13, 24)
    # (cells are defined by an integer between 0..volume-1 as described above),
    # Then possibilities[2], possibilities[5], possibilities[13] and possibilities[24]
    # will all contain the key (2, 5, 13, 24) each.
    # possibilities[2][(2, 5, 13, 24)] will be some data structure helping later
    # for generating the actual result data.
    #
    # UPDATE:
    # Because of the way the dfs is run, we don't need to store the set at all
    # positions, but only in the smallest position. In the above example, when
    # we are trying to fill cell 5, there is no need to check set (2, 5, 13, 24),
    # because we know cell 2 is already covered!

    possibilities = [{} for i in range(volume)]

    # create all rotations
    def next_rotation():
        for rot in ROTATIONS:
            b = block
            b = norm0(rotate(rot, b))
            yield b

    point = [0] * dim
    def next_point(n):
        if n == dim:
            yield tuple(point)
            return
        for i in range(offsets[n][n]):
            point[n] = i
            yield from next_point(n + 1)

    t = time.time()
    all_rotations = list(next_rotation())
    all_points = list(next_point(0))
    cache = {}
    for b in all_rotations:
        for p in all_points:
            for j in range(len(block)):
                bitmap = get_cover(cache, b, p, j, offsets)
                if len(bitmap) != len(block):
                    continue
                found = True
                s = min(bitmap)
                possibilities[s][bitmap] = (b, p, j)
    times[0] += time.time() - t

    covered = set()
    solution = [None] * volume

    def dfs(idx):
        if idx == volume:
            return True
        if idx in covered:
            return dfs(idx + 1)
        for bitmap, data in possibilities[idx].items():
            collision = False
            for cell in bitmap:
                if cell in covered:
                    collision = True
                    break
            if collision:
                continue
            for cell in bitmap:
                covered.add(cell)
            solution[idx] = data
            if dfs(idx + 1):
                return True
            for cell in bitmap:
                covered.remove(cell)
            solution[idx] = None

            if idx == 0:
                # Since we try all possible permutations of the offsets matrix,
                # we can start the dfs by positioning a random block in
                # a random position to avoid checking symmetrically equivalent
                # cases many times.
                # This is easy done here by breaking the iteration after the first
                # loop in the first step of recursion (idx == 0)
                break
        return False

    if len(possibilities[0]) == 0:
        # Nothing fits in the first cell...
        return None

    t = time.time()
    if dfs(0):
        times[1] += time.time() - t
        result = {
                "original_block": block,
                "base_blocks": [],
                "offsets": offsets
        }
        for s in solution:
            if s is None:
                continue
            b, p, i = s
            result["base_blocks"].append(place(b, p, i, offsets))
        return result
    times[1] += time.time() - t
    return None


def tile_space(block, N0, k0):
    """
    :param block: tuple of coordinates-tuples defining the block to use for tiling
                  e.g. ((0,0,0), (1,0,0), (1,1,0), (2,1,0), (2,2,0)) would be the
                  W-Pentomino laying flat on the z=0 plane.
                  The block can be anywhere in space, does not need to be at
                  the origin.
    :return: a dictionary with keys "base_blocks": a list of non overlapping blocks
                  in the same format as the input, and "offsets": a list of 3 vectors
                  indicating the directions to repeat the base_blocks to tile the
                  whole space. That means, there is a block at every position:
                  base_blocks[i] + a * offsets[0] + b * offsets[1] + c * offsets[2]
                  (vector operations).
    """
    N = len(block)
    dim = len(block[0])
    init_rotations(dim)
    tiling = None
    num_blocks_in_pattern = 1
    while tiling is None:
        num_blocks_in_pattern += 1
        if (num_blocks_in_pattern > 2):
            return
        sys.stderr.write(f"N:{num_blocks_in_pattern}\n")
        k = 0
        for offsets in next_grid(dim, num_blocks_in_pattern * N):
            k += 1
            if (num_blocks_in_pattern, k) < (N0, k0):
                continue
            if k % 1000 == 0:
                sys.stderr.write(f"k:{k}\n")
            tiling = tile_for_offsets(block, offsets)
            if tiling is not None:
                break

    return tiling


if __name__ == "__main__":
    data = json.loads(sys.stdin.read().strip())
    if len(sys.argv) > 1:
        N0 = int(sys.argv[1])
        k0 = int(sys.argv[2])
    else:
        N0 = 0
        k0 = 0
    tiling = tile_space(data, N0, k0)
    sys.stderr.write(f"times: {times}\n")
    print(json.dumps(tiling))
