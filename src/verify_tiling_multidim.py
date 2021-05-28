import json
import sys


ROTATIONS = []


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


def is_same_block(block1, block2):
    block1n = norm0(block1)
    block2n = norm0(block2)
    for rot in ROTATIONS:
        block2n = norm0(rotate(rot, block2))
        if set(block1n) == set(block2n):
            return True
    return False


def get_determinant(matrix, used_cols):
    dim = len(matrix)
    n = len(used_cols)
    if n == dim:
        return 1
    det = 0
    sign = -1
    for i in range(dim):
        if i in used_cols:
            continue
        sign = -sign
        used_cols2 = used_cols + [i]
        det += sign * matrix[n][i] * get_determinant(matrix, used_cols2)
    return det


def multiply(matrix, vector):
    product = [0] * len(matrix)
    for i in range(len(matrix)):
        for j in range(len(vector)):
            product[i] += matrix[i][j] * vector[j]
    return tuple(product)


def main():
    data = json.loads(sys.stdin.read())
    dim = len(data["original_block"][0])
    init_rotations(dim)
    for b in data["base_blocks"]:
        if not is_same_block(b, data["original_block"]):
            print(f"Block does not match original: { b }")
            sys.exit(1)

    offsets = data["offsets"]

    determinant = get_determinant(offsets, [])
    base_volume = len(data["original_block"]) * len(data["base_blocks"])
    if abs(determinant) != base_volume:
        # this can not fit
        print(f"determinant:{determinant}, expected:{base_volume}")
        sys.exit(1)

    # For each cell in base_blocks, calculate its coordinates in the base
    # defined by offsets, these will be some integer / determinant.
    # If we get only the numerator modulo determinant, these should all
    # be unique triplets.

    # inverse of transposed of offsets matrix, multiplied by determinant
    inverse = [[None] * dim for i in range(dim)]
    for i in range(dim):
        for j in range(dim):
            sign = 1 if (i + j) % 2 == 0 else -1
            matrix = [[None] * (dim - 1) for i in range(dim - 1)]
            for x in range(dim - 1):
                for y in range(dim - 1):
                    x1 = x if x < i else x + 1
                    y1 = y if y < j else y + 1
                    matrix[x][y] = offsets[x1][y1]
            inverse[i][j] = sign * get_determinant(matrix, [])

    occupied = {}
    for b in data["base_blocks"]:
        for cube in b:
            cube1 = multiply(inverse, cube)
            if cube1 in occupied:
                print(f"collision between {cube} and ({occupied[cube1]})")
                sys.exit(1)
            occupied[cube1] = cube

    print("OK")
    sys.exit(0)

if __name__ == "__main__":
    main()
