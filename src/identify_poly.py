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


def main(filename):
    f = open(filename)
    all_polys = json.loads(f.read())
    f.close()
    data = json.loads(sys.stdin.read())
    dim = len(data[0])
    init_rotations(dim)
    for i, b in enumerate(all_polys):
        if is_same_block(b, data):
            print(i)
            return
    print("Not found")

if __name__ == "__main__":
    main(sys.argv[1])
