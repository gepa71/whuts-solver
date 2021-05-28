import itertools
import json
import sys


ROTATIONS = []
ROTATIONS_MIRRORS = []


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
        ROTATIONS_MIRRORS.append(r)
        if p == 0:
            ROTATIONS.append(r)


def norm(block):
    dim = len(block[0])
    block = norm0(block)
    best = tuple(sorted(block))
    rot0 = [[0] * dim for i in range(dim)]
    used0 = [False] * dim
    for rot in ROTATIONS_MIRRORS:
        block = norm0(rotate(rot, block))
        sblock = tuple(sorted(block))
        if sblock < best:
            best = sblock
        block = norm0(rotate_inv(rot, block))

    return tuple(best)


def generate(dim, n):
    items = [None] * (n + 1)

    items[1] = set([(tuple([0] * dim),)])

    for k in range(2, n + 1):
        items[k] = set()
        for a in items[k-1]:
            for b in a:
                for i in range(dim):
                    for delta in [-1, 1]:
                        p = list(b)
                        p[i] += delta
                        p = tuple(p)
                        if p in a:
                            continue
                        q = tuple(list(a) + [p])
                        q = norm(q)
                        items[k].add(q)
        f = open(f"polyhypercubes_{dim}_{k}.json", 'w+')
        f.write(json.dumps(sorted(items[k])))
        f.close()
        print(k, len(items[k]))
    return None


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <Dimension> <N>\n")
        sys.exit(1)
    dim = int(sys.argv[1])
    n = int(sys.argv[2])
    init_rotations(dim)
    generate(dim, n)
