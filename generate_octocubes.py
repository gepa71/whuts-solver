import itertools
import json

def rotate_x(block):
    return tuple([(x, -z, y) for x, y, z in block])


def rotate_y(block):
    return tuple([(z, y, -x) for x, y, z in block])


def rotate_z(block):
    return tuple([(-y, x, z) for x, y, z in block])


def mirror(block):
    return tuple([(-x, y, z) for x, y, z in block])

def norm0(block):
    minx = block[0][0]
    miny = block[0][1]
    minz = block[0][2]
    for x, y, z in block:
        minx = min(minx, x)
        miny = min(miny, y)
        minz = min(minz, z)
    return tuple(([(x-minx, y-miny, z-minz) for x, y, z in block]))

def norm(block):
    block = norm0(block)
    best = sorted(block)
    for m in 0, 1:
        block = norm0(mirror(block))
        for i in range(3):
            for j in range(4):
                block = norm0(rotate_y(block))
                if sorted(block) < best:
                    best = sorted(block)
            block = norm0(rotate_x(block))
            for j in range(4):
                block = norm0(rotate_y(block))
                if sorted(block) < best:
                    best = sorted(block)
            block = norm0(rotate_z(block))
    return tuple(best)


def generate(n):
    items = [None] * (n + 1)

    items[1] = set([((0,0,0),)])

    for n in range(2, 9):
        items[n] = set()
        for a in items[n-1]:
            for x, y, z in a:
                for dx, dy, dz in ([-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]):
                    p = (x + dx, y + dy, z + dz)
                    if p in a:
                        continue
                    q = tuple(list(a) + [p])
                    q = norm(q)
                    items[n].add(q)
    return sorted(items[n])


if __name__ == "__main__":
    print(json.dumps(generate(8), indent=4))
