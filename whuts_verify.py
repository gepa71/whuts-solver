import json
import sys


def rotate_x(block):
    return tuple([(x, -z, y) for x, y, z in block])


def rotate_y(block):
    return tuple([(z, y, -x) for x, y, z in block])


def rotate_z(block):
    return tuple([(-y, x, z) for x, y, z in block])


def norm0(block):
    minx = block[0][0]
    miny = block[0][1]
    minz = block[0][2]
    for x, y, z in block:
        minx = min(minx, x)
        miny = min(miny, y)
        minz = min(minz, z)
    return tuple([(x-minx, y-miny, z-minz) for x, y, z in block])


def is_same_block(block1, block2):
    block1n = norm0(block1)
    block2n = norm0(block2)
    for i in range(3):
        for j in range(4):
            block2n = norm0(rotate_y(block2n))
            if set(block1n) == set(block2n):
                return True
        block2n = norm0(rotate_x(block2n))
        for j in range(4):
            block2n = norm0(rotate_y(block2n))
            if set(block1n) == set(block2n):
                return True
        block2n = norm0(rotate_z(block2n))
    return False


def main():
    data = json.loads(sys.stdin.read())
    for b in data["base_blocks"]:
        if not is_same_block(b, data["original_block"]):
            print(f"Block does not match original: { b }")

    N = 50
    box = [[[None] * N for i in range(N)] for j in range(N)]
    offsets = data["offsets"]
    for i in range(-N, N):
        for j in range(-N, N):
            for k in range(-N, N):
                for b in data["base_blocks"]:
                    for x, y, z in b:
                        x1 = x + i * offsets[0][0] + j * offsets[1][0] + k * offsets[2][0]
                        y1 = y + i * offsets[0][1] + j * offsets[1][1] + k * offsets[2][1]
                        z1 = z + i * offsets[0][2] + j * offsets[1][2] + k * offsets[2][2]
                        if x1 < 0 or y1 < 0 or z1 < 0 or x1 >= N or y1 >= N or z1 >= N:
                            continue
                        if box[x1][y1][z1] is not None:
                            print(f"Duplicate at [{x1},{y1},{z1}] for i={i}, j={j}, k={k}")
                            sys.exit(1)
                        box[x1][y1][z1] = True
    for i in range(N // 2):
        for j in range(N // 2):
            for k in range(N // 2):
                if box[i][j][k] is None:
                    print(f"Empty at [{i},{j},{k}]")
                    sys.exit(1)


if __name__ == "__main__":
    main()
