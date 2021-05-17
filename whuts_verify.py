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
            sys.exit(1)

    offsets = data["offsets"]

    determinant = offsets[0][0] * offsets[1][1] * offsets[2][2] \
            + offsets[0][1] * offsets[1][2] * offsets[2][0] \
            + offsets[0][2] * offsets[1][0] * offsets[2][1] \
            - offsets[0][0] * offsets[1][2] * offsets[2][1] \
            - offsets[0][1] * offsets[1][0] * offsets[2][2] \
            - offsets[0][2] * offsets[1][1] * offsets[2][0]

    base_volume = len(data["original_block"]) * len(data["base_blocks"])
    if determinant != base_volume:
        # this can not fit
        print(f"determinant:{determinant}, expected:{base_volume}")
        sys.exit(1)

    # For each cell in base_blocks, calculate its coordinates in the base
    # defined by offsets, these will be some integer / determinant.
    # If we get only the numerator modulo determinant, these should all
    # be unique triplets.

    # inverse of transposed of offsets matrix, multiplied by determinant
    inverse = [
            [
                offsets[1][1] * offsets[2][2] - offsets[1][2] * offsets[2][1],
                offsets[1][2] * offsets[2][0] - offsets[1][0] * offsets[2][2],
                offsets[1][0] * offsets[2][1] - offsets[1][1] * offsets[2][0]
            ],
            [
                offsets[0][2] * offsets[2][1] - offsets[0][1] * offsets[2][2],
                offsets[0][0] * offsets[2][2] - offsets[0][2] * offsets[2][0],
                offsets[0][1] * offsets[2][0] - offsets[0][0] * offsets[2][1]

            ],
            [
                offsets[0][1] * offsets[1][2] - offsets[0][2] * offsets[1][1],
                offsets[0][2] * offsets[1][0] - offsets[0][0] * offsets[1][2],
                offsets[0][0] * offsets[1][1] - offsets[0][1] * offsets[1][0]
            ]
        ]

    occupied = {}
    for b in data["base_blocks"]:
        for x, y, z in b:
            x1 = (inverse[0][0] * x + inverse[0][1] * y + inverse[0][2] * z) % determinant
            y1 = (inverse[1][0] * x + inverse[1][1] * y + inverse[1][2] * z) % determinant
            z1 = (inverse[2][0] * x + inverse[2][1] * y + inverse[2][2] * z) % determinant
            if (x1, y1, z1) in occupied:
                print(f"Collision between ({x}, {y}, {z}) and ({occupied[(x1, y1, z1)]})")
                sys.exit(1)
            occupied[(x1, y1, z1)] = (x, y, z)

    print("OK")
    sys.exit(0)

if __name__ == "__main__":
    main()
