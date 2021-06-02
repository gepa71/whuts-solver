#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int dim;
int poly_size;
int cnt_bases;
int volume;
int num_polyominos;

#define MAX_POLY_SIZE 15
#define MAX_DIM 4
// MAX_ROTATIONS = (MAX_DIM)! * 2 ^(MAX_DIM - 1)
#define MAX_ROTATIONS 192

#define MAX_VOLUME 200
// upper bound for volume of a rectangular prism of vorume MAX_VOLUME
// padded by MAX_POLY_SIZE-1 in each of the MAX_DIM dimensions.
#define MAX_PADDED_VOLUME 10000000

#define MAX_ITEMS (MAX_VOLUME * MAX_ROTATIONS * MAX_POLY_SIZE + 1)

static int INT_SQRT[MAX_VOLUME + 1];

static int rotations[MAX_ROTATIONS][MAX_DIM][MAX_DIM];
static int offsets[MAX_DIM][MAX_DIM];
static int rotated_polyominos[MAX_ROTATIONS][MAX_POLY_SIZE][MAX_DIM];
static int basic_polyomino[MAX_POLY_SIZE][MAX_DIM];

static int num_rotations;
static int num_poly_rotations;

int current_rotation[MAX_DIM][MAX_DIM];
int used[MAX_DIM];

static int cache_origins[MAX_PADDED_VOLUME];
static int check_bitmap[MAX_VOLUME];

static int next_up[MAX_ITEMS][MAX_VOLUME + 1];
static int next_down[MAX_ITEMS][MAX_VOLUME + 1];
static int next_left[MAX_ITEMS][MAX_VOLUME + 1];
static int next_right[MAX_ITEMS][MAX_VOLUME + 1];
static int next_count[MAX_VOLUME + 1];

static int solution[MAX_VOLUME];

void calculate_rotations(int n, int parity) {
    if (n == dim) {
        if (parity % 2 == 0) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    rotations[num_rotations][i][j] = current_rotation[i][j];
                }
            }
            num_rotations++;
        }
        return;
    }
    int p = parity;
    for (int i = 0; i < dim; i++) {
        if (used[i]) {
            continue;
        }
        used[i] = 1;
        current_rotation[n][i] = 1;
        calculate_rotations(n + 1, p);
        current_rotation[n][i] = -1;
        calculate_rotations(n + 1, p + 1);
        current_rotation[n][i] = 0;
        used[i] = 0;
        p++;
    }
}

void rotate(const int poly[][MAX_DIM], const int rot[][MAX_DIM], int dest[][MAX_DIM]) {
    for (int i = 0; i < poly_size; i++) {
        for (int j = 0; j < dim; j++) {
            dest[i][j] = 0;
            for (int k = 0; k < dim; k++) {
                dest[i][j] += rot[j][k] * poly[i][k];
            }
        }
    }
}

int is_less(const int *cube1, const int *cube2) {
    for (int i = 0; i < dim; i++) {
        if (cube1[i] != cube2[i]) {
            return cube1[i] < cube2[i];
        }
    }
    return 0;
}

void poly_quick_sort(int poly[][MAX_DIM], int start, int end) {
    if (start + 1 >= end) {
        return;
    }
    int pivot[MAX_DIM];
    for (int i = 0; i < dim; i++) {
        pivot[i] = poly[end-1][i];
    }
    int front = start;
    int back = end - 1;
    while (back > front) {
        if (is_less(poly[back-1], pivot)) {
            for (int j = 0; j < dim; j++) {
                int tmp = poly[front][j];
                poly[front][j] = poly[back-1][j];
                poly[back-1][j] = tmp;
            }
            front++;
        } else {
            for (int j = 0; j < dim; j++) {
                poly[back][j] = poly[back-1][j];
            }
            back--;
        }

    }
    for (int j = 0; j < dim; j++) {
        poly[front][j] = pivot[j];
    }
    poly_quick_sort(poly, start, front);
    poly_quick_sort(poly, front + 1, end);
}

void norm(int poly[][MAX_DIM]) {
    int dmin[MAX_DIM];
    for (int i = 0; i < dim; i++) {
        dmin[i] = poly[0][i];
    }
    for (int j = 0; j < poly_size; j++) {
        for (int i = 0; i < dim; i++) {
            if (poly[j][i] < dmin[i]) {
                dmin[i] = poly[j][i];
            }
        }
    }
    for (int j = 0; j < poly_size; j++) {
        for (int i = 0; i < dim; i++) {
            poly[j][i] -= dmin[i];
        }
    }
    poly_quick_sort(poly, 0, poly_size);
}

void init_sqrt() {
    int a = 1;
    for (int i = 1; i <= MAX_VOLUME; i++) {
        if ((a + 1) * (a + 1) <= i) {
            a++;
        }
        INT_SQRT[i] = a;
    }
}

void init_rotations() {
    num_rotations = 0;
    for (int i = 0; i < dim; i++) {
        used[i] = 0;
        for (int j = 0; j < dim; j++) {
            current_rotation[i][j] = 0;
        }
    }
    calculate_rotations(0, 0);
    num_poly_rotations = 0;
    for (int r = 0; r < num_rotations; r++) {
        rotate(basic_polyomino, rotations[r], rotated_polyominos[num_poly_rotations]);
        norm(rotated_polyominos[num_poly_rotations]);
        int same = 0;
        for (int j = 0; j < num_poly_rotations; j++) {
            same = 1;
            for (int k = 0; same && (k < poly_size); k++) {
                for (int d = 0; d < dim; d++) {
                    if (rotated_polyominos[j][k][d] != rotated_polyominos[num_poly_rotations][k][d]) {
                        same = 0;
                        break;
                    }
                }
            }
            if (same) {
                break;
            }
        }
        if (!same) {
            num_poly_rotations++;
        }
    }
}

void init_offsets() {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            offsets[i][j] = 0;
        }
    }
}

void int_to_coords(int p, int coords[]) {
    for (int i = dim - 1; i >= 0; i--) {
        coords[i] = p % offsets[i][i];
        p /= offsets[i][i];
    }
}

int coords_to_int(const int coords[]) {
    int p = 0;
    for (int i = 0; i < dim; i++) {
        p *= offsets[i][i];
        p += coords[i];
    }
    return p;
}

int padded_coords_to_int(const int coords[]) {
    int p = 0;
    for (int i = 0; i < dim; i++) {
        p *= offsets[i][i] + 2 * poly_size - 2;
        p += coords[i] + poly_size - 1;
    }
    return p;
}

int transform_to_origin(int coords[]) {
    // transform coords to the rectangular prism at the origin
    // of the coordinate system.
    // Returns an integer representation of one of the <volume> points.
    int padded_location = padded_coords_to_int(coords);
    if (padded_location < 0 || padded_location >= MAX_PADDED_VOLUME) {
        printf("padded_location:%d\n", padded_location);
    }
    int *cache = &cache_origins[padded_location];
    if (*cache >= 0) {
        return *cache;
    }
    for (int i = 0; i < dim; i++) {
        int mi = coords[i] % offsets[i][i];
        int di;
        if (mi < 0) {
            di = coords[i] / offsets[i][i] - 1;
        } else {
            di = coords[i] / offsets[i][i];
        }
        for (int j = i; j < dim; j++) {
            coords[j] -= di * offsets[i][j];
        }
    }
    *cache = coords_to_int(coords);
    if (*cache < 0 || *cache >= volume) {
        printf("bad location: %d\n", *cache);
    }
    return *cache;
}

void print_block(const int poly[][MAX_DIM], const int origin_coords[], int idx, int is_last) {
    printf("[");
    for (int i = 0; i < poly_size; i++) {
        int cube_coords[MAX_DIM];
        printf("[");
        for (int j = 0; j < dim; j++) {
            cube_coords[j] = origin_coords[j] + poly[i][j] - poly[idx][j];
            printf("%d", cube_coords[j]);
            if (j + 1 < dim) {
                printf(",");
            } else {
                printf("]");
            }
        }
        if (i + 1 < poly_size) {
            printf(",");
        }
    }
    if (is_last) {
        printf("]\n");
    } else {
        printf("],\n");
    }
}

void print_offsets() {
    printf("[");
    for (int i = 0; i < dim; i++) {
        printf("[");
        for (int j = 0; j < dim; j++) {
            printf("%d", offsets[i][j]);
            if (j + 1 < dim) {
                printf(",");
            }
        }
        printf("]");
        if (i + 1 < dim) {
            printf(",");
        }
    }
    printf("]");
}

void print_solution() {
    int primitive_size = volume / poly_size;
    printf("{\n    \"original_block\": ");
    print_block(basic_polyomino, basic_polyomino[0], 0, 0);
    printf("    \"base_blocks\": [\n");
    for (int i = 0; i < primitive_size; i++) {
        int r = solution[i];
        int idx = r % poly_size;
        r /= poly_size;
        int origin = r % volume;
        r /= volume;
        int origin_coords[MAX_DIM];
        int_to_coords(origin, origin_coords);
        printf("        ");
        print_block(rotated_polyominos[r], origin_coords, idx, i + 1 == primitive_size);
    }
    printf("    ],\n    \"offsets\": ");
    print_offsets();
    printf("\n}\n");
}

int get_covered_blocks(const int poly[][MAX_DIM], const int origin_coords[], int idx, int covered_blocks[], int cnt) {
    // position the <idx>-th hypercube of <poly> at <origin>
    for (int i = 0; i < poly_size; i++) {
        // cube_coords: where the i-th hypercube of poly ends up
        int cube_coords[MAX_DIM];
        for (int j = 0; j < dim; j++) {
            cube_coords[j] = origin_coords[j] + poly[i][j] - poly[idx][j];
        }
        int b = transform_to_origin(cube_coords);
        if (check_bitmap[b] == cnt) {
            return 0;
        }
        covered_blocks[i] = b;
        check_bitmap[b] = cnt;
    }
    return 1;
}

void delete_up_down(int r, int c) {
    next_down[next_up[r][c]][c] = next_down[r][c];
    next_up[next_down[r][c]][c] = next_up[r][c];
}

void delete_left_right(int r, int c) {
    next_right[r][next_left[r][c]] = next_right[r][c];
    next_left[r][next_right[r][c]] = next_left[r][c];
}

void recover_up_down(int r, int c) {
    next_down[next_up[r][c]][c] = r;
    next_up[next_down[r][c]][c] = r;
}

void recover_left_right(int r, int c) {
    next_right[r][next_left[r][c]] = c;
    next_left[r][next_right[r][c]] = c;
}


int solve(int n) {
    if (next_right[0][0] == 0) {
        return 1;
    }
    int c = next_right[0][0];
    int cmin = next_right[0][0];
    int cnt_min = next_count[cmin];
    for (c = next_right[0][0]; c != 0; c = next_right[0][c]) {
        if (next_count[c] < cnt_min) {
            cnt_min = next_count[c];
            cmin = c;
        }
    }
    for (int row = next_down[0][cmin]; row != 0; row = next_down[row][cmin]) {
        solution[n] = row - 1;
        int del_rows[MAX_ROTATIONS * MAX_POLY_SIZE * MAX_POLY_SIZE];
        int num_del_rows = 0;
        for (int col = next_right[row][0]; col != 0; col = next_right[row][col]) {
            for (int drow = next_down[row][col]; drow != row; drow = next_down[drow][col]) {
                if (drow == 0) {
                    continue;
                }
                del_rows[num_del_rows++] = drow;
                for (int dcol = next_right[drow][0]; dcol != 0; dcol = next_right[drow][dcol]) {
                    delete_up_down(drow, dcol);
                    next_count[dcol]--;
                }
            }
            delete_left_right(0, col);
        }

        if (solve(n + 1)) {
            return 1;
        }

        for (int col = next_right[row][0]; col != 0; col = next_right[row][col]) {
            // bring back col
            recover_left_right(0, col);
        }
        for (int i = num_del_rows - 1; i >= 0; i--) {
            int drow = del_rows[i];
            // bring back drow
            for (int dcol = next_right[drow][0]; dcol != 0; dcol = next_right[drow][dcol]) {
                recover_up_down(drow, dcol);
                next_count[dcol]++;
            }
        }
        if (n == 0) {
            return 0;
        }
    }
    return 0;
}

int find_tiling() {
    // c version of "yield"...
    // Here, offsets[][] has been setup by the
    // recursive calls to check_next_diagonal() and check_next_offsets()
#ifdef DEBUG
    cnt_bases++;
    if (cnt_bases % 1000 == 0) {
        printf("cnt_bases:%d offsets:", cnt_bases);
        print_offsets();
        printf("\n");
    }
#endif
    volume = 1;
    int volume_padded = 1;
    int covered_blocks[MAX_POLY_SIZE];
    for (int i = 0; i < dim; i++) {
        volume *= offsets[i][i];
        volume_padded *= offsets[i][i] + 2 * MAX_POLY_SIZE - 2;
    }
    if (volume_padded >= MAX_PADDED_VOLUME) {
        fprintf(stderr, "volume_padded exceeded max allowed:%d\n", volume_padded);
        fprintf(stderr, "Please change MAX_PADDED_VOLUME in source and recompile.\n");
        exit(1);
    }
    memset(cache_origins, -1, volume_padded * sizeof(int));
    memset(check_bitmap, 0, volume * sizeof(int));
    for (int i = 1; i <= volume; i++) {
        next_up[0][i] = 0;
        next_down[0][i] = 0;
        next_count[i] = 0;
        next_left[0][i] = i - 1;
        next_right[0][i] = (i + 1) % (volume + 1);
    }
    next_left[0][0] = volume;
    next_right[0][0] = 1;
    int cnt = 0;
    for (int r = 0; r < num_poly_rotations; r++) {
        for (int origin = 0; origin < volume; origin++) {
            int origin_coords[MAX_DIM];
            int_to_coords(origin, origin_coords);
            for (int idx = 0; idx < poly_size; idx++) {
                cnt++;
                if (cnt >= MAX_ITEMS) {
                    printf("cnt reached MAX_ITEMS!\n");
                    return 0;
                }
                if (!get_covered_blocks(rotated_polyominos[r], origin_coords, idx, covered_blocks, cnt)) {
                    continue;
                }
                next_left[cnt][0] = 0;
                next_right[cnt][0] = 0;
                for (int i = 0; i < poly_size; i++) {
                    int c = covered_blocks[i] + 1;
                    next_count[c]++;

                    next_up[cnt][c] = next_up[0][c];
                    next_down[cnt][c] = 0;
                    next_down[next_up[0][c]][c] = cnt;
                    next_up[0][c] = cnt;

                    next_left[cnt][c] = next_left[cnt][0];
                    next_right[cnt][c] = 0;
                    next_right[cnt][next_left[cnt][0]] = c;
                    next_left[cnt][0] = c;

                }
            }
        }
    }

    if (solve(0)) {
        print_solution();
        return 1;
    }
    return 0;
}

int check_next_offsets(int n, int k) {
    if (n == dim) {
        return find_tiling();
    }
    if (k == dim) {
        return check_next_offsets(n + 1, n + 2);
    }
    for (int i = 0; i < offsets[k][k]; i++) {
        offsets[n][k] = i;
        if (check_next_offsets(n, k + 1)) {
            return 1;
        }
    }
    return 0;
}

int check_next_diagonal(int remaining_volume, int n) {
    if (n == dim - 1) {
        offsets[n][n] = remaining_volume;
        return check_next_offsets(0, 1);
    }
    for (int i = INT_SQRT[remaining_volume]; i > 0; i--) {
        if (remaining_volume % i != 0) {
            continue;
        }
        offsets[n][n] = i;
        if (check_next_diagonal(remaining_volume / i, n + 1)) {
            return 1;
        }
        if (i * i != remaining_volume) {
            offsets[n][n] = remaining_volume / i;
            if (check_next_diagonal(i, n + 1)) {
                return 1;
            }
        }
    }
    return 0;
}

void usage(char *cmd) {
    fprintf(stderr, "Usage: %s <number of polyominoes in basic pattern>\n", cmd);
    fprintf(stderr, "    Standard input needs to supply dimension, polyomino size\n");
    fprintf(stderr, "    and coordinates of the polyomino separated by whitespace\n");
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        usage(argv[0]);
        exit(1);
    }
    cnt_bases = 0;
    num_polyominos = atoi(argv[1]);
    if (!scanf("%d", &dim)) {
        fprintf(stderr, "First element in input needs to be the dimension\n");
        usage(argv[0]);
        exit(1);
    }
    if (dim < 2 || dim > MAX_DIM) {
        fprintf(stderr, "Dimension was %d. Must be between 2 and %d\n", dim, MAX_DIM);
        usage(argv[0]);
        exit(1);
    }
    if (!scanf("%d", &poly_size)) {
        fprintf(stderr, "Second element in input needs to be the polyomino size\n");
        usage(argv[0]);
        exit(1);
    }
    if (poly_size > MAX_POLY_SIZE) {
        fprintf(stderr, "poly_size (second element in input) must be <%d, was %d\n", MAX_POLY_SIZE, poly_size);
        fprintf(stderr, "Please update MAX_POLY_SIZE in source and recompile\n");
        exit(1);
    }

    for (int i = 0; i < poly_size; i++) {
        for (int j = 0; j < dim; j++) {
            int p;
            if (!scanf("%d", &p)) {
                fprintf(stderr, "Unexpected end of input, needs %d integers for the polyomino, got %d\n",
                        poly_size * dim, i * dim + j);
                usage(argv[0]);
                exit(1);
            }
            basic_polyomino[i][j] = p;
        }
    }
    init_rotations();
    init_sqrt();
    init_offsets();
    if (num_polyominos * poly_size > MAX_VOLUME) {
        fprintf(stderr, "volume of basic pattern (num_polyominos * poly_size) must be less than MAX_VOLUME\n");
        fprintf(stderr, "num_polyominos=%d poly_size=%d MAX_VOLUME=%d\n", num_polyominos, poly_size, MAX_VOLUME);
        fprintf(stderr, "Please update MAX_VOLUME in source code and recompile\n");
        exit(1);
    }
    check_next_diagonal(num_polyominos * poly_size, 0);
}
