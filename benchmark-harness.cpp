#include <cstdio>
#include <cstring>
#include <vector>
#include "Strategy.h"

struct BenchCase {
    char name[64];
    int M;
    int N;
    int noX;
    int noY;
    std::vector<int> board;
    std::vector<int> top;
};

static void skip_forbidden_top(BenchCase& c, int col) {
    if (col == c.noY && c.top[col] > 0 && c.top[col] - 1 == c.noX) {
        c.top[col] -= 1;
    }
}

static bool drop_stone(BenchCase& c, int col, int player) {
    if (col < 0 || col >= c.N) {
        return false;
    }
    skip_forbidden_top(c, col);
    if (c.top[col] <= 0) {
        return false;
    }
    const int row = c.top[col] - 1;
    if (row < 0 || row >= c.M || c.board[static_cast<size_t>(row * c.N + col)] != 0) {
        return false;
    }
    c.board[static_cast<size_t>(row * c.N + col)] = player;
    c.top[col] -= 1;
    skip_forbidden_top(c, col);
    return true;
}

static void apply_columns(BenchCase& c, const std::vector<int>& cols) {
    int player = 1;
    for (int col : cols) {
        if (drop_stone(c, col, player)) {
            player = (player == 1) ? 2 : 1;
        }
    }
}

static BenchCase make_empty_case(int M, int N) {
    BenchCase c{};
    c.M = M;
    c.N = N;
    c.noX = M / 2;
    c.noY = N / 2;
    c.board.assign(static_cast<size_t>(M * N), 0);
    c.top.assign(static_cast<size_t>(N), M);
    skip_forbidden_top(c, c.noY);
    return c;
}

static BenchCase make_opening_case(int M, int N) {
    BenchCase c = make_empty_case(M, N);
    std::snprintf(c.name, sizeof(c.name), "%dx%d-opening", M, N);
    const int center = N / 2;
    apply_columns(c, {center, center - 1, center, center + 1, center, center - 1, center + 1});
    return c;
}

static BenchCase make_midgame_case(int M, int N) {
    BenchCase c = make_empty_case(M, N);
    std::snprintf(c.name, sizeof(c.name), "%dx%d-midgame", M, N);

    std::vector<int> cols;
    const int target_moves = (M * N) / 3;
    for (int i = 0; i < target_moves; i++) {
        int col = (i * 3 + 1) % N;
        if (col == c.noY) {
            col = (col + 1) % N;
        }
        cols.push_back(col);
    }
    apply_columns(c, cols);
    return c;
}

static BenchCase make_blocked_case(int M, int N) {
    BenchCase c = make_empty_case(M, N);
    std::snprintf(c.name, sizeof(c.name), "%dx%d-blocked", M, N);

    const int center = N / 2;
    std::vector<int> cols;
    for (int k = 0; k < M; k++) {
        cols.push_back(center - 1);
        cols.push_back(center + 1);
        if (center - 2 >= 0) {
            cols.push_back(center - 2);
        }
        if (center + 2 < N) {
            cols.push_back(center + 2);
        }
    }
    apply_columns(c, cols);
    return c;
}

static std::vector<BenchCase> build_all_cases() {
    std::vector<BenchCase> cases;
    for (int size = 9; size <= 12; size++) {
        cases.push_back(make_opening_case(size, size));
        cases.push_back(make_midgame_case(size, size));
        cases.push_back(make_blocked_case(size, size));
    }
    return cases;
}

static void print_csv_header() {
    std::printf(
        "case,M,N,noX,noY,chosen_x,chosen_y,init_ms,copy_ms,search_ms,selection_ms,"
        "nodes_before,edges_before,nodes_after,edges_after,root_count\n");
}

static void print_csv_row(const BenchCase& bench, const SearchBenchmarkResult& r) {
    std::printf(
        "%s,%d,%d,%d,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%u,%u,%u,%u,%u\n",
        bench.name, r.M, r.N, r.noX, r.noY, r.chosen_x, r.chosen_y, r.init_ms, r.copy_ms,
        r.search_ms, r.selection_ms, r.nodes_before, r.edges_before, r.nodes_after,
        r.edges_after, r.root_count);
}

int main() {
    const std::vector<BenchCase> cases = build_all_cases();
    print_csv_header();

    SearchContext ctx;
    int prev_m = -1;
    int prev_n = -1;

    for (const BenchCase& bench : cases) {
        if (bench.M != prev_m || bench.N != prev_n) {
            ctx.reset();
            prev_m = bench.M;
            prev_n = bench.N;
        }

        const SearchBenchmarkResult result = ctx.benchmark_position(
            bench.M, bench.N, bench.noX, bench.noY, bench.top.data(), bench.board.data());
        print_csv_row(bench, result);
    }

    return 0;
}
