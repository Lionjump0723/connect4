/********************************************************
 * Strategy.h : 策略接口与 MCGS 搜索图声明
 *********************************************************/

#ifndef STRATEGY_H_
#define STRATEGY_H_

#include <chrono>
#include <cstdint>
#include <random>
#include <vector>
#include "Point.h"
#include "Config.h"
#include "Logic.h"

extern "C" Point* getPoint(const int M, const int N, const int* top, const int* _board,
                           const int lastX, const int lastY, const int noX, const int noY);

extern "C" void clearPoint(Point* p);

struct Node {
    double Q = -2;
    float U = 0;
    uint32_t e_offset = REF_NULL;
    Move e_end = 0;
    Move e_beg = 0;
    Status S = 0;
};

struct Edge {
    uint32_t pt = REF_NULL;
    uint32_t N = 0;
};

// MCGS graph: nodes/edges in vectors, Zobrist-keyed linear probing hash
struct GTable {
    std::vector<uint32_t> hash_table;
    std::vector<Node> node_info;
    std::vector<Plate::Key> node_key;
    std::vector<Edge> edge_info;
    std::vector<Move> edge_act;
    uint32_t root_id;

    uint32_t memory_Byte() const;
    explicit GTable(uint32_t cap);

    uint32_t find_slot(const Plate::Key& key);
    void copy_node(uint32_t dst_id, const GTable& src_table, uint32_t src_id);
    void copy_table(const GTable& src_table);
    uint32_t update(uint32_t nid);
    uint32_t get_best_edge(uint32_t nid) const;
    uint32_t search(uint32_t nid, const BoardConfig& cfg, const Plate& root, std::mt19937& rnd);
};

struct SearchBenchmarkResult {
    int M = 0;
    int N = 0;
    int noX = -1;
    int noY = -1;
    int chosen_x = -1;
    int chosen_y = -1;
    double init_ms = 0;
    double copy_ms = 0;
    double search_ms = 0;
    double selection_ms = 0;
    uint32_t nodes_before = 0;
    uint32_t edges_before = 0;
    uint32_t nodes_after = 0;
    uint32_t edges_after = 0;
    uint32_t root_count = 0;
};

struct SearchContext {
    BoardConfig board;
    Plate root;
    GTable graph;
    std::mt19937 rnd;

    SearchContext();
    void reset();
    void load_position(int M, int N, int noX, int noY, const int* raw_board);
    Point* get_point(const int* top,
                     const std::chrono::high_resolution_clock::time_point& t0);
    SearchBenchmarkResult benchmark_position(int M, int N, int noX, int noY, const int* top,
                                             const int* raw_board);
};

#endif
