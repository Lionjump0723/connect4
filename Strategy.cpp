#include <chrono>
#include <cmath>
#include <iostream>
#include <tuple>
#include <cassert>
#include "Strategy.h"

// --- GTable ---

uint32_t GTable::memory_Byte() const {
    uint32_t s = sizeof(GTable);
    s += sizeof(uint32_t) * hash_table.capacity();
    s += sizeof(Node) * node_info.capacity();
    s += sizeof(Plate::Key) * node_key.capacity();
    s += sizeof(Edge) * edge_info.capacity();
    s += sizeof(Move) * edge_act.capacity();
    return s;
}

GTable::GTable(uint32_t cap) : root_id(REF_NULL), hash_table(cap, REF_NULL) {
    node_info.reserve(kNodeReserve);
    node_key.reserve(kNodeReserve);
    edge_info.reserve(kEdgeReserve);
    edge_act.reserve(kEdgeReserve);
}

uint32_t GTable::find_slot(const Plate::Key& key) {
    uint32_t h = std::hash<Plate::Key>{}(key) & (hash_table.size() - 1);
    while (hash_table[h] != REF_NULL) {
        if (node_key[hash_table[h]] == key) {
            break;
        }
        h += 1;
        h &= (hash_table.size() - 1);
    }
    return h;
}

void GTable::copy_node(uint32_t dst_id, const GTable& src_table, uint32_t src_id) {
    node_info[dst_id] = src_table.node_info[src_id];
    node_key[dst_id] = src_table.node_key[src_id];
    const uint8_t e_beg = src_table.node_info[src_id].e_beg;
    const uint8_t e_end = src_table.node_info[src_id].e_end;
    if (e_end != 0) {
        uint32_t src_offset = src_table.node_info[src_id].e_offset;
        uint32_t dst_offset = edge_info.size();
        node_info[dst_id].e_offset = dst_offset;
        for (Move i = 0; i < e_end; i++) {
            edge_info.push_back(src_table.edge_info[src_offset + i]);
            edge_act.push_back(src_table.edge_act[src_offset + i]);
            uint32_t src_ch_id = edge_info.back().pt;
            if (src_ch_id != REF_NULL) {
                const Plate::Key& src_ch_key = src_table.node_key[src_ch_id];
                uint32_t dst_ch_slot = find_slot(src_ch_key);
                if (hash_table[dst_ch_slot] == REF_NULL) {
                    hash_table[dst_ch_slot] = node_info.size();
                    node_info.emplace_back();
                    node_key.emplace_back();
                    node_info.back().S = STALE;
                }
                edge_info.back().pt = hash_table[dst_ch_slot];
            }
        }
        for (Move i = 0; i < e_end; i++) {
            uint32_t ch_id = edge_info[dst_offset + i].pt;
            if (ch_id != REF_NULL && node_info[ch_id].S == STALE) {
                copy_node(ch_id, src_table, src_table.edge_info[src_offset + i].pt);
            }
        }
    }
}

void GTable::copy_table(const GTable& src_table) {
    assert(root_id == REF_NULL);
    if (src_table.root_id != REF_NULL) {
        uint32_t root_slot = find_slot(src_table.node_key[src_table.root_id]);
        assert(hash_table[root_slot] == REF_NULL);
        root_id = node_info.size();
        hash_table[root_slot] = node_info.size();
        node_info.emplace_back();
        node_key.emplace_back();
        copy_node(root_id, src_table, src_table.root_id);
    }
}

uint32_t GTable::update(uint32_t nid) {
    Node& node = node_info[nid];
    if (node.S == CERTAIN) {
        return 0;
    }
    assert(node.e_end > 0);
    uint32_t sum_edge_N = 0;
    double sum_edge_NxQ = 0;
    double sum_edge_P = 0;
    double sum_edge_PxQ = 0;
    for (Move i = node.e_beg; i < node.e_end; i++) {
        Edge& e = edge_info[node.e_offset + i];
        assert(e.pt != REF_NULL);
        Node& ch = node_info[e.pt];
        if (ch.S == CERTAIN) {
            if (ch.Q < 0) {
                std::swap(edge_info[node.e_offset + i], edge_info[node.e_offset]);
                std::swap(edge_act[node.e_offset + i], edge_act[node.e_offset]);
                node.e_beg = node.e_end;
                break;
            } else {
                std::swap(edge_info[node.e_offset + i], edge_info[node.e_offset + node.e_beg]);
                std::swap(edge_act[node.e_offset + i], edge_act[node.e_offset + node.e_beg]);
                if (node.e_beg > 0) {
                    Edge& best_e = edge_info[node.e_offset];
                    Node& best_ch = node_info[best_e.pt];
                    if (ch.Q < best_ch.Q) {
                        std::swap(edge_info[node.e_offset], edge_info[node.e_offset + node.e_beg]);
                        std::swap(edge_act[node.e_offset], edge_act[node.e_offset + node.e_beg]);
                    }
                }
                node.e_beg += 1;
            }
        } else {
            sum_edge_N += e.N;
            sum_edge_NxQ += e.N * ch.Q;
            sum_edge_P += 1.0 / node.e_end;
            sum_edge_PxQ += 1.0 / node.e_end * ch.Q;
        }
    }
    if (node.e_beg == node.e_end) {
        node.S = CERTAIN;
        Edge& best_edge = edge_info[node.e_offset];
        Node& best_ch = node_info[best_edge.pt];
        node.Q = -best_ch.Q;
        return 0;
    } else {
        if (sum_edge_N == 0) {
            node.U = -sum_edge_PxQ / sum_edge_P;
        }
        Edge& best_edge = edge_info[node.e_offset];
        Node& best_ch = node_info[best_edge.pt];
        if (node.e_beg > 0) {
            node.U = std::max(node.U, (float)-best_ch.Q);
        }
        node.Q = (-sum_edge_NxQ + node.U) / (sum_edge_N + 1);
        if (node.e_beg > 0) {
            node.Q = std::max(node.Q, -best_ch.Q);
        }
        return sum_edge_N + 1;
    }
}

uint32_t GTable::get_best_edge(uint32_t nid) const {
    double max_gain = -2;
    uint32_t ret = REF_NULL;
    const Node& node = node_info[nid];
    uint32_t N = 1;
    for (Move i = node.e_beg; i < node.e_end; i++) {
        const Edge& e = edge_info[node.e_offset + i];
        N += e.N;
    }
    double param = std::sqrt(N) * c_puct / node.e_end;
    for (Move i = node.e_beg; i < node.e_end; i++) {
        const Edge& e = edge_info[node.e_offset + i];
        assert(e.pt != REF_NULL);
        const Node& ch = node_info[e.pt];
        double gain;
        if (ch.S == CERTAIN) {
            gain = -ch.Q;
        } else {
            gain = param / (e.N + 1) + (ch.e_end == 0 ? (node.U - ch.Q) / 2 : -ch.Q);
        }
        if (max_gain < gain) {
            max_gain = gain;
            ret = node.e_offset + i;
        }
    }
    return ret;
}

uint32_t GTable::search(uint32_t nid, const BoardConfig& cfg, const Plate& root) {
    Move all[16];
    Node& node = node_info[nid];
    if (node.S == CERTAIN) {
        return 0;
    }
    if (node.e_end > 0) {
        uint32_t eid = get_best_edge(nid);
        uint32_t ch_n = search(edge_info[eid].pt, cfg, root);
        if (edge_info[eid].N < ch_n) {
            edge_info[eid].N += 1;
        }
    } else {
        Plate pl = Plate::get_plate(node_key[nid], cfg, root);
        ExInfo exi = pl.build(cfg);
        double v = pl.get_move(exi, all, node.e_end, cfg);
        if (v != kMoveUnknown) {
            node.S = CERTAIN;
            node.Q = v;
            node.e_offset = all[0];
            node.e_beg = 0;
            node.e_end = 0;
        } else {
            node.e_offset = edge_info.size();
            const Move range = node.e_end;
            for (Move i = 0; i < range; i++) {
                edge_info.emplace_back();
                edge_act.push_back(all[i]);
                Plate tmp_pl = pl;
                ExInfo tmp_exi = exi;
                tmp_pl.step_exi(all[i], tmp_exi, cfg);
                Plate::Key ch_key = tmp_pl.get_key();
                uint32_t ch_slot = find_slot(ch_key);
                if (hash_table[ch_slot] == REF_NULL) {
                    hash_table[ch_slot] = node_info.size();
                    node_info.emplace_back();
                    node_key.push_back(ch_key);
                    Move mv = MV_NULL;
                    std::tie(mv, node_info.back().Q) = tmp_pl.forward_check(tmp_exi, cfg);
                    if (mv != MV_NULL) {
                        node_info.back().S = CERTAIN;
                        node_info.back().e_offset = mv;
                    }
                }
                edge_info.back().pt = hash_table[ch_slot];
            }
        }
    }
    return update(nid);
}

// --- Search helpers ---

using Clock = std::chrono::high_resolution_clock;

static double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

struct RootPick {
    int y = 0;
    uint32_t root_count = 0;
};

static void init_search_graph(GTable& g, const Plate& root) {
    Plate::Key root_key = root.get_key();
    g.root_id = g.hash_table[g.find_slot(root_key)];
    GTable new_graph(kHashTableCapacity);
    if (g.root_id != REF_NULL) {
        new_graph.copy_table(g);
    } else {
        uint32_t slot = new_graph.find_slot(root_key);
        new_graph.root_id = 0;
        new_graph.hash_table[slot] = 0;
        new_graph.node_info.emplace_back();
        new_graph.node_key.push_back(root_key);
    }
    PRINT_CERR << "memory peak: "
               << (g.memory_Byte() + new_graph.memory_Byte() + 0.0) / 1048576 << "MB"
               << std::endl;
    std::swap(g, new_graph);
}

static void run_search_loop(GTable& g, const BoardConfig& cfg, const Plate& root,
                            const std::chrono::high_resolution_clock::time_point& t0) {
    while (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count() <
               kSearchTimeSec &&
           g.node_info.size() <= kMaxGraphNodes && g.search(g.root_id, cfg, root)) {
    }
}

static RootPick select_root_column(GTable& g) {
    Node& root_node = g.node_info[g.root_id];
    const uint32_t root_count = g.update(g.root_id);
    int y = 0;
    if (root_node.S == CERTAIN) {
        if (root_node.e_end > 0) {
            y = g.edge_act[root_node.e_offset];
        } else {
            y = root_node.e_offset;
        }
        PRINT_CERR << "certain!" << root_node.Q << " " << y << std::endl;
    } else {
        uint32_t max_n = 0;
        double ch_Q = -2;
        for (Move i = root_node.e_beg; i < root_node.e_end; i++) {
            Edge& e = g.edge_info[root_node.e_offset + i];
            if (e.N >= max_n / 4 && -g.node_info[e.pt].Q > ch_Q) {
                max_n = e.N;
                y = g.edge_act[root_node.e_offset + i];
                ch_Q = -g.node_info[e.pt].Q;
            }
        }
        if (root_node.e_beg > 0) {
            Edge& best_e = g.edge_info[root_node.e_offset];
            double best_Q = -g.node_info[best_e.pt].Q;
            if (ch_Q < best_Q) {
                y = g.edge_act[root_node.e_offset];
                ch_Q = best_Q;
            }
        }
        PRINT_CERR << ch_Q << " " << y << std::endl;
    }
    return {y, root_count};
}

// --- SearchContext ---

SearchContext::SearchContext() : graph(kHashTableCapacity) {}

void SearchContext::reset() {
    board.reset();
    root = Plate();
    graph = GTable(kHashTableCapacity);
}

void SearchContext::load_position(int M, int N, int noX, int noY, const int* raw_board) {
    board.configure(M, N, noX, noY);
    root = make_root_plate(board, raw_board);
    init_search_graph(graph, root);
}

Point* SearchContext::get_point(const int* top,
                                const std::chrono::high_resolution_clock::time_point& t0) {
    PRINT_CERR << "init_time: "
               << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0)
                      .count()
               << std::endl;
    PRINT_CERR << "begin graph size: node " << graph.node_info.size() << ", edge "
               << graph.edge_info.size() << std::endl;

    run_search_loop(graph, board, root, t0);

    PRINT_CERR << "end graph size: node " << graph.node_info.size() << ", edge "
               << graph.edge_info.size() << std::endl;

    const RootPick pick = select_root_column(graph);
    PRINT_CERR << "root count: " << pick.root_count << std::endl;
    return new Point(top[pick.y] - 1, pick.y);
}

SearchBenchmarkResult SearchContext::benchmark_position(int M, int N, int noX, int noY,
                                                        const int* top, const int* raw_board) {
    SearchBenchmarkResult result;
    result.M = M;
    result.N = N;
    result.noX = noX;
    result.noY = noY;

    const auto t_init_start = Clock::now();
    board.configure(M, N, noX, noY);
    root = make_root_plate(board, raw_board);
    const auto t_init_end = Clock::now();
    result.init_ms = elapsed_ms(t_init_start, t_init_end);

    const auto t_copy_start = Clock::now();
    init_search_graph(graph, root);
    const auto t_copy_end = Clock::now();
    result.copy_ms = elapsed_ms(t_copy_start, t_copy_end);

    result.nodes_before = static_cast<uint32_t>(graph.node_info.size());
    result.edges_before = static_cast<uint32_t>(graph.edge_info.size());

    const auto t_search_start = Clock::now();
    run_search_loop(graph, board, root, t_search_start);
    const auto t_search_end = Clock::now();
    result.search_ms = elapsed_ms(t_search_start, t_search_end);

    result.nodes_after = static_cast<uint32_t>(graph.node_info.size());
    result.edges_after = static_cast<uint32_t>(graph.edge_info.size());

    const auto t_sel_start = Clock::now();
    const RootPick pick = select_root_column(graph);
    const auto t_sel_end = Clock::now();
    result.selection_ms = elapsed_ms(t_sel_start, t_sel_end);
    result.root_count = pick.root_count;
    result.chosen_y = pick.y;
    result.chosen_x = top[pick.y] - 1;

    return result;
}

extern "C" Point* getPoint(const int M, const int N, const int* top, const int* _board,
                           const int lastX, const int lastY, const int noX, const int noY) {
    (void)lastX;
    (void)lastY;

    thread_local SearchContext ctx;
    auto t0 = std::chrono::high_resolution_clock::now();
    ctx.load_position(M, N, noX, noY, _board);
    return ctx.get_point(top, t0);
}

extern "C" void clearPoint(Point* p) {
    delete p;
}
