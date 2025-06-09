/********************************************************
*	Strategy.h : 策略接口文件                           *
*	张永锋                                              *
*	zhangyf07@gmail.com                                 *
*	2014.5                                              *
*********************************************************/

#ifndef STRATEGY_H_
#define	STRATEGY_H_
#include <cstdint>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <tuple>
#include <random>
#include <iostream>
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


struct GTable {
    std::vector<uint32_t> hash_table;
    std::vector<Node> node_info;
    std::vector<Plate::Key> node_key;
    std::vector<Edge> edge_info;
    std::vector<Move> edge_act;
    uint32_t root_id;
    uint32_t memory_Byte() const {
        uint32_t s = sizeof(GTable);
        s += sizeof(uint32_t) * hash_table.capacity();
        s += sizeof(Node) * node_info.capacity();
        s += sizeof(Plate::Key) * node_key.capacity();
        s += sizeof(Edge) * edge_info.capacity();
        s += sizeof(Move) * edge_act.capacity();
        return s;
    }

    GTable(const uint32_t cap): root_id(REF_NULL), hash_table(cap, REF_NULL) {
        node_info.reserve((1<<14) + 10);
        node_key.reserve((1<<14) + 10);
        edge_info.reserve(46875);
        edge_act.reserve(46875);
    }

    uint32_t find_slot(const Plate::Key& key) {
        uint32_t h = std::hash<Plate::Key>{}(key) & (hash_table.size() - 1);
        while (hash_table[h] != REF_NULL) {
            if (node_key[hash_table[h]] == key) break;
            h += 1;
            h &= (hash_table.size() - 1);
        }
        return h;
    }

    void copy_node(uint32_t dst_id, const GTable & src_table, uint32_t src_id) {
        node_info[dst_id] = src_table.node_info[src_id];
        node_key[dst_id] = src_table.node_key[src_id];
        const uint8_t e_beg = src_table.node_info[src_id].e_beg;
        const uint8_t e_end = src_table.node_info[src_id].e_end;
        if (e_end != 0) {
            uint32_t src_offset = src_table.node_info[src_id].e_offset;
            uint32_t dst_offset = edge_info.size();
            node_info[dst_id].e_offset = dst_offset;
            for (Move i = 0; i < e_end; i++) {
                edge_info.push_back(src_table.edge_info[src_offset+i]);
                edge_act.push_back(src_table.edge_act[src_offset+i]);
                uint32_t src_ch_id = edge_info.back().pt;
                if (src_ch_id != REF_NULL) {
                    const Plate::Key & src_ch_key = src_table.node_key[src_ch_id];
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
                uint32_t ch_id = edge_info[dst_offset+i].pt;
                if (ch_id != REF_NULL && node_info[ch_id].S == STALE) {
                    copy_node(ch_id, src_table, src_table.edge_info[src_offset+i].pt);
                }
            }
        }
        
    }
    
    void copy_table(const GTable & src_table) {
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

    uint32_t update(uint32_t nid) {
        Node &node = node_info[nid];
        if (node.S == CERTAIN) return 0;
        assert(node.e_end > 0);
        uint32_t sum_edge_N = 0;
        double sum_edge_NxQ = 0;
        double sum_edge_P = 0;
        double sum_edge_PxQ = 0;
        for (Move i = node.e_beg; i < node.e_end; i++) {
            Edge &e = edge_info[node.e_offset+i];
            assert(e.pt != REF_NULL);
            Node &ch = node_info[e.pt];
            if (ch.S == CERTAIN) {
                if (ch.Q < 0) { 
                    // win cut
                    std::swap(edge_info[node.e_offset+i], edge_info[node.e_offset]);
                    std::swap(edge_act[node.e_offset+i], edge_act[node.e_offset]);
                    node.e_beg = node.e_end;
                    break;
                }
                else {
                    // normal cut
                    std::swap(edge_info[node.e_offset+i], edge_info[node.e_offset+node.e_beg]);
                    std::swap(edge_act[node.e_offset+i], edge_act[node.e_offset+node.e_beg]);
                    if (node.e_beg > 0) {
                        Edge & best_e = edge_info[node.e_offset];
                        Node & best_ch = node_info[best_e.pt];
                        if (ch.Q < best_ch.Q) {
                            std::swap(edge_info[node.e_offset], edge_info[node.e_offset+node.e_beg]);
                            std::swap(edge_act[node.e_offset], edge_act[node.e_offset+node.e_beg]);
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
            Edge & best_edge = edge_info[node.e_offset];
            Node & best_ch = node_info[best_edge.pt];
            node.Q = -best_ch.Q;
            return 0;
        } else {
            if (sum_edge_N == 0) {
                node.U = - sum_edge_PxQ / sum_edge_P;
            }
            // double u = - sum_edge_PxQ / sum_edge_P;
            Edge & best_edge = edge_info[node.e_offset];
            Node & best_ch = node_info[best_edge.pt];
            if (node.e_beg > 0)
                node.U = std::max(node.U, (float)-best_ch.Q);
            // if (node.e_beg > 0)
            //     u = std::max(u, -best_ch.Q);
            node.Q = (-sum_edge_NxQ + node.U) / (sum_edge_N + 1);
            if (node.e_beg > 0)
                node.Q = std::max(node.Q, -best_ch.Q);
            return sum_edge_N + 1;
        }
    }
    
    uint32_t get_best_edge (uint32_t nid) const {
        double max_gain = -2;
        uint32_t ret = REF_NULL;
        const Node & node = node_info[nid];
        uint32_t N = 1;
        for (Move i = node.e_beg; i < node.e_end; i++) {
            const Edge &e = edge_info[node.e_offset+i];
            N += e.N;
        }
        double param = std::sqrt(N) * c_puct / node.e_end;
        for (Move i = node.e_beg; i < node.e_end; i++) {
            const Edge &e = edge_info[node.e_offset+i];
            assert(e.pt != REF_NULL);
            const Node &ch = node_info[e.pt];
            double gain; 
            if (ch.S == CERTAIN) {
                gain = -ch.Q;
            } else {
                gain = param / (e.N+1) + (ch.e_end == 0 ? (node.U-ch.Q)/2 : -ch.Q);
            }
            if (max_gain < gain) {
                max_gain = gain;
                ret = node.e_offset + i;
            }
        }
        return ret;
    }

    uint32_t search(uint32_t nid) {
        static Move all[16];
        // FIXME reference invalid when increase vector size
        Node & node = node_info[nid];
        if (node.S == CERTAIN) return 0;
        if (node.e_end > 0) {
            uint32_t eid = get_best_edge(nid);
            uint32_t ch_n = search(edge_info[eid].pt);
            if (edge_info[eid].N < ch_n) {
                edge_info[eid].N += 1;
            }
        } else {
            // expand
            Plate pl = Plate::get_plate(node_key[nid]);
            ExInfo exi = pl.build();
            double v = pl.get_move(exi, all, node.e_end);
            if (v != -2) {
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
                    tmp_pl.step_exi(all[i], tmp_exi);
                    Plate::Key ch_key = tmp_pl.get_key();
                    uint32_t ch_slot = find_slot(ch_key);
                    if (hash_table[ch_slot] == REF_NULL) {
                        hash_table[ch_slot] = node_info.size();
                        node_info.emplace_back();
                        node_key.push_back(ch_key);
                        Move mv = MV_NULL;
                        std::tie(mv, node_info.back().Q) = tmp_pl.forward_check(tmp_exi);
                        if (mv != MV_NULL) {
                            node_info.back().S = CERTAIN;
                            node_info.back().e_offset = mv;
                        }
                    }
                    edge_info.back().pt = hash_table[ch_slot];
                }
            }
        }
        int ret = update(nid);
        return ret;
    }
};


#endif