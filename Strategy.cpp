#include <random>
#include <chrono>
#include "Strategy.h"
using namespace std;

int Plate::M = -1, 
    Plate::N = -1, 
    Plate::noX = -1, 
    Plate::noY = -1;

BitBoard Plate::VALID_POINT = {{0,0,0}};
BitBoard Plate::NO_POINT = {{0,0,0}};
uint32_t Plate::random_z_hash [192][2] = {};
bool Plate::sflag = true;
Plate Plate::current_plate;

GTable graph(1<<22);

extern "C" Point* getPoint(const int M, const int N, const int* top, const int* _board, 
	const int lastX, const int lastY, const int noX, const int noY
) {
    // std::cerr << "enter" << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    
    Plate::init(M, N, noX, noY, _board);
    // std::cerr << "plate_init_success" << std::endl;
    {
        Plate::Key root_key = Plate::current_plate.get_key();
        graph.root_id = graph.hash_table[graph.find_slot(root_key)];
        GTable new_graph(1<<22);
        if (graph.root_id != REF_NULL) {
            new_graph.copy_table(graph);
        } else {
            uint32_t slot = new_graph.find_slot(root_key);
            new_graph.root_id = 0;
            new_graph.hash_table[slot] = 0;
            new_graph.node_info.emplace_back();
            new_graph.node_key.push_back(root_key);
        }

        std::swap(graph, new_graph);
        std::cerr << "memory peak: " << (graph.memory_Byte() + new_graph.memory_Byte() + 0.0) / 1048576 << "MB" << std::endl;
    }
    std::cerr << "init_time: " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t1).count() << std::endl;
    std::cerr << "begin graph size: node " << graph.node_info.size() << ", edge " << graph.edge_info.size() << std::endl;
    
    while(
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t1).count() < 2.7
        && graph.node_info.size() <= (1<<21)
        && graph.search(graph.root_id)
    ) {
    }

    std::cerr << "end graph size: node " << graph.node_info.size() << ", edge " << graph.edge_info.size() << std::endl;
    int y = 0;
    Node & root_node = graph.node_info[graph.root_id];
    std::cerr << "root count: " << graph.update(graph.root_id) << std::endl;
    if (root_node.S == CERTAIN) {
        if (root_node.e_end > 0) {
            y = graph.edge_act[root_node.e_offset];
        } else {
            y = root_node.e_offset;
        }
        std::cerr << "certain!" << root_node.Q << " " << y << std::endl;
    } else {
        uint32_t max_n = 0;
        double ch_Q = -2;
        for (Move i = root_node.e_beg; i < root_node.e_end; i++) {
            Edge & e = graph.edge_info[root_node.e_offset+i];
            if (e.N >= max_n / 4 && -graph.node_info[e.pt].Q > ch_Q) {
                max_n = e.N;
                y = graph.edge_act[root_node.e_offset+i];
                ch_Q = -graph.node_info[e.pt].Q;
            }
        }
        if (root_node.e_beg > 0) {
            Edge & best_e = graph.edge_info[root_node.e_offset];
            double best_Q = -graph.node_info[best_e.pt].Q;
            if (ch_Q < best_Q) {
                y = graph.edge_act[root_node.e_offset];
                ch_Q = best_Q;
            }
        }
        std::cerr << ch_Q << " " << y << std::endl;
    }
    
    return new Point(top[y]-1, y);
}

extern "C" void clearPoint(Point* p) {
    delete p;
}