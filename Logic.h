#ifndef LOGIC_H_
#define LOGIC_H_

#include <tuple>
#include <cstdint>
#include <utility>
#include <functional>
#include "Config.h"

// Bitboard: 3 x uint64, 16 bits per column (up to 4 columns per word).
// Cell (x, y) maps to bit (x + (y&3)<<4) in arr[y>>2].
struct BitBoard {
    uint64_t arr[3];

    bool operator==(const BitBoard& other) const;
    void set(int x, int y);
    void unset(int x, int y);
    bool test(int x, int y) const;
    bool test_col(int y) const;

    BitBoard operator&(const BitBoard& other) const;
    BitBoard operator|(const BitBoard& other) const;
    BitBoard operator^(const BitBoard& other) const;
    void operator&=(const BitBoard& other);
    void operator|=(const BitBoard& other);
    void operator^=(const BitBoard& other);

    // Gravity: lowest empty cell per column becomes playable (011..1 + 1 trick)
    BitBoard getTop() const;
    BitBoard operator>>(int val) const;
    BitBoard operator<<(int val) const;

    // Open-three threat patterns (horizontal / diagonal)
    BitBoard h3Win() const;
    // Three stacked vertically in one column
    BitBoard v1Win() const;
    int count() const;
};

// Cached threat and drop masks for rollout / get_move
struct ExInfo {
    BitBoard my_h3, en_h3;  // horizontal/diag open threes
    BitBoard my_v1, en_v1;  // vertical triples
    BitBoard utop, ftop;    // 2nd / 3rd playable cell per column
};

struct Plate {
    static int M, N, noX, noY;
    static BitBoard VALID_POINT, NO_POINT;
    static uint32_t random_z_hash[192][2];
    static bool sflag;
    static Plate current_plate;

    BitBoard my_point;  // our stones (xor with bo_point for opponent)
    BitBoard bo_point;  // all stones on board
    BitBoard atop;      // first playable cell per column
    uint32_t h;         // Zobrist hash
    Move vacant;
    uint8_t role;

    typedef std::pair<uint32_t, std::pair<BitBoard, BitBoard>> Key;

    Key get_key() const;
    static Plate get_plate(const Key& key);
    ExInfo build() const;
    static void init(int _M, int _N, int _noX, int _noY, const int* _board);
    void print() const;
    void step(int y);
    void step_exi(int y, ExInfo& exi);

    // Returns kMoveUnknown (-2) if undecided; else value in [-1, 1] with Rdelta tie-break.
    // Fills data[] with legal columns; count = number of moves.
    double get_move(const ExInfo& exi, Move* data, Move& count) const;

    // Rollout from this position. MV_NULL move => unknown leaf value.
    std::tuple<Move, double> forward_check(const ExInfo& exi) const;
};

template <>
struct std::hash<Plate::Key> {
    size_t operator()(const Plate::Key& k) const { return k.first; }
};

#endif
