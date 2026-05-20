#ifndef LOGIC_H_
#define LOGIC_H_

#include <tuple>
#include <cstdint>
#include <utility>
#include <functional>
#include <random>
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

inline bool BitBoard::operator==(const BitBoard& other) const {
    return arr[0] == other.arr[0] && arr[1] == other.arr[1] && arr[2] == other.arr[2];
}

inline void BitBoard::set(int x, int y) {
    arr[y >> 2] |= 1ull << (x + ((y & 3) << 4));
}

inline void BitBoard::unset(int x, int y) {
    arr[y >> 2] &= ~(1ull << (x + ((y & 3) << 4)));
}

inline bool BitBoard::test(int x, int y) const {
    return arr[y >> 2] & (1ull << (x + ((y & 3) << 4)));
}

inline bool BitBoard::test_col(int y) const {
    return arr[y >> 2] & (0xffffull << ((y & 3) << 4));
}

inline BitBoard BitBoard::operator&(const BitBoard& other) const {
    return {{arr[0] & other.arr[0], arr[1] & other.arr[1], arr[2] & other.arr[2]}};
}

inline BitBoard BitBoard::operator|(const BitBoard& other) const {
    return {{arr[0] | other.arr[0], arr[1] | other.arr[1], arr[2] | other.arr[2]}};
}

inline BitBoard BitBoard::operator^(const BitBoard& other) const {
    return {{arr[0] ^ other.arr[0], arr[1] ^ other.arr[1], arr[2] ^ other.arr[2]}};
}

inline void BitBoard::operator&=(const BitBoard& other) {
    arr[0] &= other.arr[0];
    arr[1] &= other.arr[1];
    arr[2] &= other.arr[2];
}

inline void BitBoard::operator|=(const BitBoard& other) {
    arr[0] |= other.arr[0];
    arr[1] |= other.arr[1];
    arr[2] |= other.arr[2];
}

inline void BitBoard::operator^=(const BitBoard& other) {
    arr[0] ^= other.arr[0];
    arr[1] ^= other.arr[1];
    arr[2] ^= other.arr[2];
}

inline BitBoard BitBoard::getTop() const {
    return {{arr[0] + 0x1000100010001ull, arr[1] + 0x1000100010001ull, arr[2] + 0x1000100010001ull}};
}

inline BitBoard BitBoard::operator>>(int val) const {
    return {{(arr[0] >> val) | (arr[1] << (64 - val)), (arr[1] >> val) | (arr[2] << (64 - val)), arr[2] >> val}};
}

inline BitBoard BitBoard::operator<<(int val) const {
    return {{(arr[0] << val), (arr[1] << val) | (arr[0] >> (64 - val)), (arr[2] << val) | (arr[1] >> (64 - val))}};
}

inline BitBoard BitBoard::h3Win() const {
    BitBoard ret = {{0, 0, 0}}, amask, omask;
#define help_h3Win(x, op0, op1)           \
    amask = (*this) op0 x;                \
    omask = (*this) op1(3 * x);           \
    omask |= amask;                       \
    amask &= (*this);                     \
    ret |= omask & (amask op1(2 * x));    \
    ret |= (omask op0(2 * x)) & (amask op0 x)

    help_h3Win(16, <<, >>);
    help_h3Win(15, <<, >>);
    help_h3Win(17, <<, >>);
#undef help_h3Win
    return ret;
}

inline BitBoard BitBoard::v1Win() const {
    return {{
        (arr[0] << 3) & (arr[0] << 2) & (arr[0] << 1),
        (arr[1] << 3) & (arr[1] << 2) & (arr[1] << 1),
        (arr[2] << 3) & (arr[2] << 2) & (arr[2] << 1),
    }};
}

inline int BitBoard::count() const {
    return __builtin_popcountll(arr[0]) + __builtin_popcountll(arr[1]) + __builtin_popcountll(arr[2]);
}

// Per-board geometry, masks, and Zobrist table (no global static state).
struct BoardConfig {
    int M = -1;
    int N = -1;
    int noX = -1;
    int noY = -1;
    BitBoard VALID_POINT = {{0, 0, 0}};
    BitBoard NO_POINT = {{0, 0, 0}};
    uint32_t random_z_hash[192][2] = {};
    bool configured = false;

    void configure(int M, int N, int noX, int noY);
    void reset();
};

// Cached threat and drop masks for rollout / get_move
struct ExInfo {
    BitBoard my_h3, en_h3;  // horizontal/diag open threes
    BitBoard my_v1, en_v1;  // vertical triples
    BitBoard utop, ftop;    // 2nd / 3rd playable cell per column
};

struct Plate {
    BitBoard my_point;  // our stones (xor with bo_point for opponent)
    BitBoard bo_point;  // all stones on board
    BitBoard atop;      // first playable cell per column
    uint32_t h;         // Zobrist hash
    Move vacant;
    uint8_t role;

    typedef std::pair<uint32_t, std::pair<BitBoard, BitBoard>> Key;

    Key get_key() const;
    static Plate get_plate(const Key& key, const BoardConfig& cfg, const Plate& root);
    ExInfo build(const BoardConfig& cfg) const;
    void print(const BoardConfig& cfg) const;
    void step(int y, const BoardConfig& cfg);
    void step_exi(int y, ExInfo& exi, const BoardConfig& cfg);

    // Returns kMoveUnknown (-2) if undecided; else value in [-1, 1] with Rdelta tie-break.
    double get_move(const ExInfo& exi, Move* data, Move& count, const BoardConfig& cfg) const;

    // Rollout from this position. MV_NULL move => unknown leaf value.
    std::tuple<Move, double> forward_check(const ExInfo& exi, const BoardConfig& cfg,
                                           std::mt19937& rnd) const;
};

Plate make_root_plate(const BoardConfig& cfg, const int* raw_board);

template <>
struct std::hash<Plate::Key> {
    size_t operator()(const Plate::Key& k) const { return k.first; }
};

#endif
