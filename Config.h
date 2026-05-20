#ifndef CONFIG_H_
#define CONFIG_H_
#include <cstdint>

// Define PRINT (for example, compile with -DPRINT) to enable debug logs.
#ifdef PRINT
#define PRINT_CERR std::cerr
#else
#define PRINT_CERR if (true) {} else std::cerr
#endif

const uint32_t REF_NULL = UINT32_MAX;

typedef uint8_t Move;

const Move MV_NULL = 0x7fu;

typedef uint8_t Status;

const Status STALE = 3;
const Status CERTAIN = 4;

// PUCT exploration constant
const double c_puct = 2.3;

// Tie-break on piece count; avoids uniform UCT and "giving up" when losing
const double Rdelta = 1.0 / 8192;

// Search time budget (seconds) per move
const double kSearchTimeSec = 2.7;

// Max graph nodes before stopping expansion
const uint32_t kMaxGraphNodes = 1u << 21;

// Linear-probing hash table capacity (power of 2)
const uint32_t kHashTableCapacity = 1u << 22;

// Initial vector reserves for graph storage
const uint32_t kNodeReserve = (1u << 14) + 10;
const uint32_t kEdgeReserve = 46875;

// get_move: unknown outcome
const double kMoveUnknown = -2.0;

#endif
