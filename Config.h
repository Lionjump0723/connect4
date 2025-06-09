#ifndef CONFIG_H_
#define CONFIG_H_
#include <cstdint>

const uint32_t REF_NULL = UINT32_MAX;

typedef uint8_t Move;

const Move MV_NULL = 0x7fu;

typedef uint8_t Status;

const Status STALE = 3;
const Status CERTAIN = 4; 
const double c_puct = 2.3;

const double Rdelta = 1.0/8192;
#endif