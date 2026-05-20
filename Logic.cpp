#include "Logic.h"
#include <random>
#include <iostream>
#include <cassert>
#include <algorithm>

// --- BoardConfig ---

void BoardConfig::reset() {
    M = -1;
    N = -1;
    noX = -1;
    noY = -1;
    VALID_POINT = {{0, 0, 0}};
    NO_POINT = {{0, 0, 0}};
    configured = false;
}

void BoardConfig::configure(int _M, int _N, int _noX, int _noY) {
    if (configured && M == _M && N == _N && noX == _noX && noY == _noY) {
        return;
    }
    M = _M;
    N = _N;
    noX = _noX;
    noY = _noY;
    VALID_POINT = {{0, 0, 0}};
    NO_POINT = {{0, 0, 0}};

    std::mt19937 my_rnd(42);
    for (int i = 0; i < 192; i++) {
        random_z_hash[i][0] = my_rnd();
        random_z_hash[i][1] = my_rnd();
    }
    for (int y = 0; y < N; y++) {
        VALID_POINT.arr[y >> 2] |= ((1ull << M) - 1) << ((y & 3) << 4);
    }
    if (noX >= 0 && noX < M && noY >= 0 && noY < N) {
        NO_POINT.set(M - 1 - noX, noY);
        VALID_POINT.unset(M - 1 - noX, noY);
    }
    configured = true;
}

Plate make_root_plate(const BoardConfig& cfg, const int* raw_board) {
    assert(cfg.configured);
    uint32_t h = 0;
    BitBoard my_point = {{0, 0, 0}}, bo_point = {{0, 0, 0}};
    for (int i = 0; i < cfg.M; i++) {
        for (int j = 0; j < cfg.N; j++) {
            int v = raw_board[i * cfg.N + j];
            if (v > 0) {
                bo_point.set(cfg.M - 1 - i, j);
                if (v - 1) {
                    my_point.set(cfg.M - 1 - i, j);
                }
                h ^= cfg.random_z_hash[cfg.M - 1 - i + (j << 4)][v - 1];
            }
        }
    }
    Plate ref;
    ref.role = 1;
    Plate ret = Plate::get_plate({h, {my_point, bo_point}}, cfg, ref);
    ret.role = 1;
    return ret;
}

// --- Plate ---

Plate::Key Plate::get_key() const {
    return {h, {my_point, bo_point}};
}

Plate Plate::get_plate(const Key& key, const BoardConfig& cfg, const Plate& root) {
    Plate ret;
    ret.my_point = key.second.first;
    ret.bo_point = key.second.second;
    ret.atop = (key.second.second | cfg.NO_POINT).getTop() & cfg.VALID_POINT;
    ret.h = key.first;
    ret.vacant = (uint8_t)(cfg.VALID_POINT ^ key.second.second).count();
    ret.role = ((root.bo_point.count() + key.second.second.count()) & 1)
                   ? (uint8_t)(1 - root.role)
                   : root.role;
    return ret;
}

ExInfo Plate::build(const BoardConfig& cfg) const {
    BitBoard base = bo_point | atop | cfg.NO_POINT;
    BitBoard utop = base.getTop() & cfg.VALID_POINT;
    base |= utop;
    BitBoard ftop = base.getTop() & cfg.VALID_POINT;
    return {
        my_point.h3Win(),
        (my_point ^ bo_point).h3Win(),
        my_point.v1Win(),
        (my_point ^ bo_point).v1Win(),
        utop,
        ftop,
    };
}

void Plate::print(const BoardConfig& cfg) const {
    assert(cfg.configured);
    for (int i = 0; i < cfg.M; i++) {
        for (int j = 0; j < cfg.N; j++) {
            if (bo_point.test(cfg.M - 1 - i, j)) {
                PRINT_CERR << (my_point.test(cfg.M - 1 - i, j) ? "2 " : "1 ");
            } else {
                PRINT_CERR << "0 ";
            }
        }
        PRINT_CERR << std::endl;
    }
}

void Plate::step(int y, const BoardConfig& cfg) {
    my_point ^= bo_point;
    const uint64_t mask = 0xffffull << ((y & 3) << 4);
    const int hhh = __builtin_ctzll(atop.arr[y >> 2] & mask);
    h ^= cfg.random_z_hash[((y >> 2) << 6) + hhh][role];
    bo_point.arr[y >> 2] |= atop.arr[y >> 2] & mask;
    atop = (bo_point | cfg.NO_POINT).getTop() & cfg.VALID_POINT;
    vacant -= 1;
    role ^= 1;
}

void Plate::step_exi(int y, ExInfo& exi, const BoardConfig& cfg) {
    step(y, cfg);
    exi.my_h3 = exi.en_h3;
    const BitBoard en_point = (my_point ^ bo_point);
    exi.en_h3 = en_point.h3Win();
    std::swap(exi.my_v1, exi.en_v1);
    exi.en_v1.arr[y >> 2] =
        (en_point.arr[y >> 2] << 3) & (en_point.arr[y >> 2] << 2) & (en_point.arr[y >> 2] << 1);
    BitBoard base = bo_point | atop | cfg.NO_POINT;
    exi.utop = base.getTop() & cfg.VALID_POINT;
    base |= exi.utop;
    exi.ftop = base.getTop() & cfg.VALID_POINT;
}

double Plate::get_move(const ExInfo& exi, Move* data, Move& count, const BoardConfig& cfg) const {
    const BitBoard win1 = (exi.my_h3 | exi.my_v1) & atop;
    const BitBoard lose1 = (exi.en_h3 | exi.en_v1) & atop;
    const BitBoard lose2 = (exi.en_h3 | exi.en_v1) & exi.utop;
    const BitBoard win2 = (exi.my_h3 | exi.my_v1) & exi.utop;
    const BitBoard win3 = (exi.my_h3 | exi.my_v1) & exi.ftop;

    Move force = MV_NULL;
    Move must = MV_NULL;
    Move good = MV_NULL;
    Move must_cnt = 0;
    count = 0;
    for (int y = 0; y < cfg.N; y++) {
        if (atop.test_col(y)) {
            force = y;
            if (win1.test_col(y)) {
                count = 1;
                data[0] = y;
                return 1 - (bo_point.count() + 1) * Rdelta;
            }
            if (lose1.test_col(y)) {
                must_cnt += 1;
                must = y;
                if (lose2.test_col(y)) {
                    must_cnt += 1;
                }
                data[count++] = y;
            } else {
                if (!lose2.test_col(y)) {
                    data[count++] = y;
                    if (win2.test_col(y) && win3.test_col(y)) {
                        good = y;
                    }
                }
            }
        }
    }
    if (vacant == 1) {
        count = 1;
        data[0] = force;
        return 0;
    } else if (must_cnt > 1) {
        count = 1;
        data[0] = must;
        return -(1 - (bo_point.count() + 2) * Rdelta);
    } else if (must_cnt == 1) {
        count = 1;
        data[0] = must;
        return kMoveUnknown;
    } else if (count == 0) {
        count = 1;
        data[0] = force;
        return -(1 - (bo_point.count() + 2) * Rdelta);
    } else if (good != MV_NULL) {
        count = 1;
        data[0] = good;
        return (1 - (bo_point.count() + 3) * Rdelta);
    } else {
        return kMoveUnknown;
    }
}

std::tuple<Move, double> Plate::forward_check(const ExInfo& _exi, const BoardConfig& cfg,
                                              std::mt19937& rnd) const {
    Move all[16];
    Move candidate[16];
    Move aggressive[16];

    bool my_flg = true, op_flg = true;
    Move first_mv = MV_NULL;
    ExInfo exi = _exi;
    Plate local_p = *this;
    while (true) {
        Move cnt = 0;
        Move ccnt = 0;
        Move acnt = 0;
        double v = local_p.get_move(exi, all, cnt, cfg);
        if (v != kMoveUnknown) {
            if (first_mv == MV_NULL) {
                first_mv = all[0];
            }
            if ((v < 0 || op_flg) && (v > 0 || my_flg)) {
                return std::make_tuple(first_mv, local_p.role == role ? v : -v);
            } else {
                return std::make_tuple(MV_NULL, local_p.role == role ? v : -v);
            }
        } else {
            assert(cnt > 0);
            BitBoard win2 = (exi.my_h3 | exi.my_v1) & exi.utop;
            for (int i = 0; i < cnt; i++) {
                if (win2.test_col(all[i])) {
                    aggressive[acnt++] = all[i];
                } else {
                    candidate[ccnt++] = all[i];
                }
            }
            Move chosen_y = MV_NULL;
            if (ccnt > 0) {
                if (ccnt > 1) {
                    std::swap(candidate[0], candidate[rnd() % ccnt]);
                }
                chosen_y = candidate[0];
            } else {
                if (acnt > 1) {
                    std::swap(aggressive[0], aggressive[rnd() % acnt]);
                }
                chosen_y = aggressive[0];
            }
            if (first_mv == MV_NULL) {
                first_mv = chosen_y;
            }
            local_p.step_exi(chosen_y, exi, cfg);
            my_flg &= (cnt == 1);
            std::swap(my_flg, op_flg);
        }
    }
}
