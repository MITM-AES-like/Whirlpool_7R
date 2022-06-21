#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdint.h>
#include <immintrin.h>
#include <array>
#include <tuple>
#include <vector>
#include <map>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <csetjmp>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "omp.h"

using namespace std;

// Cipher Parameters
#define CELL_InBits (8U)
#define CELL_Range (1 << (CELL_InBits))
#define SB_n (CELL_InBits)
#define SB_a (CELL_Range)
#define MP_Mode (0)


#define ROW_N (4)
#define COL_N (4)
#define STATE_N (ROW_N * COL_N)
#define STATE_InBytes (ROW_N * COL_N)
#define STATE_InBits (ROW_N * COL_N * SB_n)

// Attack configurations
#ifndef DoFF
#define DoFF (28)
#endif

#if (DoFF == 28)
#define DoFB (CELL_InBits>>1)
#else
#define DoFB (CELL_InBits)
#endif

const __m128i one_128 = _mm_set_epi64x(0xffffffffffffffffULL, 0xffffffffffffffffULL);
const __m128i zero_128 = _mm_setzero_si128();

#ifndef State_to_Parameter_Order
// 00, 04, 08, 12,
// 01, 05, 09, 13,
// 02, 06, 10, 14,
// 03, 07, 11, 15,
#define State_to_Parameter_Order( \
    x00, x04, x08, x12, \
    x01, x05, x09, x13, \
    x02, x06, x10, x14, \
    x03, x07, x11, x15) \
    x15, x14, x13, x12, x11, x10, x09, x08, x07, x06, x05, x04, x03, x02, x01, x00
#endif

#ifndef _
#define _ (0x00)
#endif

#ifndef PARTIAL_InBits
#define PARTIAL_InBits (40)
#endif
#define PARTIAL_InBytes (PARTIAL_InBits >> 3)

#if (PARTIAL_InBits == 64)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF,   _ , 0xFF,   _ ,
                                                        0xFF,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL_InBits == 56)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF,   _ , 0xFF,   _ ,
                                                        0xFF,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL_InBits == 48)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF,   _ , 0xFF,   _ ,
                                                          _ ,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL_InBits == 44)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0x0F,   _ , 0xFF,   _ ,
                                                          _ ,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL_InBits == 40)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                          _ ,   _ , 0xFF,   _ ,
                                                          _ ,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL_InBits == 36)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0x0F, 0xFF,   _ ,   _ ,
                                                          _ ,   _ , 0xFF,   _ ,
                                                          _ ,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL_InBits == 32)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                          _ , 0xFF,   _ ,   _ ,
                                                          _ ,   _ , 0xFF,   _ ,
                                                          _ ,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL_InBits == 28)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                          _ , 0xFF,   _ ,   _ ,
                                                          _ ,   _ , 0xFF,   _ ,
                                                          _ ,   _ ,   _ , 0x0F
                                                        ));
#else //(PARTIAL_InBits == 32)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                          _ , 0xFF,   _ ,   _ ,
                                                          _ ,   _ , 0xFF,   _ ,
                                                          _ ,   _ ,   _ , 0xFF
                                                        ));
#endif

const __m128i GUESSED_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                          _ ,   _ ,   _ ,   _ ,
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF,   _ ,   _ ,   _ 
                                                        ));


#define PRINT_COUNT 1
#define DEBUG 0

#ifndef omp_nb_threads
#define omp_nb_threads (4)
#endif

typedef uint8_t cell_t;
typedef uint16_t cell2_t;
typedef uint32_t cell3_t;
typedef uint32_t cell4_t;
typedef array<cell_t, ROW_N> column_t;
typedef array<column_t, COL_N> matrix_t;
typedef __m128i state_t;

template<typename T>
inline T operator ^ (const T& lhs, const T& rhs)
{
  T tmp;
  for (size_t i = 0; i < lhs.size(); i++)
  {
    tmp[i] = lhs[i] ^ rhs[i];
  }
  return tmp;
}

#define EQU(w, u)  (_mm_testz_si128(one_128, _mm_xor_si128(w, u)) == 1)

inline void printstate(state_t &state)
{
    cell_t * state_pt = (cell_t *)&state;

    cout << hex << setfill('0');
    for (size_t ri = 0; ri < ROW_N; ri++)
    {
        for (size_t ci = 0; ci < COL_N; ci++)
        {

            cout << setw(sizeof(cell_t) * 2) << state_pt[ci * ROW_N + ri] + '\0' << "  ";
        }
        cout << endl;
    }
    cout << endl;
    cout << dec << setfill(' ');
}