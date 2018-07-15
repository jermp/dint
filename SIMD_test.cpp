#include <iostream>
#include <vector>
#include <cmath>
#include <x86intrin.h>

#include "tables.hpp"

void push_entry(std::vector<uint32_t>& dict, int size, int l, int& pos, int& j)
{
    for (int i = 0; i < size; ++i) {
        dict[pos + i] = j + i;
    }
    pos += l + 1;
    j += l + 1;
    dict[pos - 1] = size;
}

int main()
{
    const int N = 1000;
    std::vector<uint32_t> out(N, 0);


    // NOTE: 256-bit version without scatter and with a uint64_t-word dictionary
    // std::vector<uint16_t> vec(16, 0);
    // for (int i = 0; i < 16; ++i) {
    //     vec[i] = i;
    // }

    // std::vector<uint64_t> dict(16, 0);
    // for (int i = 0; i < 16; ++i) {
    //     dict[i] = std::pow(2, i);
    // }

    // __m256i win;
    // __m256i wout;
    // __m256i* pout = (__m256i *) out.data();
    // const __m256i mask = _mm256_set1_epi32(65535);
    // auto dict_ptr = reinterpret_cast<long long const*>(dict.data());
    // win = _mm256_lddqu_si256((__m256i*) vec.data());
    // wout = _mm256_and_si256(mask, win);
    // __m256i w0 = _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8);
    // __m256i w1 = _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8);

    // wout = _mm256_srli_epi32(win, 16);
    // __m256i w2 = _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8);
    // __m256i w3 = _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8);

    // _mm256_storeu_si256(pout + 0, _mm256_xor_si256(w0, _mm256_slli_epi64(w2, 32)));
    // _mm256_storeu_si256(pout + 1, _mm256_xor_si256(w1, _mm256_slli_epi64(w3, 32)));

    // for (int i = 0; i < 64; ++i) {
    //     std::cout << out[i] << " ";
    // }
    // std::cout << std::endl;




    // NOTE: 128-bit version with scatter and with a uint32_t-word dictionary
    // std::vector<uint16_t> vec(16, 0);
    // for (int i = 0; i < 16; ++i) {
    //     vec[i] = i;
    // }

    // std::vector<uint32_t> dict(16, 0);
    // for (int i = 0; i < 16; ++i) {
    //     dict[i] = std::pow(2, i);
    // }

    // __m128i win;
    // __m128i wout;
    // __m128i vindex;
    // __m128i* pout = (__m128i *) out.data();

    // const static __m128i mask = _mm_set1_epi32(65535);
    // const static uint32_t vindex0[4] = {0, 2, 4, 6};
    // const static uint32_t vindex1[4] = {1, 3, 5, 7};

    // auto dict_ptr = reinterpret_cast<int const*>(dict.data());
    // win = _mm_loadu_si128((__m128i*) vec.data());
    // wout = _mm_and_si128(mask, win);

    // __m128i w0 = _mm_i32gather_epi32(dict_ptr, wout, 4);
    // vindex = _mm_loadu_si128((__m128i*) vindex0);
    // _mm_i32scatter_epi32(pout + 0, vindex, w0, 4);

    // wout = _mm_srli_epi32(win, 16);
    // __m128i w1 = _mm_i32gather_epi32(dict_ptr, wout, 4);
    // vindex = _mm_loadu_si128((__m128i*) vindex1);
    // _mm_i32scatter_epi32(pout + 0, vindex, w1, 4);


    // for (int i = 0; i < 64; ++i) {
    //     std::cout << out[i] << " ";
    // }
    // std::cout << std::endl;

    const int l = 4;

    // NOTE: write groups of 8 codewords at a time, in an interleaved way, i.e.:
    // say, the codewords we want to write are just 0 1 2 3 4 5 6 7, then we
    // write them as 0 4 1 5 2 6 3 7.
    std::vector<uint16_t> vec(16, 0);
    int j = 0;
    for (int i = 0; i < 8; i += 2, ++j) {
        vec[i] = j;
    }
    for (int i = 1; i < 8; i += 2, ++j) {
        vec[i] = j;
    }

    for (int i = 8; i < 16; i += 2, ++j) {
        vec[i] = j;
    }
    for (int i = 9; i < 16; i += 2, ++j) {
        vec[i] = j;
    }

    std::vector<uint32_t> dict(16 * (l + 1), 0);
    int pos = 0;
    j = 1;

    // 4 4 2 1
    push_entry(dict, 4, l, pos, j);
    push_entry(dict, 4, l, pos, j);
    push_entry(dict, 2, l, pos, j);
    push_entry(dict, 1, l, pos, j);

    // 1 2 1 1
    push_entry(dict, 1, l, pos, j);
    push_entry(dict, 2, l, pos, j);
    push_entry(dict, 1, l, pos, j);
    push_entry(dict, 1, l, pos, j);

    // 2 2 2 2
    push_entry(dict, 2, l, pos, j);
    push_entry(dict, 2, l, pos, j);
    push_entry(dict, 2, l, pos, j);
    push_entry(dict, 2, l, pos, j);

    // 2 1 2 4
    push_entry(dict, 2, l, pos, j);
    push_entry(dict, 1, l, pos, j);
    push_entry(dict, 2, l, pos, j);
    push_entry(dict, 4, l, pos, j);

    // print dict for debug
    pos = 0;
    for (int i = 0; i < 16; ++i, pos += l + 1) {
        std::cout << dict[pos + l] << ": ["
                  << dict[pos + 0] << "|"
                  << dict[pos + 1] << "|"
                  << dict[pos + 2] << "|"
                  << dict[pos + 3] << "]\n";
    }

    __m128i win;
    __m128i wout;

    uint32_t* pout = out.data();

    const static __m128i mask = _mm_set1_epi32(65535);
    const static __m128i words_offsets = _mm_set1_epi32(l + 1);
    const static __m128i lengths_offsets = _mm_set1_epi32(l);
    const static __m128i increments = _mm_set1_epi32(1);
    static const uint32_t lengths_to_coefficients[] = {0, 0, 1, 0, 2};
                                                    // 0  1  2  3  4

    // const static __m128i powers = _mm_set_epi32(27, 9, 3, 1);

    // const static int p3 = 3 * 3 * 3;
    // const static int p2 = 3 * 3;
    // const static int p1 = 3;
    // int index = (p3 * (3 - 1) + p2 * (3 - 1) + p1 * (2 - 1) + (1 - 1)) * l;
    // std::cout << "index " << index << std::endl;

    auto dict_ptr = reinterpret_cast<int const*>(dict.data());
    uint64_t _128bit_words = vec.size() * 16 / 128;
    __m128i const* in = reinterpret_cast<__m128i const*>(vec.data());

    auto decode4 = [&]() // decode 4 codewords at a time
    {
        static uint32_t lengths[] = {0, 0, 0, 0};
        _mm_store_si128(reinterpret_cast<__m128i*>(lengths),
                        _mm_i32gather_epi32(dict_ptr, _mm_add_epi32(wout, lengths_offsets), 4));
        int decoded_ints = lengths[0] + lengths[1] + lengths[2] + lengths[3];
        int index = lengths_to_coefficients[lengths[0]] * 27
                  + lengths_to_coefficients[lengths[1]] *  9
                  + lengths_to_coefficients[lengths[2]] *  3
                  + lengths_to_coefficients[lengths[3]] *  1;
        index *= l;
        // std::cout << "index " << index << std::endl;

        for (int i = 0; i < 4; ++i) {
            __m128i words = _mm_i32gather_epi32(dict_ptr, wout, 4);
            __m128i vindex = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ds2i::tables::indices[index + i]));
            _mm_i32scatter_epi32(pout, vindex, words, 4);
            wout = _mm_add_epi32(wout, increments);
        }
        pout += decoded_ints;
    };

    for (uint64_t w = 0; w != _128bit_words; ++w)
    {
        win = _mm_loadu_si128(in);

        wout = _mm_mullo_epi32(_mm_and_si128(mask, win), words_offsets);
        decode4();

        wout = _mm_mullo_epi32(_mm_srli_epi32(win, 16), words_offsets);
        decode4();

        ++in;
    }

    for (int i = 0; i < 64; ++i) {
        std::cout << out[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
