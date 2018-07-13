#include <iostream>
#include <vector>

#include <x86intrin.h>

int main()
{
    const int N = 1000;
    std::vector<uint32_t> out(N, -1);

    std::vector<uint16_t> vec(16, 0);
    for (int i = 0; i < 16; ++i) {
        vec[i] = i;
    }

    std::vector<uint64_t> dict(16, 0);
    for (int i = 0; i < 16; ++i) {
        dict[i] = i + 100;
    }
    // const __m128i mask = _mm_set1_epi32(65535);
    // __m128i indexes = _mm_loadu_si128((__m128i*) vec.data());
    // __m128i wout = _mm_and_si128(mask, indexes);
    // // indexes = _mm_mul_epi32(indexes, lengths); // obtain the proper offsets
    // // __m128i words = _mm256_i32gather_epi32(base, index, 4);
    // _mm_storeu_si128((__m128i*) out.data(), wout);
    // for (int i = 0; i < 32; ++i) {
    //     std::cout << out[i] << " ";
    // }
    // std::cout << std::endl;


    __m256i win;
    __m256i wout;
    __m256i* pout = (__m256i *) out.data();
    const __m256i mask = _mm256_set1_epi32(65535);
    auto dict_ptr = reinterpret_cast<long long const*>(dict.data());
    win = _mm256_lddqu_si256((__m256i*) vec.data());
    wout = _mm256_and_si256(mask, win);



    // _mm256_storeu_si256(pout + 5, wout);
    // _mm256_storeu_si256(pout + 0, _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8));
    // _mm256_storeu_si256(pout + 1, _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8));
    // wout = _mm256_srli_epi32(win, 16);
    // _mm256_storeu_si256(pout + 6, wout);
    // _mm256_storeu_si256(pout + 2, _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8));
    // _mm256_storeu_si256(pout + 3, _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8));




    __m256i w0 = _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8);
    __m256i w1 = _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8);

    wout = _mm256_srli_epi32(win, 16);
    __m256i w2 = _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8);
    __m256i w3 = _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8);

    _mm256_storeu_si256(pout + 0, _mm256_xor_si256(w0, _mm256_slli_epi64(w2, 32)));
    _mm256_storeu_si256(pout + 1, _mm256_xor_si256(w1, _mm256_slli_epi64(w3, 32)));

    for (int i = 0; i < 64; ++i) {
        std::cout << out[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
