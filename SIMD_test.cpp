#include <iostream>
#include <vector>
#include <cmath>

#include <x86intrin.h>
#include <immintrin.h> // AVX2

#include "dint_configuration.hpp"
#include "dictionary.hpp"
#include "codecs.hpp"
#include "util.hpp"
#include "hash_utils.hpp"
#include "binary_collection.hpp"
#include "tables.hpp"

// NOTE: stolen from Daniel Lemire code
#define RDTSC_START(cycles)                                                   \
    do {                                                                      \
        unsigned cyc_high, cyc_low;                                           \
        __asm volatile(                                                       \
            "cpuid\n\t"                                                       \
            "rdtsc\n\t"                                                       \
            "mov %%edx, %0\n\t"                                               \
            "mov %%eax, %1\n\t"                                               \
            : "=r"(cyc_high), "=r"(cyc_low)::"%rax", "%rbx", "%rcx", "%rdx"); \
        (cycles) = ((uint64_t)cyc_high << 32) | cyc_low;                      \
    } while (0)

#define RDTSC_FINAL(cycles)                                                   \
    do {                                                                      \
        unsigned cyc_high, cyc_low;                                           \
        __asm volatile(                                                       \
            "rdtscp\n\t"                                                      \
            "mov %%edx, %0\n\t"                                               \
            "mov %%eax, %1\n\t"                                               \
            "cpuid\n\t"                                                       \
            : "=r"(cyc_high), "=r"(cyc_low)::"%rax", "%rbx", "%rcx", "%rdx"); \
        (cycles) = ((uint64_t)cyc_high << 32) | cyc_low;                      \
    } while (0)

const static uint64_t GiB = 1073741824;

using namespace ds2i;

// NO EXCEPTIONS AND NO RUNS
uint64_t encode(uint32_t const* in,
                uint32_t n,
                std::vector<uint8_t>& out,
                dictionary_type::builder* builder)
{
    uint32_t const* begin = in;
    uint32_t const* end = begin + n;

    uint64_t written = 0;
    while (begin < end) {
        uint32_t index = EXCEPTIONS;
        for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
            uint32_t sub_block_size = constants::target_sizes[s];
            uint32_t len = std::min<uint32_t>(sub_block_size, end - begin);
            index = builder->lookup(begin, len);
            if (index != dictionary_type::invalid_index) {
                auto ptr = reinterpret_cast<uint8_t const*>(&index);
                // out.insert(out.end(), ptr, ptr + 2); // b = 16
                out.insert(out.end(), ptr, ptr + 1); // b = 8
                begin += len;
                written += len;
                break;
            }
        }

        if (index == dictionary_type::invalid_index) {
            begin += 1;
        }
    }

    return written;
}

uint64_t encode(char const* collection_name,
                char const* dictionary_filename,
                uint64_t n, std::vector<uint8_t>& encoded,
                dictionary_type::builder& builder)
{
    binary_collection input(collection_name);

    auto it = input.begin();
    uint64_t num_processed_lists = 0;
    uint64_t num_total_ints = 0;

    uint64_t total_progress = n;
    bool docs = true;
    boost::filesystem::path collection_path(collection_name);
    if (collection_path.extension() == ".freqs") {
        docs = false;
        logger() << "encoding freqs..." << std::endl;
    } else if (collection_path.extension() == ".docs") {
        // skip first singleton sequence, containing num. of docs
        ++it;
        logger() << "encoding docs..." << std::endl;
    } else {
        throw std::runtime_error("unsupported file format");
    }

    std::vector<uint32_t> buf;

    for (; it != input.end(); ++it)
    {
        auto const& list = *it;
        uint32_t size = list.size();
        if (size > constants::min_size)
        {
            buf.reserve(size);
            uint32_t prev = docs ? -1 : 0;
            for (auto b = list.begin(); b != list.end(); ++b) {
                buf.push_back(*b - prev - 1);
                if (docs) {
                    prev = *b;
                }
            }
            assert(buf.size() == size);

            uint64_t written = encode(buf.data(), size, encoded, &builder);
            buf.clear();

            ++num_processed_lists;
            num_total_ints += written;

            if (num_total_ints >= n) break;
        }
    }

    std::cerr << std::endl;

    double GiB_space = encoded.size() * 1.0 / GiB;
    double bpi_space = encoded.size() * sizeof(encoded.front()) * 8.0 / num_total_ints;

    logger() << "encoded " << num_processed_lists << " lists" << std::endl;
    logger() << "encoded " << num_total_ints << " integers" << std::endl;
    logger() << GiB_space << " [GiB]" << std::endl;
    logger() << "bits x integer: " << bpi_space << std::endl;

    // if (output_filename) {
    //     logger() << "writing encoded data..." << std::endl;
    //     std::ofstream output_file(output_filename);
    //     output_file.write(reinterpret_cast<char const*>(output.data()),
    //                       output.size() * sizeof(output[0]));
    //     output_file.close();
    //     logger() << "DONE" << std::endl;
    // }

    return num_total_ints;
}

uint64_t decode_scalar(std::vector<uint8_t> const& encoded,
                       dictionary_type const& dictionary, uint32_t* out)
{
    // uint16_t const* ptr = reinterpret_cast<uint16_t const*>(encoded.data());
    // uint64_t _16bit_words = encoded.size() * sizeof(encoded.front()) * 8 / 16;
    // uint64_t decoded_ints = 0;
    // for (size_t w = 0; w != _16bit_words; ++w) {
    //     uint32_t index = *ptr;

    //     // uint32_t target_size = dictionary.copy(index, out);
    //     uint32_t target_size = 1;
    //     *out = *(dictionary.data() + index * (dictionary_type::max_entry_size + 1));

    //     decoded_ints += target_size;
    //     out += target_size;
    //     ++ptr;
    // }
    // return decoded_ints;

    uint8_t const* ptr = encoded.data();
    uint64_t _8bit_words = encoded.size();
    uint64_t decoded_ints = 0;
    for (size_t w = 0; w != _8bit_words; ++w) {
        uint32_t index = *ptr;

        // uint32_t target_size = dictionary.copy(index, out);
        uint32_t target_size = 1;
        *out = *(dictionary.data() + index * (dictionary_type::max_entry_size + 1));

        decoded_ints += target_size;
        out += target_size;
        ++ptr;
    }
    return decoded_ints;
}

uint64_t _m256_decode_SIMD1(std::vector<uint8_t> const& encoded,
                            dictionary_type const& dictionary, uint32_t* out)
{
    const static __m256i mask = _mm256_set1_epi32(65535);
    // const static __m256i words_offsets = _mm256_set1_epi32(dictionary_type::max_entry_size + 1);

    uint64_t _256bit_words = encoded.size() * sizeof(encoded.front()) * 8 / 256;
    __m256i const* in = reinterpret_cast<__m256i const*>(encoded.data());
    __m256i* pout = reinterpret_cast<__m256i*>(out);

    uint64_t num_total_decoded_ints = 0;
    for (uint64_t i = 0; i != _256bit_words; ++i)
    {
        __m256i win = _mm256_lddqu_si256(in);

        __m256i wout = _mm256_and_si256(mask, win);
        __m256i w = _mm256_i32gather_epi32(dictionary.data(), wout, 4);
        _mm256_storeu_si256(pout + 0, w);
        // wout = _mm256_mullo_epi32(_mm256_and_si256(mask, win), words_offsets);
        // num_total_decoded_ints += decode8();

        wout = _mm256_srli_epi32(win, 16);
        w = _mm256_i32gather_epi32(dictionary.data(), wout, 4);
        _mm256_storeu_si256(pout + 1, w);
        // wout = _mm256_mullo_epi32(_mm256_srli_epi32(win, 16), words_offsets);
        // num_total_decoded_ints += decode8();

        num_total_decoded_ints += 16;
        pout += 2;

        ++in;
    }

    return num_total_decoded_ints;
}

uint64_t _m512_decode_SIMD1(std::vector<uint8_t> const& encoded,
                            dictionary_type const& dictionary, uint32_t* out)
{
    const static __m512i mask = _mm512_set1_epi32(65535);

    // NOTE: try to multiply by a power of 2 (l = 4)
    // const static __m512i words_offsets = _mm512_set1_epi32(dictionary_type::max_entry_size);
    // const static __m512i words_offsets = _mm512_set1_epi32(dictionary_type::max_entry_size + 1);

    uint64_t _512bit_words = encoded.size() * sizeof(encoded.front()) * 8 / 512;
    __m512i const* in = reinterpret_cast<__m512i const*>(encoded.data());
    __m512i* pout = reinterpret_cast<__m512i*>(out);

    uint64_t num_total_decoded_ints = 0;
    for (uint64_t i = 0; i != _512bit_words; ++i)
    {
        __m512i win = _mm512_loadu_si512(in);

        __m512i wout = _mm512_and_si512(mask, win);
        // wout = _mm512_mullo_epi32(wout, words_offsets);
        __m512i w = _mm512_i32gather_epi32(wout, dictionary.data(), 4);
        _mm512_storeu_si512(pout + 0, w);

        // wout = _mm512_mullo_epi32(_mm512_srli_epi32(mask, win), words_offsets);
        // num_total_decoded_ints += decode8();

        wout = _mm512_srli_epi32(win, 16);
        // wout = _mm512_mullo_epi32(wout, words_offsets);
        w = _mm512_i32gather_epi32(wout, dictionary.data(), 4);
        _mm512_storeu_si512(pout + 1, w);

        // wout = _mm512_mullo_epi32(_mm512_srli_epi32(win, 16), words_offsets);
        // num_total_decoded_ints += decode8();

        num_total_decoded_ints += 32;
        pout += 2;

        ++in;
    }

    return num_total_decoded_ints;
}

uint64_t _m512_decode_SIMD1_8bit(std::vector<uint8_t> const& encoded,
                                 dictionary_type const& dictionary, uint32_t* out)
{
    const static __m512i mask = _mm512_set1_epi32(255);

    uint64_t _512bit_words = encoded.size() * 8 / 512;
    __m512i const* in = reinterpret_cast<__m512i const*>(encoded.data());
    __m512i* pout = reinterpret_cast<__m512i*>(out);

    uint64_t num_total_decoded_ints = 0;
    for (uint64_t i = 0; i != _512bit_words; ++i)
    {
        __m512i win = _mm512_loadu_si512(in);
        __m512i wout;
        __m512i w;

        wout = _mm512_and_si512(mask, win);
        wout = _mm512_slli_epi32(wout, 2);
        w = _mm512_i32gather_epi32(wout, dictionary.data(), 4);
        _mm512_storeu_si512(pout + 0, w);

        win = _mm512_srli_epi32(win, 8);
        wout = _mm512_and_si512(mask, win);
        wout = _mm512_slli_epi32(wout, 2);
        w = _mm512_i32gather_epi32(wout, dictionary.data(), 4);
        _mm512_storeu_si512(pout + 1, w);

        win = _mm512_srli_epi32(win, 8);
        wout = _mm512_and_si512(mask, win);
        wout = _mm512_slli_epi32(wout, 2);
        w = _mm512_i32gather_epi32(wout, dictionary.data(), 4);
        _mm512_storeu_si512(pout + 2, w);

        win = _mm512_srli_epi32(win, 8);
        wout = _mm512_and_si512(mask, win);
        wout = _mm512_slli_epi32(wout, 2);
        w = _mm512_i32gather_epi32(wout, dictionary.data(), 4);
        _mm512_storeu_si512(pout + 3, w);

        num_total_decoded_ints += 64;
        pout += 4;

        ++in;
    }

    return num_total_decoded_ints;
}

uint64_t decode_SIMD(std::vector<uint8_t> const& encoded,
                     dictionary_type const& dictionary, uint32_t* out)
{
    const static __m128i mask = _mm_set1_epi32(65535);
    const static __m128i words_offsets = _mm_set1_epi32(dictionary_type::max_entry_size + 1);
    const static __m128i lengths_offsets = _mm_set1_epi32(dictionary_type::max_entry_size);
    const static __m128i increments = _mm_set1_epi32(1);
    static const uint32_t lengths_to_coefficients[] = {0, 0, 1, 0, 2};
                                                    // 0  1  2  3  4

    __m128i win, wout;
    uint64_t _128bit_words = encoded.size() * sizeof(encoded.front()) * 8 / 128;
    __m128i const* in = reinterpret_cast<__m128i const*>(encoded.data());

    __m128i* pout = reinterpret_cast<__m128i*>(out);

    auto decode4 = [&]() // decode 4 codewords at a time
    {
        static uint32_t lengths[] = {0, 0, 0, 0};
        _mm_store_si128(reinterpret_cast<__m128i*>(lengths),
                        _mm_i32gather_epi32(dictionary.data(), _mm_add_epi32(wout, lengths_offsets), 4));
        uint32_t decoded_ints = lengths[0] + lengths[1] + lengths[2] + lengths[3];
        uint32_t index = lengths_to_coefficients[lengths[0]] * 27
                       + lengths_to_coefficients[lengths[1]] *  9
                       + lengths_to_coefficients[lengths[2]] *  3
                       + lengths_to_coefficients[lengths[3]] *  1;
        index *= dictionary_type::max_entry_size;

        for (uint32_t i = 0; i < 4; ++i) {
            __m128i words = _mm_i32gather_epi32(dictionary.data(), wout, 4);
            __m128i vindex = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ds2i::tables::indices[index + i]));
            _mm_i32scatter_epi32(out, vindex, words, 4);
            wout = _mm_add_epi32(wout, increments);
        }
        pout += 4;

        return decoded_ints;
    };

    uint64_t num_total_decoded_ints = 0;
    for (uint64_t w = 0; w != _128bit_words; ++w)
    {
        win = _mm_loadu_si128(in);

        wout = _mm_mullo_epi32(_mm_and_si128(mask, win), words_offsets);
        num_total_decoded_ints += decode4();

        wout = _mm_mullo_epi32(_mm_srli_epi32(win, 16), words_offsets);
        num_total_decoded_ints += decode4();

        ++in;
    }

    return num_total_decoded_ints;
}

int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t<collection_name> <dictionary_filename> <n> <scalar|simd>"
                  << std::endl;
        return 1;
    }

    char const* collection_name = argv[1];
    char const* dictionary_filename = argv[2];
    uint64_t n = std::stoull(argv[3]);  // number of postings to encode
    std::string type = argv[4];

    dictionary_type::builder builder;
    if (dictionary_filename) {
        std::ifstream dictionary_file(dictionary_filename);
        builder.load(dictionary_file);
        logger() << "preparing for encoding..." << std::endl;
        builder.prepare_for_encoding();
        // builder.print();
    }

    std::vector<uint8_t> encoded;
    uint64_t bytes = 1 * GiB;
    encoded.reserve(bytes);
    n = encode(collection_name, dictionary_filename, n, encoded, builder);

    dictionary_type dictionary;
    builder.build(dictionary);

    logger() << type << std::endl;

    std::vector<uint32_t> out(n);
    uint64_t decoded_ints = 0;
    uint64_t cycles_start, cycles_final, cycles_diff;
    const int runs = 10;

    if (type == std::string("scalar"))
    {
        // auto start = clock_type::now();
        RDTSC_START(cycles_start);
        for (int run = 0; run != runs; ++run) {
            decoded_ints += decode_scalar(encoded, dictionary, out.data());
        }
        RDTSC_FINAL(cycles_final);
        cycles_diff = cycles_final - cycles_start;

        logger() << "total cycles: " << cycles_diff << std::endl;
        logger() << "decoded ints: " << decoded_ints << std::endl;
        logger() << "avg. cycles per integer: "
                 << double(cycles_diff) / decoded_ints << std::endl;

        // auto end = clock_type::now();
        // std::chrono::nanoseconds elapsed = end - start;
        // logger() << "total nanoseconds: " << elapsed.count() << std::endl;
        // logger() << "decoded ints: " << decoded_ints << std::endl;
        // logger() << "avg. nanoseconds per integer: "
        //          << double(elapsed.count()) / decoded_ints << std::endl;
    }
    else if (type == std::string("simd"))
    {
        // auto start = clock_type::now();

        RDTSC_START(cycles_start);
        for (int run = 0; run != runs; ++run) {
            // decoded_ints += _m256_decode_SIMD1(encoded, dictionary, out.data());
            // decoded_ints += _m512_decode_SIMD1(encoded, dictionary, out.data());
            decoded_ints += _m512_decode_SIMD1_8bit(encoded, dictionary, out.data());
        }
        RDTSC_FINAL(cycles_final);
        cycles_diff = cycles_final - cycles_start;

        logger() << "total cycles: " << cycles_diff << std::endl;
        logger() << "decoded ints: " << decoded_ints << std::endl;
        logger() << "avg. cycles per integer: "
                 << double(cycles_diff) / decoded_ints << std::endl;

        // auto end = clock_type::now();
        // std::chrono::nanoseconds elapsed = end - start;
        // logger() << "total nanoseconds: " << elapsed.count() << std::endl;
        // logger() << "decoded ints: " << decoded_ints << std::endl;
        // logger() << "avg. nanoseconds per integer: "
        //          << double(elapsed.count()) / decoded_ints << std::endl;
    }
    else {
        std::cerr << "unsupported type" << std::endl;
        return 1;
    }

    return 0;
}

// void push_entry(std::vector<uint32_t>& dict, int size, int l, int& pos, int& j)
// {
//     for (int i = 0; i < size; ++i) {
//         dict[pos + i] = j + i;
//     }
//     // pos += l + 1;
//     // j += l + 1;
//     pos += l;
//     j += l;
//     // dict[pos - 1] = size;
// }

// int main()
// {
//     const int N = 1000;
//     std::vector<uint32_t> out(N, 0);


//     // NOTE: 256-bit version without scatter and with a uint64_t-word dictionary
//     // std::vector<uint16_t> vec(16, 0);
//     // for (int i = 0; i < 16; ++i) {
//     //     vec[i] = i;
//     // }

//     // std::vector<uint64_t> dict(16, 0);
//     // for (int i = 0; i < 16; ++i) {
//     //     dict[i] = std::pow(2, i);
//     // }

//     // __m256i win;
//     // __m256i wout;
//     // __m256i* pout = (__m256i *) out.data();
//     // const __m256i mask = _mm256_set1_epi32(65535);
//     // auto dict_ptr = reinterpret_cast<long long const*>(dict.data());
//     // win = _mm256_lddqu_si256((__m256i*) vec.data());
//     // wout = _mm256_and_si256(mask, win);
//     // __m256i w0 = _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8);
//     // __m256i w1 = _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8);

//     // wout = _mm256_srli_epi32(win, 16);
//     // __m256i w2 = _mm256_i32gather_epi64(dict_ptr, _mm256_castsi256_si128(wout), 8);
//     // __m256i w3 = _mm256_i32gather_epi64(dict_ptr, _mm256_extractf128_si256(wout, 1), 8);

//     // _mm256_storeu_si256(pout + 0, _mm256_xor_si256(w0, _mm256_slli_epi64(w2, 32)));
//     // _mm256_storeu_si256(pout + 1, _mm256_xor_si256(w1, _mm256_slli_epi64(w3, 32)));

//     // for (int i = 0; i < 64; ++i) {
//     //     std::cout << out[i] << " ";
//     // }
//     // std::cout << std::endl;




//     // NOTE: 128-bit version with scatter and with a uint32_t-word dictionary
//     // std::vector<uint16_t> vec(16, 0);
//     // for (int i = 0; i < 16; ++i) {
//     //     vec[i] = i;
//     // }

//     // std::vector<uint32_t> dict(16, 0);
//     // for (int i = 0; i < 16; ++i) {
//     //     dict[i] = std::pow(2, i);
//     // }

//     // __m128i win;
//     // __m128i wout;
//     // __m128i vindex;
//     // __m128i* pout = (__m128i *) out.data();

//     // const static __m128i mask = _mm_set1_epi32(65535);
//     // const static uint32_t vindex0[4] = {0, 2, 4, 6};
//     // const static uint32_t vindex1[4] = {1, 3, 5, 7};

//     // auto dict_ptr = reinterpret_cast<int const*>(dict.data());
//     // win = _mm_loadu_si128((__m128i*) vec.data());
//     // wout = _mm_and_si128(mask, win);

//     // __m128i w0 = _mm_i32gather_epi32(dict_ptr, wout, 4);
//     // vindex = _mm_loadu_si128((__m128i*) vindex0);
//     // _mm_i32scatter_epi32(pout + 0, vindex, w0, 4);

//     // wout = _mm_srli_epi32(win, 16);
//     // __m128i w1 = _mm_i32gather_epi32(dict_ptr, wout, 4);
//     // vindex = _mm_loadu_si128((__m128i*) vindex1);
//     // _mm_i32scatter_epi32(pout + 0, vindex, w1, 4);


//     const int l = 4;

//     // NOTE: write groups of 8 codewords at a time, in an interleaved way, i.e.:
//     // say, the codewords we want to write are just 0 1 2 3 4 5 6 7, then we
//     // write them as 0 4 1 5 2 6 3 7.
//     std::vector<uint16_t> vec(16, 0);
//     int j = 0;
//     for (int i = 0; i < 8; i += 2, ++j) {
//         vec[i] = j;
//     }
//     for (int i = 1; i < 8; i += 2, ++j) {
//         vec[i] = j;
//     }

//     for (int i = 8; i < 16; i += 2, ++j) {
//         vec[i] = j;
//     }
//     for (int i = 9; i < 16; i += 2, ++j) {
//         vec[i] = j;
//     }

//     // NOTE: just l, not l + 1
//     std::vector<uint32_t> dict(16 * l, 0);
//     int pos = 0;
//     j = 1;

//     // 4 4 4 4
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);

//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);

//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);

//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);
//     push_entry(dict, 4, l, pos, j);

//     // // 1 2 1 1
//     // push_entry(dict, 1, l, pos, j);
//     // push_entry(dict, 2, l, pos, j);
//     // push_entry(dict, 1, l, pos, j);
//     // push_entry(dict, 1, l, pos, j);

//     // // 2 2 2 2
//     // push_entry(dict, 2, l, pos, j);
//     // push_entry(dict, 2, l, pos, j);
//     // push_entry(dict, 2, l, pos, j);
//     // push_entry(dict, 2, l, pos, j);

//     // // 2 1 2 4
//     // push_entry(dict, 2, l, pos, j);
//     // push_entry(dict, 1, l, pos, j);
//     // push_entry(dict, 2, l, pos, j);
//     // push_entry(dict, 4, l, pos, j);

//     // print dict for debug
//     pos = 0;
//     for (int i = 0; i < 16; ++i, pos += l) {
//         std::cout // << dict[pos + l] << ":
//                   << "["
//                   << dict[pos + 0] << "|"
//                   << dict[pos + 1] << "|"
//                   << dict[pos + 2] << "|"
//                   << dict[pos + 3] << "]\n";
//     }

//     __m128i win;
//     __m128i wout;
//     __m128i vindex;
//     __m128i* pout = (__m128i *) out.data();

//     const static __m128i mask = _mm_set1_epi32(65535);
//     const static __m128i ints = _mm_set_epi32(8, 7, 6, 5);

//     auto dict_ptr = reinterpret_cast<int const*>(dict.data());
//     win = _mm_loadu_si128((__m128i*) vec.data());
//     wout = _mm_and_si128(mask, win);

//     _mm_storeu_si128(pout + 3, wout);
//     _mm_storeu_si128(pout + 5, _mm_slli_epi32(wout, 2));

//     wout = _mm_slli_epi32(wout, 2);
//     __m128i w0 = _mm_i32gather_epi32(dict_ptr, wout, 4);
//     _mm_storeu_si128(pout + 0, w0);


//     // wout = _mm_srli_epi32(win, 16);
//     // __m128i w1 = _mm_i32gather_epi32(dict_ptr, wout, 4);
//     // vindex = _mm_loadu_si128((__m128i*) vindex1);
//     // _mm_i32scatter_epi32(pout + 0, vindex, w1, 4);


//     // __m512i o = _mm512_broadcast_i32x4(ints);
//     __m512i o = _mm512_broadcastd_epi32(ints);
//     _mm512_storeu_si512(out.data(), o);

//     for (int i = 0; i < 64; ++i) {
//         std::cout << out[i] << " ";
//     }
//     std::cout << std::endl;

//     return 0;
// }
