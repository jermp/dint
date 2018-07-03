#include <iostream>
#include <fstream>
#include <algorithm>

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/filesystem.hpp>

#include <boost/iostreams/device/mapped_file.hpp>
#include <sys/mman.h>

#include "codecs.hpp"
#include "util.hpp"
#include "dictionary.hpp"

using namespace ds2i;

template<typename Decoder>
void decode(std::string const& type,
            char const* encoded_data_filename,
            char const* dictionary_filename)
{
    boost::iostreams::mapped_file_source file;
    file.open(encoded_data_filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening index file");
    }

    uint8_t const* begin = (uint8_t const*) file.data();
    uint64_t size = file.size() / sizeof(uint8_t);
    uint8_t const* end = begin + size;
    auto ret = posix_madvise((void*) begin, size, POSIX_MADV_SEQUENTIAL);
    if (ret) {
        logger() << "Error calling madvice: " << errno << std::endl;
    }

    dictionary_type dict;
    uint64_t dict_size = 0;
    if (dictionary_filename) {
        typename dictionary_type::builder builder;
        std::ifstream dictionary_file(dictionary_filename);
        builder.load(dictionary_file);
        dict_size = builder.size();
        builder.build(dict);
    }

    const static uint64_t MAX_SIZE = 50000000; // ensure enough space for the largest sequence
    std::vector<uint32_t> decoded;
    decoded.resize(MAX_SIZE, 0);

    logger() << "decoding..." << std::endl;

    uint64_t num_decoded_ints = 0;
    uint64_t num_decoded_lists = 0;
    std::vector<double> timings;

    dint_statistics stats(dictionary_type::num_entries);

    // bool emit_selectors = false;

    while (begin != end) {
        uint32_t n, universe;
        begin = header::read(begin, &n, &universe);
        // if (n < 8)
        //     std::cout << "n = " << n << "; universe = " << universe << std::endl;
        auto start = clock_type::now();
        begin = Decoder::decode(begin, decoded.data(), universe, n, &dict
                                // , stats
                                // , emit_selectors
                                );
        auto finish = clock_type::now();
        std::chrono::duration<double> elapsed = finish - start;
        timings.push_back(elapsed.count());
        num_decoded_ints += n;
        ++num_decoded_lists;
    }

    double tot_elapsed = std::accumulate(timings.begin(), timings.end(), double(0.0));
    double ns_x_int = tot_elapsed * 1000000000 / num_decoded_ints;
    uint64_t ints_x_sec = uint64_t(1 / ns_x_int * 1000000000);

    logger() << "elapsed time " << tot_elapsed << " [sec]" << std::endl;
    logger() << ns_x_int << " [ns] x int" << std::endl;
    logger() << ints_x_sec << " ints x [sec]" << std::endl;

    // stats to std output
    // std::cout << "{";
    // std::cout << "\"filename\": \"" << encoded_data_filename << "\", ";
    // std::cout << "\"num_sequences\": \"" << num_decoded_lists << "\", ";
    // std::cout << "\"num_integers\": \"" << num_decoded_ints << "\", ";
    // std::cout << "\"type\": \"" << type << "\", ";
    // std::cout << "\"tot_elapsed_time\": \"" << tot_elapsed << "\", ";
    // std::cout << "\"ns_x_int\": \"" << ns_x_int << "\", ";
    // std::cout << "\"ints_x_sec\": \"" << ints_x_sec << "\"";
    // std::cout << "}" << std::endl;

    // uint64_t total_codewords = 0;
    // uint64_t total_decoded_ints = 0;
    // for (int i = 0; i < stats.codewords_distr.size(); ++i) {
    //     total_codewords += stats.codewords_distr[i];
    //     total_decoded_ints += stats.ints_distr[i];
    // }

    // std::cout << "total_codewords " << total_codewords << std::endl;
    // std::cout << "total_decoded_ints " << total_decoded_ints << std::endl;

    // for (int i = 0; i < stats.codewords_distr.size(); ++i) {
    //     if (i == 0) std::cout << "freq:\n";
    //     else if (i == stats.codewords_distr.size() - 1) std::cout << "rare:\n";
    //     else std::cout << (uint32_t(1) << (i - 1)) << ":\n";
    //     std::cout << "\t codewords: " << stats.codewords_distr[i] * 100.0 / total_codewords << "%" << std::endl;
    //     std::cout << "\t integers: " << stats.ints_distr[i] * 100.0 / total_decoded_ints << "%" << std::endl;
    // }






    // NOTE: old stuff
    // std::cout << "{";
    // std::cout << "\"filename\": \"" << encoded_data_filename << "\", ";
    // std::cout << "\"num_sequences\": \"" << num_decoded_lists << "\", ";
    // std::cout << "\"num_integers\": \"" << total_decoded_ints << "\", ";
    // std::cout << "\"num_codewords\": \"" << total_codewords << "\", ";

    // std::cout << "\"ints_runs\": \"" << stats.ints[0] << "\", ";
    // std::cout << "\"ints_table\": \"" << stats.ints[1] << "\", ";
    // std::cout << "\"ints_exception\": \"" << stats.ints[2] << "\", ";

    // std::cout << "\"codewords_runs\": \"" << stats.codewords[0] << "\", ";
    // std::cout << "\"codewords_table\": \"" << stats.codewords[1] << "\", ";
    // std::cout << "\"codewords_exceptions\": \"" << stats.codewords[2] << "\", ";

    // std::cout << "\"codewords_16_ints\": \"" << stats.codewords_distr[0] << "\", ";
    // std::cout << "\"codewords_8_ints\": \"" << stats.codewords_distr[1] << "\", ";
    // std::cout << "\"codewords_4_ints\": \"" << stats.codewords_distr[2] << "\", ";
    // std::cout << "\"codewords_2_ints\": \"" << stats.codewords_distr[3] << "\", ";
    // std::cout << "\"codewords_1_ints\": \"" << stats.codewords_distr[4] << "\", ";
    // std::cout << "\"codewords_ex_ints\": \"" << stats.codewords_distr[5] << "\", ";

    // std::cout << "\"exceptions_1_bytes\": \"" << stats.exceptions[0] << "\", ";
    // std::cout << "\"exceptions_2_bytes\": \"" << stats.exceptions[1] << "\", ";
    // std::cout << "\"exceptions_3_bytes\": \"" << stats.exceptions[2] << "\", ";
    // std::cout << "\"exceptions_4_bytes\": \"" << stats.exceptions[3] << "\"";

    // std::cout << "}" << std::endl;

    // sorting in decreasing frequency order
    // for (uint32_t i = 0; i < constants::max_fractal_steps; ++i) {
    //     std::sort(stats.codewords_freqs[i].begin(),
    //               stats.codewords_freqs[i].end(),
    //               [](auto const& p_x, auto const& p_y) {
    //                     return p_x.second > p_y.second;
    //               });
    // }

    // (index, frequency)
    // std::vector<std::pair<uint32_t, uint64_t>> v;
    // v.reserve(stats.freqs.size());

    // uint32_t index = 0;
    // for (auto f: stats.freqs) {
    //     v.emplace_back(index, f);
    //     ++index;
    // }

    // std::sort(v.begin(), v.end(),
    //     [](auto const& f_x, auto const& f_y) {
    //         return f_x.second > f_y.second;
    //     });

    // logger() << "computing space usage" << std::endl;
    // index = 0;
    // uint64_t bits = 0;
    // uint64_t sum_freqs = 0;
    // for (auto const& p: v) {
    //     uint32_t frequency = p.second;
    //     uint32_t codeword_bits = floor_log2(index + 2);
    //     bits += codeword_bits * frequency;
    //     ++index;
    //     sum_freqs += frequency;
    // }

    // logger() << "total integers: " << num_decoded_ints << std::endl;
    // logger() << "total codewords: " << total_codewords << std::endl;
    // // logger() << "sum frequencies: " << sum_freqs << std::endl;
    // // logger() << "total bits: " << bits << std::endl;

    // double BPI_cw_fixed = total_codewords * 16.0 / num_decoded_ints;
    // double BPI_cw_variable = double(bits) / num_decoded_ints;
    // double BPI_cw_len = 4.0 * sum_freqs / num_decoded_ints;
    // double BPI_ex = stats.ints[2] * 32.0 / num_decoded_ints;

    // logger() << "BPI for variable-len. codewords: " << BPI_cw_variable << std::endl;
    // logger() << "BPI for codewords' lengths: " << BPI_cw_len << std::endl;
    // logger() << "BPI exceptions: " << BPI_ex << std::endl;

    // logger() << "total BPI for fixed-len. codewords: " << BPI_cw_fixed << std::endl;
    // logger() << "total BPI for variable-len. codewords: " << BPI_cw_variable + BPI_cw_len + BPI_ex << std::endl;

    // substitute the freq value in freqs with their rank to derive the codeword
    // index = 0;
    // for (auto const& p: v) {
    //     stats.freqs[p.first] = index;
    //     ++index;
    // }

    // uint64_t total_covered_codewords = 0;
    // for (uint32_t i = 0; i < constants::max_fractal_steps; ++i) {
    //     auto const& freqs = stats.codewords_freqs[i];
    //     uint64_t codewords = 0;
    //     for (uint32_t k = 0; k < constants::top_k; ++k) {
    //         codewords += freqs[k].second;
    //         // std::cout << "index: " << freqs[k].first << "; freq: " << freqs[k].second << std::endl;
    //     }
    //     total_covered_codewords += codewords;
    //     std::cout << "covering " << codewords * 100.0 / total_codewords // stats.codewords[1]
    //               << "% of codewords "
    //               << "with the " << constants::top_k << " most frequent codewords for "
    //               << "targets of length " << (uint32_t(1) << i) << std::endl;
    // }
    // std::cout << "total covered codewords: "
    //           << total_covered_codewords * 100.0 / total_codewords // stats.codewords[1]
    //           << "%" << std::endl;

    // for (uint32_t i = 0; i < constants::max_fractal_steps; ++i) {
    //     std::sort(stats.codewords_freqs[i].begin(),
    //               stats.codewords_freqs[i].begin() + constants::top_k,
    //               [](auto const& p_x, auto const& p_y) {
    //                     return p_x.first < p_y.first;
    //               });
    // }

    // for (uint32_t i = 0; i < constants::max_fractal_steps; ++i) {
    //     auto const& freqs = stats.codewords_freqs[i];
    //     for (uint32_t k = 0; k < constants::top_k; ++k) {
    //         std::cout << "index: " << freqs[k].first << "; freq: " << freqs[k].second << std::endl;
    //     }
    // }

    // begin = (uint8_t const*) file.data();
    // emit_selectors = true;
    // while (begin != end) {
    //     uint32_t n, universe;
    //     begin = header::read(begin, &n, &universe);
    //     begin = Decoder::decode(begin, decoded.data(), universe, n, &dict
    //                             , stats, emit_selectors
    //                             );
    // }

    // {
    //     std::vector<std::pair<uint32_t, uint64_t>> v;
    //     v.reserve(dictionary_type::num_entries);

    //     for (uint32_t i = dictionary_type::reserved;
    //                   i < dict_size; ++i)
    //     {
    //         v.emplace_back(i, stats.freqs[i]);
    //     }
    //     std::sort(v.begin(), v.end(),
    //         [](auto const& x, auto const& y) {
    //             return x.second > y.second;
    //         });
    //     for (auto const& p: v) {
    //         // std::cout << std::setw( 6) << p.first
    //         //           << std::setw(23) << "freq: " << p.second << "; ";
    //         // std::cout << "entry: ";
    //         // dict.print(p.first);
    //         std::cout << p.first << std::endl;
    //     }
    //     // std::cout << std::endl;
    // }

    // {
    //     std::cout << "decoded " << stats.exceptions_freqs.size() << " distinct exceptions" << std::endl;
    //     std::vector<std::pair<uint32_t, uint64_t>> v;
    //     v.reserve(stats.exceptions_freqs.size());
    //     for (auto const& p: stats.exceptions_freqs) {
    //         v.emplace_back(p.first, p.second);
    //     }
    //     std::sort(v.begin(), v.end(),
    //         [](auto const& x, auto const& y) {
    //             return x.second > y.second;
    //         }
    //     );

    //     uint32_t i = 0;
    //     for (auto const& p: v) {
    //         std::cout << std::setw( 6) << p.first
    //                   << std::setw(23) << "freq: " << p.second << "\n";
    //         ++i;
    //         if (i == 50) break;
    //     }
    //     std::cout << std::endl;
    // }

    file.close();
}

int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t<type> <encoded_data_filename> [--dict <dictionary_filename>]"
                  << std::endl;
        return 1;
    }

    using namespace ds2i;
    std::string type = argv[1];
    char const* encoded_data_filename = argv[2];
    char const* dictionary_filename = nullptr;

    std::string cmd(std::string(argv[0]) + " " + type + " " + std::string(encoded_data_filename));

    for (int i = 3; i < argc; ++i) {
        if (argv[i] == std::string("--dict")) {
            ++i;
            dictionary_filename = argv[i];
            cmd += " --dict " + std::string(dictionary_filename);
        } else {
            throw std::runtime_error("unknown parameter");
        }
    }

    logger() << cmd << std::endl;

    // decode<dint>(type, encoded_data_filename, dictionary_filename);

    if (false) {
#define LOOP_BODY(R, DATA, T)                                          \
        } else if (type == BOOST_PP_STRINGIZE(T)) {                    \
            decode<BOOST_PP_CAT(T, )>                                  \
                (type, encoded_data_filename, dictionary_filename);    \
            /**/

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, CODECS);
#undef LOOP_BODY
    } else {
        logger() << "ERROR: unknown type '"
                 << type << "'" << std::endl;
    }

    return 0;
}
