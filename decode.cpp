#include <iostream>
#include <fstream>
#include <algorithm>

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <sys/mman.h>

#include <succinct/mapper.hpp>

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

    std::vector<large_dictionary_type> large_dicts(constants::num_selectors);
    // std::vector<small_dictionary_type> small_dicts(constants::num_selectors);

    // NOTE: single dictionary
    // dictionary_type dict;
    // uint64_t dict_size = 0;

    if (dictionary_filename) {
        // typename dictionary_type::builder builder;
        // std::ifstream dictionary_file(dictionary_filename);
        // builder.load(dictionary_file);
        // dict_size = builder.size();
        // logger() << "dictionary with "
        //          << dict_size << " entries" << std::endl;
        // builder.print_usage();
        // builder.build(dict);

        // NOTE: contexts
        std::string prefix(dictionary_filename);
        for (int s = 0; s != constants::num_selectors; ++s)
        {
            std::string large_dict_filename = prefix + "."
                + std::to_string(constants::selector_codes[s]) + ".large";
            typename large_dictionary_type::builder large_dict_builder;
            large_dict_builder.load_from_file(large_dict_filename);
            large_dict_builder.build(large_dicts[s]);

            // std::string small_dict_filename = prefix + "."
            //     + std::to_string(constants::selector_codes[s]) + ".small";
            // typename small_dictionary_type::builder small_dict_builder;
            // small_dict_builder.load_from_file(small_dict_filename);
            // small_dict_builder.build(small_dicts[s]);
        }
    }

    std::vector<uint32_t> decoded;
    decoded.resize(constants::max_size, 0);

    logger() << "decoding..." << std::endl;

    uint64_t num_decoded_ints = 0;
    uint64_t num_decoded_lists = 0;
    std::vector<double> timings;

    dint_statistics stats;

    while (begin != end) {
        uint32_t n, universe;
        begin = header::read(begin, &n, &universe);
        auto start = clock_type::now();
        begin = Decoder::decode(large_dicts,
                                // small_dicts,
                                begin, decoded.data(), universe, n
                                , &dict
                                // , stats
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

    logger() << "avg. # of decoded integers x codeword: " << double(stats.decoded_ints_from_dict) / stats.dict_codewords << std::endl;

    // stats to std output
    std::cout << "{";
    std::cout << "\"filename\": \"" << encoded_data_filename << "\", ";
    std::cout << "\"num_sequences\": \"" << num_decoded_lists << "\", ";
    std::cout << "\"num_integers\": \"" << num_decoded_ints << "\", ";
    std::cout << "\"type\": \"" << type << "\", ";
    std::cout << "\"tot_elapsed_time\": \"" << tot_elapsed << "\", ";
    std::cout << "\"ns_x_int\": \"" << ns_x_int << "\", ";
    std::cout << "\"ints_x_sec\": \"" << ints_x_sec << "\"";
    std::cout << "}" << std::endl;

    // binary entropy
    {
        uint64_t total_exceptions = 0;

        uint64_t total_codewords = 0;
        for (uint64_t i = 0; i != constants::num_entries; ++i) {
            total_codewords += stats.occs[i];
        }

        for (auto const& pair: stats.exceptions) {
            total_exceptions += pair.second;
        }

        logger() << "total_exceptions " << total_exceptions << std::endl;

        double exceptions_entropy = 0.0;
        for (auto const& pair: stats.exceptions) {
            double p = double(pair.second) / total_exceptions;
            // std::cout << "p " << p << "; p * std::log2(1.0 / p) = " << p * std::log2(1.0 / p) << std::endl;
            exceptions_entropy += p * std::log2(1.0 / p);
            // std::cout << "entropy " << exceptions_entropy << std::endl;
        }

        // std::cout << "total_codewords " << total_codewords << std::endl;
        double entropy = 0.0;
        for (uint64_t i = 0; i != constants::num_entries; ++i) {
            double x = stats.occs[i];
            if (x != 0) {
                double p = x / total_codewords;
                entropy += p * std::log2(1.0 / p);
            }
        }

        logger() << "total_ints " << stats.total_ints << std::endl;
        logger() << "binary entropy of the source: " << (entropy * total_codewords + exceptions_entropy * total_exceptions) / 8.0 / constants::GiB << " [GiB]" << std::endl;
        logger() << "binary entropy x codeword: " << entropy << " bits" << std::endl;
        logger() << "binary entropy x exception: " << exceptions_entropy << " bits" << std::endl;
        logger() << "binary entropy x integer: " << (entropy * total_codewords + exceptions_entropy * total_exceptions) / stats.total_ints << " bits" << std::endl;
    }

    uint64_t total_codewords = 0;
    uint64_t total_decoded_ints = 0;
    for (uint64_t i = 0; i < stats.codewords_distr.size(); ++i) {
        total_codewords += stats.codewords_distr[i];
        total_decoded_ints += stats.ints_distr[i];
    }

    std::cout << "total_codewords " << total_codewords << std::endl;
    std::cout << "total_decoded_ints " << total_decoded_ints << std::endl;

    for (uint64_t i = 0; i < stats.codewords_distr.size(); ++i) {
        if (i == 0) std::cout << "freq:\n";
        else if (i == stats.codewords_distr.size() - 1) std::cout << "rare:\n";
        else std::cout << (uint32_t(1) << (i - 1)) << ":\n";
        std::cout << "\t codewords: " << stats.codewords_distr[i] * 100.0 / total_codewords << "%" << std::endl;
        std::cout << "\t integers: " << stats.ints_distr[i] * 100.0 / total_decoded_ints << "%" << std::endl;
    }

    file.close();
}

void decode_pef(char const* encoded_data_filename, bool freqs)
{
    succinct::bit_vector bv;
    boost::iostreams::mapped_file_source m(encoded_data_filename);
    succinct::mapper::map(bv, m);
    uint64_t num_bits = bv.size();

    std::vector<uint32_t> decoded;
    decoded.resize(constants::max_size, 0);
    std::vector<double> timings;

    uint64_t num_decoded_ints = 0;
    uint64_t num_decoded_lists = 0;
    uint64_t offset = 0;
    uint64_t universe, n;

    logger() << "decoding..." << std::endl;
    // std::cout << "num_bits " << num_bits << std::endl;
    while (offset < num_bits) {

        uint64_t next_offset = bv.get_bits(offset, 64);
        offset += 64;
        universe = bv.get_bits(offset, 32);
        // std::cout << "universe " << universe << "; n ";
        offset += 32;
        n = bv.get_bits(offset, 32);
        // std::cout << n << std::endl;
        offset += 32;

        auto start = clock_type::now();
        pef::decode(
            bv, decoded.data(), offset, universe, n, freqs
        );
        auto finish = clock_type::now();
        std::chrono::duration<double> elapsed = finish - start;
        timings.push_back(elapsed.count());
        num_decoded_ints += n;
        ++num_decoded_lists;
        offset = next_offset;

        // std::cout << num_decoded_lists << ": offset " << offset << "/" << num_bits << std::endl;
    }

    double tot_elapsed = std::accumulate(timings.begin(), timings.end(), double(0.0));
    double ns_x_int = tot_elapsed * 1000000000 / num_decoded_ints;
    uint64_t ints_x_sec = uint64_t(1 / ns_x_int * 1000000000);

    logger() << "elapsed time " << tot_elapsed << " [sec]" << std::endl;
    logger() << ns_x_int << " [ns] x int" << std::endl;
    logger() << ints_x_sec << " ints x [sec]" << std::endl;
}

int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t<type> <encoded_data_filename> [--dict <dictionary_filename>] [--freqs]"
                  << std::endl;
        return 1;
    }

    using namespace ds2i;
    std::string type = argv[1];
    char const* encoded_data_filename = argv[2];
    char const* dictionary_filename = nullptr;
    bool freqs = false;

    std::string cmd(std::string(argv[0]) + " " + type + " " + std::string(encoded_data_filename));

    for (int i = 3; i < argc; ++i) {
        if (argv[i] == std::string("--dict")) {
            ++i;
            dictionary_filename = argv[i];
            cmd += " --dict " + std::string(dictionary_filename);
        }
        else
        if (argv[i] == std::string("--freqs")) {
            freqs = true;
            ++i;
            cmd += " --freqs";
        }
        else {
            throw std::runtime_error("unknown parameter");
        }
    }

    logger() << cmd << std::endl;

    // TODO: refactor this later
    if (type == std::string("greedy_dint")) {
        decode<greedy_dint>(type, encoded_data_filename, dictionary_filename);
    }

    // if (type == std::string("opt_dint")) {
    //     decode<opt_dint>(type, encoded_data_filename, dictionary_filename);
    // }

    if (type == std::string("pef")) {
        decode_pef(encoded_data_filename, freqs);
    } else {

//     if (false) {
// #define LOOP_BODY(R, DATA, T)                                          \
//         } else if (type == BOOST_PP_STRINGIZE(T)) {                    \
//             decode<BOOST_PP_CAT(T, )>                                  \
//                 (type, encoded_data_filename, dictionary_filename);    \
//             /**/

//         BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, CODECS);
// #undef LOOP_BODY
//     } else {
//         logger() << "ERROR: unknown type '"
//                  << type << "'" << std::endl;
//     }
    }

    return 0;
}
