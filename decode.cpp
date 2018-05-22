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

    dictionary dict;
    if (dictionary_filename) {
        dictionary::builder builder;
        std::ifstream dictionary_file(dictionary_filename);
        builder.load(dictionary_file);
        builder.build(dict);
    }

    const static uint64_t MAX_SIZE = 50000000; // ensure enough space for the largest sequence
    std::vector<uint32_t> decoded;
    decoded.resize(MAX_SIZE, 0);

    logger() << "decoding..." << std::endl;

    uint64_t num_decoded_ints = 0;
    uint64_t num_decoded_lists = 0;
    std::vector<double> timings;

    dint_statistics stats;

    while (begin != end) {
        uint32_t n, universe;
        begin = header::read(begin, &n, &universe);
        auto start = clock_type::now();
        begin = Decoder::decode(begin, decoded.data(), universe, n, &dict
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

    uint64_t total_codewords = stats.codewords[0] +
                               stats.codewords[1] +
                               stats.codewords[2] ;
    uint64_t total_decoded_ints = stats.ints[0] +
                                  stats.ints[1] +
                                  stats.ints[2] ;
    std::cout << "{";
    std::cout << "\"filename\": \"" << encoded_data_filename << "\", ";
    std::cout << "\"num_sequences\": \"" << num_decoded_lists << "\", ";
    std::cout << "\"num_integers\": \"" << total_decoded_ints << "\", ";
    std::cout << "\"num_codewords\": \"" << total_codewords << "\", ";

    std::cout << "\"ints_runs\": \"" << stats.ints[0] << "\", ";
    std::cout << "\"ints_table\": \"" << stats.ints[1] << "\", ";
    std::cout << "\"ints_exception\": \"" << stats.ints[2] << "\", ";

    std::cout << "\"codewords_runs\": \"" << stats.codewords[0] << "\", ";
    std::cout << "\"codewords_table\": \"" << stats.codewords[1] << "\", ";
    std::cout << "\"codewords_exception\": \"" << stats.codewords[2] << "\", ";

    std::cout << "\"codewords_16_ints\": \"" << stats.codewords_distr[0] << "\", ";
    std::cout << "\"codewords_8_ints\": \"" << stats.codewords_distr[1] << "\", ";
    std::cout << "\"codewords_4_ints\": \"" << stats.codewords_distr[2] << "\", ";
    std::cout << "\"codewords_2_ints\": \"" << stats.codewords_distr[3] << "\", ";
    std::cout << "\"codewords_1_ints\": \"" << stats.codewords_distr[4] << "\"";

    std::cout << "}" << std::endl;

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

    std::string cmd(std::string(argv[0]) + ": " + type + " " + std::string(encoded_data_filename));

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

    decode<dint>(type, encoded_data_filename, dictionary_filename);

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

    return 0;
}
