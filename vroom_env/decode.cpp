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
#include "dint_codecs.hpp"
#include "util.hpp"
#include "statistics.hpp"

using namespace ds2i;

template<typename Decoder>
void decode(std::string const& type,
            char const* encoded_data_filename)
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

    std::vector<uint32_t> decoded;
    decoded.resize(constants::max_size, 0);

    logger() << "decoding..." << std::endl;

    uint64_t num_decoded_ints = 0;
    uint64_t num_decoded_lists = 0;
    std::vector<double> timings;

    while (begin != end) {
        uint32_t n, universe;
        begin = header::read(begin, &n, &universe);

        auto start = clock_type::now();
        begin = Decoder::decode(begin, decoded.data(), universe, n);
        auto finish = clock_type::now();
        std::chrono::duration<double> elapsed = finish - start;
        timings.push_back(elapsed.count());
        num_decoded_ints += n;
        ++num_decoded_lists;

        // if (n > 4096) {
        //     auto start = clock_type::now();
        //     begin = Decoder::decode(begin, decoded.data(), universe, n);
        //     auto finish = clock_type::now();
        //     std::chrono::duration<double> elapsed = finish - start;
        //     timings.push_back(elapsed.count());
        //     num_decoded_ints += n;
        //     ++num_decoded_lists;
        // } else {
        //     begin = Decoder::decode(begin, decoded.data(), universe, n);
        // }
    }

    file.close();

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
}

template<typename Decoder, typename Dictionary>
void decode_dint(std::string const& type,
                 char const* encoded_data_filename,
                 char const* dictionary_filename)
{
    if (!dictionary_filename) {
        throw std::runtime_error("dictionary_filename must be specified");
    }

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

    Dictionary dict;
    size_t dictionaries_bytes = 0;
    typename Dictionary::builder builder;
    std::ifstream dictionary_file(dictionary_filename);
    dictionaries_bytes += builder.load(dictionary_file);
    builder.print_usage();
    builder.build(dict);
    dictionary_file.close();
    logger() << "Dictionary memory: " << double(dictionaries_bytes) / constants::MiB << " [MiB]" << std::endl;

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
        begin = Decoder::decode(dict, begin, decoded.data(), universe, n /*,stats*/);
        auto finish = clock_type::now();
        std::chrono::duration<double> elapsed = finish - start;
        timings.push_back(elapsed.count());
        num_decoded_ints += n;
        ++num_decoded_lists;
    }

    file.close();
    print_statistics(type, encoded_data_filename,
                     timings, num_decoded_ints, num_decoded_lists /*,stats*/);
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
    while (offset < num_bits)
    {
        uint64_t next_offset = bv.get_bits(offset, 64);
        offset += 64;
        universe = bv.get_bits(offset, 32);
        offset += 32;
        n = bv.get_bits(offset, 32);
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


    if (type == std::string("single_rect_dint")) {
        decode_dint<single_opt_dint, single_dictionary_rectangular_type>(
            type, encoded_data_filename, dictionary_filename
        );
    } else
    if (type == std::string("single_packed_dint")) {
        decode_dint<single_opt_dint, single_dictionary_packed_type>(
            type, encoded_data_filename, dictionary_filename
        );
    } else
    if (type == std::string("multi_packed_dint")) {
        decode_dint<multi_opt_dint, multi_dictionary_packed_type>(
            type, encoded_data_filename, dictionary_filename
        );
    } else
    if (type == std::string("pef")) {
        decode_pef(encoded_data_filename, freqs);
    }
    else {
        if (false) {
    #define LOOP_BODY(R, DATA, T)                                \
            } else if (type == BOOST_PP_STRINGIZE(T)) {          \
                decode<BOOST_PP_CAT(T, )>                        \
                    (type, encoded_data_filename);               \
                /**/

            BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, CODECS);
    #undef LOOP_BODY
        } else {
            logger() << "ERROR: unknown type '"
                     << type << "'" << std::endl;
        }
    }

    return 0;
}
