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
void decode(char const* encoded_data_filename,
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

    const static uint64_t MAX_SIZE = 30000000;
    std::vector<uint32_t> decoded;
    decoded.resize(MAX_SIZE);

    logger() << "decoding..." << std::endl;

    uint64_t total_decoded_ints = 0;
    std::vector<double> timings;

    while (begin != end) {
        uint32_t n, universe;
        begin = header::read(begin, &n, &universe);
        auto start = clock_type::now();
        begin = Decoder::decode(begin, decoded.data(), universe, n, &dict);
        auto finish = clock_type::now();
        std::chrono::duration<double> elapsed = finish - start;
        timings.push_back(elapsed.count());
        total_decoded_ints += n;
    }

    double tot_elapsed = std::accumulate(timings.begin(), timings.end(), double(0.0));
    logger() << "elapsed time " << tot_elapsed << " [sec]" << std::endl;
    double ns_x_int = tot_elapsed * 1000000000 / total_decoded_ints;
    logger() << ns_x_int << " [ns] x int" << std::endl;
    logger() << 1 / ns_x_int * 1000000000 << " ints x [sec]" << std::endl;

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

    if (false) {
#define LOOP_BODY(R, DATA, T)                                    \
        } else if (type == BOOST_PP_STRINGIZE(T)) {              \
            decode<BOOST_PP_CAT(T, )>                            \
                (encoded_data_filename, dictionary_filename);    \
            /**/

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, CODECS);
#undef LOOP_BODY
    } else {
        logger() << "ERROR: unknown type '"
                 << type << "'" << std::endl;
    }

    return 0;
}
