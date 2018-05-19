#include <iostream>
#include <fstream>
#include <algorithm>

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/filesystem.hpp>

#include <boost/iostreams/device/mapped_file.hpp>
#include <sys/mman.h>

#include "binary_collection.hpp"
#include "codecs.hpp"
#include "util.hpp"
#include "dictionary.hpp"

using namespace ds2i;

template<typename Decoder>
void check(char const* collection_filename,
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
    auto ret = posix_madvise((void*) begin, size, POSIX_MADV_SEQUENTIAL);
    if (ret) {
        logger() << "Error calling madvice: " << errno << std::endl;
    }

    binary_collection input(collection_filename);

    dictionary dict;
    if (dictionary_filename) {
        dictionary::builder builder;
        std::ifstream dictionary_file(dictionary_filename);
        builder.load(dictionary_file);
        builder.build(dict);
    }

    const static uint64_t MAX_SIZE = 30000000;
    std::vector<uint32_t> decoded;
    decoded.resize(MAX_SIZE, 1);

    bool take_gaps = true;
    boost::filesystem::path collection_path(collection_filename);
    if (collection_path.extension() == ".freqs") {
        take_gaps = false;
        logger() << "not taking d-gaps" << std::endl;
    } else if (collection_path.extension() == ".docs") {
        logger() << "taking d-gaps" << std::endl;
    } else {
        throw std::runtime_error("unsupported file format");
    }

    logger() << "checking encoded data..." << std::endl;

    uint64_t total_decoded_ints = 0;
    uint64_t sequence = 0;

    for (auto const& list: input)
    {
        uint32_t size = list.size();
        if (size > MIN_SIZE)
        {
            uint32_t n, universe;
            begin = header::read(begin, &n, &universe);
            if (n != size) {
                std::cerr << "sequence has wrong length: got "
                          << n << " but expected " << sequence << std::endl;
            }

            begin = Decoder::decode(begin,
                                    decoded.data(),
                                    universe, n, &dict);
            total_decoded_ints += n;

            uint32_t prev = 0;
            uint64_t j = 0;
            for (auto b = list.begin(); b != list.end(); ++b, ++j) {
                uint32_t expected = *b - prev;
                if (take_gaps) {
                    prev = *b;
                }
                if (decoded[j] != expected) {
                    std::cerr << "Error at position " << j << "/" << n << ": got " << decoded[j] << " but expected " << expected << std::endl;
                }
                // std::cerr << expected << " ";
            }
            std::cerr << std::endl;
            std::exit(1);
        }

        ++sequence;
    }

    logger() << "checked " << total_decoded_ints << " integers: OK!" << std::endl;

    file.close();
}

int main(int argc, char** argv) {

    if (argc < 4) {
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t<type> <collection_filename> <encoded_data_filename> [--dict <dictionary_filename>]"
                  << std::endl;
        return 1;
    }

    using namespace ds2i;
    std::string type = argv[1];
    char const* collection_filename = argv[2];
    char const* encoded_data_filename = argv[3];
    char const* dictionary_filename = nullptr;

    for (int i = 4; i < argc; ++i) {
        if (argv[i] == std::string("--dict")) {
            ++i;
            dictionary_filename = argv[i];
        } else {
            throw std::runtime_error("unknown parameter");
        }
    }

    if (false) {
#define LOOP_BODY(R, DATA, T)                                   \
        } else if (type == BOOST_PP_STRINGIZE(T)) {             \
            check<BOOST_PP_CAT(T, )>                            \
                (collection_filename, encoded_data_filename,    \
                 dictionary_filename);                          \
            /**/

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, CODECS);
#undef LOOP_BODY
    } else {
        logger() << "ERROR: unknown type '"
                 << type << "'" << std::endl;
    }

    return 0;
}
