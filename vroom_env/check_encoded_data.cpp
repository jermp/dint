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
#include "dint_codecs.hpp"
#include "util.hpp"

using namespace ds2i;

template <typename Decoder, typename Dictionary>
void check_dint(char const* collection_filename,
                char const* encoded_data_filename,
                char const* dictionary_filename) {
    if (!dictionary_filename) {
        throw std::runtime_error("dictionary_filename must be specified");
    }

    boost::iostreams::mapped_file_source file;
    file.open(encoded_data_filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening index file");
    }

    uint8_t const* begin = (uint8_t const*)file.data();
    uint64_t size = file.size() / sizeof(uint8_t);
    auto ret = posix_madvise((void*)begin, size, POSIX_MADV_SEQUENTIAL);
    if (ret) {
        logger() << "Error calling madvice: " << errno << std::endl;
    }

    binary_collection input(collection_filename);
    auto it = input.begin();

    Dictionary dict;
    size_t dictionaries_bytes = 0;
    typename Dictionary::builder builder;
    std::ifstream dictionary_file(dictionary_filename);
    dictionaries_bytes += builder.load(dictionary_file);
    builder.print_usage();
    builder.build(dict);
    dictionary_file.close();
    logger() << "Dictionary memory: "
             << double(dictionaries_bytes) / constants::MiB << " [MiB]"
             << std::endl;

    std::vector<uint32_t> decoded;
    decoded.resize(constants::max_size, 0);

    bool docs = true;
    boost::filesystem::path collection_path(collection_filename);
    if (collection_path.extension() == ".freqs") {
        docs = false;
        logger() << "checking freqs..." << std::endl;
    } else if (collection_path.extension() == ".docs") {
        ++it;  // skip first singleton sequence, containing num. of docs
        logger() << "checking docs..." << std::endl;
    } else {
        throw std::runtime_error("unsupported file format");
    }

    uint64_t total_decoded_ints = 0;
    uint64_t sequence = 0;

    dint_statistics stats;

    for (; it != input.end(); ++it) {
        auto const& list = *it;
        uint32_t size = list.size();
        if (size > constants::min_size) {
            uint32_t n, universe;
            begin = header::read(begin, &n, &universe);

            if (n != size) {
                std::cerr << "sequence has wrong length: got " << n
                          << " but expected " << sequence << std::endl;
            }

            begin = Decoder::decode(dict, begin, decoded.data(), universe, n);
            total_decoded_ints += n;

            uint32_t prev = docs ? -1 : 0;
            uint64_t j = 0;
            for (auto b = list.begin(); b != list.end(); ++b, ++j) {
                uint32_t expected = *b - prev - 1;
                if (docs) {
                    prev = *b;
                }
                if (decoded[j] != expected) {
                    std::cerr << "Sequence " << sequence
                              << ": error at position " << j << "/" << n
                              << " (got " << decoded[j] << " but expected "
                              << expected << ")" << std::endl;
                }
                decoded[j] = 0;
            }

            for (; j != n + constants::max_entry_size; ++j) {
                decoded[j] = 0;
            }
        }

        ++sequence;
    }

    logger() << "checked " << total_decoded_ints << " integers: OK!"
             << std::endl;

    file.close();
}

int main(int argc, char** argv) {
    int mandatory = 4;
    if (argc < mandatory) {
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t<type> <collection_filename> <encoded_data_filename> "
                     "[--dict <dictionary_filename>]"
                  << std::endl;
        return 1;
    }

    using namespace ds2i;
    std::string type = argv[1];
    char const* collection_filename = argv[2];
    char const* encoded_data_filename = argv[3];
    char const* dictionary_filename = nullptr;

    for (int i = mandatory; i < argc; ++i) {
        if (argv[i] == std::string("--dict")) {
            ++i;
            dictionary_filename = argv[i];
        } else {
            throw std::runtime_error("unknown parameter");
        }
    }

    if (type == std::string("single_rect_dint")) {
        check_dint<single_opt_dint, single_dictionary_rectangular_type>(
            collection_filename, encoded_data_filename, dictionary_filename);
    } else if (type == std::string("single_packed_dint")) {
        check_dint<single_opt_dint, single_dictionary_packed_type>(
            collection_filename, encoded_data_filename, dictionary_filename);
    } else if (type == std::string("multi_packed_dint")) {
        check_dint<multi_opt_dint, multi_dictionary_packed_type>(
            collection_filename, encoded_data_filename, dictionary_filename);
    }

    return 0;
}
