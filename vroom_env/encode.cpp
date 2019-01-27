#include <iostream>
#include <algorithm>
#include <fstream>

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <succinct/mapper.hpp>

#include "dint_configuration.hpp"
#include "codecs.hpp"
#include "dint_codecs.hpp"
#include "util.hpp"
#include "hash_utils.hpp"
#include "binary_collection.hpp"
#include "semiasync_queue.hpp"
#include "jobs.hpp"

using namespace ds2i;

typedef binary_collection::posting_type const* iterator_type;
const uint32_t num_jobs = 1 << 24;

void save_if(char const* output_filename,
             std::vector<uint8_t> const& output)
{
    if (output_filename) {
        logger() << "writing encoded data..." << std::endl;
        std::ofstream output_file(output_filename);
        output_file.write(reinterpret_cast<char const*>(output.data()),
                          output.size() * sizeof(output[0]));
        output_file.close();
        logger() << "DONE" << std::endl;
    }
}

void print_statistics(std::string type, char const* collection_name,
                      std::vector<uint8_t> const& output,
                      uint64_t num_total_ints,
                      uint64_t num_processed_lists)
{
    double GiB_space = output.size() * 1.0 / constants::GiB;
    double bpi_space = output.size() * sizeof(output[0]) * 8.0 / num_total_ints;

    logger() << "encoded " << num_processed_lists << " lists" << std::endl;
    logger() << "encoded " << num_total_ints << " integers" << std::endl;
    logger() << GiB_space << " [GiB]" << std::endl;
    logger() << "bits x integer: " << bpi_space << std::endl;

    // stats to std output
    std::cout << "{";
    std::cout << "\"filename\": \"" << collection_name << "\", ";
    std::cout << "\"num_sequences\": \"" << num_processed_lists << "\", ";
    std::cout << "\"num_integers\": \"" << num_total_ints << "\", ";
    std::cout << "\"type\": \"" << type << "\", ";
    std::cout << "\"GiB\": \"" << GiB_space << "\", ";
    std::cout << "\"bpi\": \"" << bpi_space << "\"";
    std::cout << "}" << std::endl;
}

template<typename Encoder>
void encode(std::string const& type,
            char const* collection_name,
            char const* output_filename)
{
    binary_collection input(collection_name);

    auto it = input.begin();
    uint64_t num_processed_lists = 0;
    uint64_t num_total_ints = 0;

    uint64_t total_progress = input.num_postings();
    bool docs = true;
    boost::filesystem::path collection_path(collection_name);
    if (collection_path.extension() == ".freqs") {
        docs = false;
        logger() << "encoding freqs..." << std::endl;
    } else if (collection_path.extension() == ".docs") {
        // skip first singleton sequence, containing num. of docs
        ++it;
        total_progress -= 2;
        logger() << "encoding docs..." << std::endl;
    } else {
        throw std::runtime_error("unsupported file format");
    }

    std::vector<uint8_t> output;
    uint64_t bytes = 5 * constants::GiB;
    output.reserve(bytes);

    std::vector<uint32_t> buf;
    boost::progress_display progress(total_progress);
    semiasync_queue jobs_queue(num_jobs);

    for (; it != input.end(); ++it) {
        auto const& list = *it;
        uint32_t n = list.size();

        // 1. sequential construction
        std::vector<uint32_t> buf;
        buf.reserve(n);
        uint64_t universe = 0;
        uint32_t prev = docs ? -1 : 0;
        auto begin = list.begin();
        for (uint64_t i = 0; i != n; ++i, ++begin) {
            buf.push_back(*begin - prev - 1);
            if (docs) prev = *begin;
            universe += buf.back();
        }
        assert(buf.size() == n);

        header::write(n, universe, output);
        Encoder::encode(
            buf.data(), universe, buf.size(), output
        );
        progress += n + 1;
        ++num_processed_lists;
        num_total_ints += n;



        // 2. parallel construction
        // std::shared_ptr<sequence_adder<iterator_type, Encoder>>
        //     ptr(new sequence_adder<iterator_type, Encoder>(
        //         list.begin(), n,
        //         progress, output, docs,
        //         num_processed_lists, num_total_ints
        //     )
        // );
        // jobs_queue.add_job(ptr, n);
    }

    jobs_queue.complete();
    print_statistics(type, collection_name, output,
                     num_total_ints, num_processed_lists);
    save_if(output_filename, output);
}

template<typename Encoder, typename Dictionary>
void encode_dint(std::string const& type,
                 char const* collection_name,
                 char const* output_filename,
                 char const* dictionary_filename)
{
    binary_collection input(collection_name);

    auto it = input.begin();
    uint64_t num_processed_lists = 0;
    uint64_t num_total_ints = 0;

    typedef typename Dictionary::builder Builder;
    Builder builder;

    if (dictionary_filename) {
        std::ifstream dictionary_file(dictionary_filename);
        builder.load(dictionary_file);
        logger() << "preparing for encoding..." << std::endl;
        builder.prepare_for_encoding();
        builder.print_usage();
    }

    uint64_t total_progress = input.num_postings();
    bool docs = true;
    boost::filesystem::path collection_path(collection_name);
    if (collection_path.extension() == ".freqs") {
        docs = false;
        logger() << "encoding freqs..." << std::endl;
    } else if (collection_path.extension() == ".docs") {
        // skip first singleton sequence, containing num. of docs
        ++it;
        total_progress -= 2;
        logger() << "encoding docs..." << std::endl;
    } else {
        throw std::runtime_error("unsupported file format");
    }

    std::vector<uint8_t> output;
    uint64_t bytes = 5 * constants::GiB;
    output.reserve(bytes);

    std::vector<uint32_t> buf;
    boost::progress_display progress(total_progress);
    semiasync_queue jobs_queue(num_jobs);

    for (; it != input.end(); ++it) {
        auto const& list = *it;
        uint32_t n = list.size();
        std::shared_ptr<single_dict_sequence_adder<iterator_type, Encoder, Builder>>
            ptr(new single_dict_sequence_adder<iterator_type, Encoder, Builder>(
                list.begin(), n,
                builder,
                progress, output, docs,
                num_processed_lists, num_total_ints
            )
        );
        jobs_queue.add_job(ptr, n);
    }

    jobs_queue.complete();
    print_statistics(type, collection_name, output,
                     num_total_ints, num_processed_lists);
    save_if(output_filename, output);
}

void encode_pef(char const* collection_name,
                char const* output_filename)
{
    binary_collection input(collection_name);

    auto it = input.begin();
    uint64_t num_processed_lists = 0;
    uint64_t num_total_ints = 0;

    uint64_t total_progress = input.num_postings();
    bool docs = true;
    boost::filesystem::path collection_path(collection_name);
    if (collection_path.extension() == ".freqs") {
        docs = false;
        logger() << "encoding freqs..." << std::endl;
    } else if (collection_path.extension() == ".docs") {
        // skip first singleton sequence, containing num. of docs
        ++it;
        total_progress -= 2;
        logger() << "encoding docs..." << std::endl;
    } else {
        throw std::runtime_error("unsupported file format");
    }

    succinct::bit_vector_builder bvb;
    boost::progress_display progress(total_progress);
    semiasync_queue jobs_queue(num_jobs);

    for (; it != input.end(); ++it)
    {
        auto const& list = *it;
        uint32_t n = list.size();

        // 1. sequential version
        // if (n > constants::min_size) {
        //     pef::encode(list.begin(), list.back(), n, bvb, not docs);
        //     ++num_processed_lists;
        //     num_total_ints += n;
        //     progress += n + 1;
        // }

        // 2. parallel version
        if (n > constants::min_size) {
            std::shared_ptr<pef_sequence_adder<iterator_type>>
                ptr(new pef_sequence_adder<iterator_type>(
                    list.begin(),
                    n, list.back() + 1, bvb,
                    progress, docs,
                    num_processed_lists, num_total_ints
                )
            );
            jobs_queue.add_job(ptr, n);
        }
    }

    jobs_queue.complete();

    double GiB_space = (bvb.size() + 7.0) / 8.0 / constants::GiB;
    double bpi_space = double(bvb.size()) / num_total_ints;

    logger() << "encoded " << num_processed_lists << " lists" << std::endl;
    logger() << "encoded " << num_total_ints << " integers" << std::endl;
    logger() << GiB_space << " [GiB]" << std::endl;
    logger() << "bits x integer: " << bpi_space << std::endl;

    // stats to std output
    std::cout << "{";
    std::cout << "\"filename\": \"" << collection_name << "\", ";
    std::cout << "\"num_sequences\": \"" << num_processed_lists << "\", ";
    std::cout << "\"num_integers\": \"" << num_total_ints << "\", ";
    std::cout << "\"type\": \"pef\", ";
    std::cout << "\"GiB\": \"" << GiB_space << "\", ";
    std::cout << "\"bpi\": \"" << bpi_space << "\"";
    std::cout << "}" << std::endl;

    if (output_filename) {
        logger() << "writing encoded data..." << std::endl;
        succinct::bit_vector bv(&bvb);
        succinct::mapper::freeze(bv, output_filename);
        logger() << "DONE" << std::endl;
    }
}

int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t<type> <collection_name> [--dict <dictionary_filename>] [--out <output_filename>]"
                  << std::endl;
        return 1;
    }

    std::string type = argv[1];
    char const* collection_name = argv[2];
    char const* dictionary_filename = nullptr;
    char const* output_filename = nullptr;

    std::string cmd(std::string(argv[0]) + " " + type + " " + std::string(collection_name));

    for (int i = 3; i < argc; ++i) {
        if (argv[i] == std::string("--dict")) {
            ++i;
            dictionary_filename = argv[i];
            cmd += " --dict " + std::string(dictionary_filename);
        } else if (argv[i] == std::string("--out")) {
            ++i;
            output_filename = argv[i];
            cmd += " --out " + std::string(output_filename);
        } else {
            throw std::runtime_error("unknown parameter");
        }
    }

    logger() << cmd << std::endl;

    if (type == std::string("single_rect_dint")) {
        encode_dint<single_opt_dint, single_dictionary_rectangular_type>(
            type, collection_name, output_filename, dictionary_filename
        );
    } else
    if (type == std::string("single_packed_dint")) {
        encode_dint<single_opt_dint, single_dictionary_packed_type>(
            type, collection_name, output_filename, dictionary_filename
        );
    } else
    if (type == std::string("multi_packed_dint")) {
        encode_dint<multi_opt_dint, multi_dictionary_packed_type>(
            type, collection_name, output_filename, dictionary_filename
        );
    } else
    if (type == std::string("pef")) {
        encode_pef(collection_name, output_filename);
    }
    else {
        if (false) {
    #define LOOP_BODY(R, DATA, T)                                \
            } else if (type == BOOST_PP_STRINGIZE(T)) {          \
                encode<BOOST_PP_CAT(T, )>                        \
                    (type, collection_name, output_filename);    \
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
