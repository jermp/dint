#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <numeric>

#include <succinct/mapper.hpp>

#include "configuration.hpp"
#include "index_types.hpp"
#include "util.hpp"
#include "verify_collection.hpp"
#include "index_build_utils.hpp"

using namespace ds2i;

template<typename Collection>
void dump_index_specific_stats(Collection const&, std::string const&)
{}

void dump_index_specific_stats(uniform_index const& coll,
                               std::string const& type)
{
    stats_line()
        ("type", type)
        ("log_partition_size", int(coll.params().log_partition_size))
        ;
}


void dump_index_specific_stats(opt_index const& coll,
                               std::string const& type)
{
    auto const& conf = configuration::get();

    double long_postings = 0;
    double docs_partitions = 0;
    double freqs_partitions = 0;

    for (size_t s = 0; s < coll.size(); ++s) {
        auto const& list = coll[s];
        if (list.size() > MIN_SIZE) {
            long_postings += list.size();
            docs_partitions += list.docs_enum().num_partitions();
            freqs_partitions += list.freqs_enum().base().num_partitions();
        }
    }

    stats_line()
        ("type", type)
        ("eps1", conf.eps1)
        ("eps2", conf.eps2)
        ("fix_cost", conf.fix_cost)
        ("docs_avg_part", long_postings / docs_partitions)
        ("freqs_avg_part", long_postings / freqs_partitions)
        ;
}

template<typename InputCollection, typename Builder>
void build_model(InputCollection const&, Builder&, const char*)
{}

template<typename InputCollection>
void build_model(InputCollection const& input,
                 block_dint_index::builder& builder, const char* output_filename)
{
    builder.build_model(input, std::string(output_filename));
}

template <typename InputCollection, typename CollectionType>
void create_collection(InputCollection const& input,
                       global_parameters const& params,
                       const char* output_filename, bool check,
                       std::string const& seq_type)
{
    DS2I_LOG << "Processing " << input.num_docs() << " documents";
    double tick = get_time_usecs();
    double user_tick = get_user_time_usecs();

    typename CollectionType::builder builder(input.num_docs(), params);
    build_model(input, builder, output_filename);

    progress_logger plog("Encoded");
    DS2I_LOG << "Encoding...";
    for (auto const& plist: input) {
        // if (plist.docs.size() > MIN_SIZE) {
            uint64_t freqs_sum = std::accumulate(plist.freqs.begin(),
                                                 plist.freqs.end(), uint64_t(0));
            builder.add_posting_list(plist.docs.size(), plist.docs.begin(),
                                     plist.freqs.begin(), freqs_sum);
            plog.done_sequence(plist.docs.size());
        // }
    }

    plog.log();
    CollectionType coll;
    builder.build(coll);
    double elapsed_secs = (get_time_usecs() - tick) / 1000000;
    double user_elapsed_secs = (get_user_time_usecs() - user_tick) / 1000000;
    DS2I_LOG << seq_type << " collection built in "
             << elapsed_secs << " seconds";

    stats_line()
        ("type", seq_type)
        ("worker_threads", configuration::get().worker_threads)
        ("construction_time", elapsed_secs)
        ("construction_user_time", user_elapsed_secs)
        ;

    dump_stats(coll, seq_type, plog.postings);
    dump_index_specific_stats(coll, seq_type);

    if (output_filename) {
        succinct::mapper::freeze(coll, output_filename);
        if (check) {
            verify_collection<InputCollection, CollectionType>(input, output_filename);
        }
    }
}

int main(int argc, const char** argv) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << ":\n"
                  << "\t<index type> <collection basename> [<output filename>] [--check]"
                  << std::endl;
        return 1;
    }

    std::string type = argv[1];
    const char* input_basename = argv[2];
    const char* output_filename = nullptr;
    if (argc > 3) {
        output_filename = argv[3];
    }

    bool check = false;
    if (argc > 4 && std::string(argv[4]) == "--check") {
        check = true;
    }

    binary_freq_collection input(input_basename);
    ds2i::global_parameters params;
    params.log_partition_size = configuration::get().log_partition_size;

    if (false) {
#define LOOP_BODY(R, DATA, T)                                   \
        } else if (type == BOOST_PP_STRINGIZE(T)) {             \
            create_collection<binary_freq_collection,           \
                              BOOST_PP_CAT(T, _index)>          \
                (input, params, output_filename, check, type);  \
            /**/

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, DS2I_INDEX_TYPES);
#undef LOOP_BODY
    } else {
        DS2I_LOG << "ERROR: Unknown type " << type;
    }

    return 0;
}
