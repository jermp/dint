#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <numeric>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include "dictionary_types.hpp"
#include "block_statistics.hpp"
#include "block_codecs.hpp"
#include "logging.hpp"

using namespace ds2i;
using ds2i::logger;


template<class block_stat_type>
block_stat_type create_block_stats(std::string input_basename,dict_type type)
{
    std::string file_name = input_basename + ".docs";
    if(type == dict_type::freqs) file_name = input_basename + ".freqs";

    std::string block_stats_file = file_name + "." + block_stat_type::type();
    if( boost::filesystem::exists(block_stats_file) ) {
        return block_stat_type(block_stats_file);
    }

    binary_collection input(file_name.c_str());
    block_stat_type block_stats(input,type == dict_type::docs);
    block_stats.try_to_store(block_stats_file);
    return block_stats;
}

template<class block_stat_type,class dict_constructor_type>
ds2i::dictionary::builder build_dict(std::ostream& os,block_stat_type& block_stats)
{
    ds2i::dictionary::builder dict_builder;
    dict_constructor_type::build(os,dict_builder,block_stats);
    dict_builder.prepare_for_encoding();
    return dict_builder;
}

void encode_lists(ds2i::dictionary::builder& dict,std::string input_basename,dict_type type,size_t block_size)
{
    logger() << "encoding lists" << std::endl;
    //encoding_stats stats(dict);

    std::string file_name = input_basename + ".docs";
    if(type == dict_type::freqs) file_name = input_basename + ".freqs";

    binary_collection input(file_name.c_str());

    boost::progress_display progress(input.data_size());
    std::vector<uint32_t> buf(block_size);
    std::vector<uint32_t> output(block_size*3);
    size_t first = true;
    for (auto const& list: input) {
        if(type == dict_type::docs && first) {
            first = false;
            continue; // skip first doc list as it contains #docs in col
        }
        size_t n = list.size();
        auto itr = list.begin();
        uint32_t prev = 0;
        size_t full_blocks = n / block_size;
        size_t left = n % block_size;
        for(size_t i=0;i<full_blocks;i++) {
            for(size_t j=0;j<block_size;j++) {
                buf[j] = *itr - prev;
                if(type == dict_type::docs) prev = *itr;
                ++itr;
            }
            size_t written_codes = dint_block::encode(dict,buf.data(),block_size,output.data());
        }
        if(left) {
            for(size_t j=0;j<left;j++) {
                buf[j] = *itr - prev;
                if(type == dict_type::docs) prev = *itr;
                ++itr;
            }
            size_t written_codes = dint_block::encode(dict,buf.data(),left,output.data());
        }
        progress += n+1;
    }
}


template<class dict_constructor_type,class block_stat_type,uint32_t encoding_block_size>
void eval_dict(std::string input_basename,std::string log_prefix)
{
    {
        auto log = ds2i::start_logging(log_prefix + "-docs-" + dict_constructor_type::type());
        auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::docs);
        auto dict = build_dict<block_stat_type,dict_constructor_type>(block_stats);
        encode_lists(dict,input_basename,dict_type::docs,encoding_block_size);
        stop_logging(log);
    }
    {
        auto log = ds2i::start_logging(log_prefix + "-freqs-" + dict_constructor_type::type());
        auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::freqs);
        auto dict = build_dict<block_stat_type,dict_constructor_type>(block_stats);
        encode_lists(dict,input_basename,dict_type::freqs,encoding_block_size);
        stop_logging(log);
    }
}

int main(int argc, const char** argv) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << ":\n"
                  << "\t<collection basename> <log_prefix>"
                  << std::endl;
        return 1;
    }

    const char* input_basename = argv[1];
    std::string log_prefix = argv[2];

    // define dict type
    const uint32_t encoding_block_size = 256;
    const uint32_t max_entry_width = 16;
    const uint32_t dict_entries = 65536;
    
    // DSF
    {
        using block_stat_type = ds2i::block_statistics<max_entry_width,ds2i::stats_geometric>;
        using dict_type = ds2i::dint_dict_type_DSF<block_stat_type,dict_entries, max_entry_width>;
        eval_dict<dict_type,block_stat_type,encoding_block_size>(input_basename,log_prefix);
    }
    // // // PDF
    // {
    //     using block_stat_type = ds2i::block_stats_full_stride_geom<max_entry_width>;
    //     using dict_type = ds2i::dint_dict_builder_PDF<block_stat_type,dict_entries, max_entry_width>;
    //     eval_dict<dict_type,block_stat_type,encoding_block_size>(input_basename,log_prefix);
    // }

    // // SDF
    // {
    //     using block_stat_type = ds2i::block_stats_full_stride_linear<max_entry_width>;
    //     using dict_type = ds2i::dint_dict_builder_SDF<block_stat_type,dict_entries, max_entry_width>;
    //     eval_dict<dict_type,block_stat_type,encoding_block_size>(input_basename,log_prefix);
    // }

    return 0;
}
