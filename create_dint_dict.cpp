#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <numeric>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include "dictionary_builders.hpp"
#include "block_stats.hpp"
#include "block_codecs.hpp"

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
ds2i::dictionary::builder build_dict(block_stat_type& block_stats)
{
    ds2i::dictionary::builder dict_builder;
    dict_constructor_type::build(dict_builder,block_stats);
    dict_builder.prepare_for_encoding();
    return dict_builder;
}

struct block_enc_stats {
    std::vector<uint64_t> dict_usage_lens;
    std::vector<uint64_t> dict_usage;
    std::vector<uint64_t> dict_entry_lens;
    std::vector<uint32_t> codes_per_block;
    std::vector<uint32_t> block_size;
    std::vector<std::string> dict_entries;
    size_t num_blocks = 0;
    size_t total_codes = 0;
    size_t postings_encoded = 0;

    void init(ds2i::dictionary::builder& dict) {
        size_t capacity = dict.capacity();
        dict_usage.resize(capacity);
        dict_entry_lens.resize(capacity);
        dict_usage_lens.resize(256+1);
        dict_entry_lens[0] = 0;
        dict_entries.push_back("[exception code]");
        dict_entry_lens[1] = 256;
        dict_entries.push_back("[1]*256");
        dict_entry_lens[2] = 128;
        dict_entries.push_back("[1]*128");
        dict_entry_lens[3] = 64;
        dict_entries.push_back("[1]*64");
        dict_entry_lens[4] = 32;
        dict_entries.push_back("[1]*32");
        for(size_t i=5;i<capacity;i++) {
            dict_entry_lens[i] = dict.size(i);
            dict_entries.push_back(dict.entry_string(i));
        }

    }

    template<class t_dict>
    void update(t_dict& dict,size_t input_size,std::vector<uint32_t>& codes,size_t num_codes)
    {
        num_blocks++;
        postings_encoded += input_size;
        total_codes += num_codes;
        block_size.push_back(input_size);
        codes_per_block.push_back(num_codes);
        for(size_t i=0;i<num_codes;i++) {
            dict_usage[codes[i]]++;
            if(codes[i] < dict.special_cases()) {
                switch(codes[i]) {
                    case 0:
                        dict_usage_lens[0]++;
                        break;
                    case 1:
                        dict_usage_lens[256]++;
                        break;
                    case 2:
                        dict_usage_lens[128]++;
                        break;
                    case 3:
                        dict_usage_lens[64]++;
                        break;
                    case 4:
                        dict_usage_lens[32]++;
                        break;
                    case 5:
                        dict_usage_lens[16]++;
                        break;
                }
            } else {
                dict_usage_lens[dict.size(codes[i])]++;
            }

        }
    }
};

std::ostream &operator<<(std::ostream &os, block_enc_stats const &stats) {
    os << "DICT LEN USAGE:\n";
    for(size_t i=0;i<stats.dict_usage_lens.size();i++) {
        if(stats.dict_usage_lens[i] != 0) {
            size_t encoded_nums = i;
            if(i == 0) encoded_nums = 1;
            os  << "\tlen = " << std::setw(3) << i
                << "\tnum_codes = "  << std::setw(10) << stats.dict_usage_lens[i]
                << "\tnum_postings = " << std::setw(10) << stats.dict_usage_lens[i] * encoded_nums
                << "\tpercent of codes = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(stats.dict_usage_lens[i]) / double(stats.total_codes) * 100
                << "\tpercent of postings = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(stats.dict_usage_lens[i]*encoded_nums) / double(stats.postings_encoded) * 100 << "\n";
        }
    }
    os << "CODE WORD USAGE:\n";
    for(size_t i=0;i<stats.dict_usage.size();i++) {
        size_t encoded_nums = stats.dict_entry_lens[i];
        if(i == 0) encoded_nums = 1;
        os  << "\tcode = " << std::setw(5) << i
            << "\tcode_len = " << std::setw(3) << stats.dict_entry_lens[i]
            << "\tfreq = " << std::setw(10) << stats.dict_usage[i]
            << "\tpercent of codes = " << std::setw(6) << std::fixed << std::setprecision(2) << double(stats.dict_usage[i]) / double(stats.total_codes) * 100
            << "\tpercent of postings = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(stats.dict_usage[i]*encoded_nums) / double(stats.postings_encoded) * 100
            << "\t" << stats.dict_entries[i] << "\n";
    }

    auto codes_per_block = stats.codes_per_block;
    std::sort(codes_per_block.begin(),codes_per_block.end());

    os << "CODES PER BLOCK:\n";
    size_t cur = codes_per_block[0];
    size_t freq = 1;
    for(size_t i=1;i<codes_per_block.size();i++) {
        if(codes_per_block[i] != cur) {
            os  << "\t\tnum_codes_in_block = " << std::setw(10) << cur
                << "\tbpi = " << std::setw(6) << std::fixed <<  std::setprecision(3) << double(cur*16) / double(stats.block_size[0])
                << "\tfreq = " << std::setw(10) << freq
                << "\tpercent of blocks = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(freq) / double(stats.num_blocks) * 100 << "\n";
            freq = 1;
        } else {
            freq++;
        }
        cur = codes_per_block[i];
    }


    return os;
}

struct encoding_stats {
    encoding_stats(ds2i::dictionary::builder& dict) {
        small_lists.init(dict);
        full_blocks.init(dict);
        nonfull_blocks.init(dict);
    }

    block_enc_stats small_lists;
    block_enc_stats full_blocks;
    block_enc_stats nonfull_blocks;
};

std::ostream &operator<<(std::ostream &os, encoding_stats const &stats) {
    return os << "small_lists: \n" << stats.small_lists << "\n"
              << "full_blocks: \n" << stats.full_blocks << "\n"
              << "nonfull_blocks: \n" << stats.nonfull_blocks << "\n";
}



encoding_stats encode_lists(ds2i::dictionary::builder& dict,std::string input_basename,dict_type type,size_t block_size)
{
    logger() << "encoding lists" << std::endl;
    encoding_stats stats(dict);

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
        if(n < block_size) { // small list
            for(size_t j=0;j<n;j++) {
                    buf[j] = *itr - prev;
                    if(type == dict_type::docs) prev = *itr;
                    ++itr;
            }
            auto written_codes = dint_block::encode(dict,buf.data(),n,output.data());
            stats.small_lists.update(dict,n,output,written_codes);
        } else {
            size_t full_blocks = n / block_size;
            size_t left = n % block_size;
            for(size_t i=0;i<full_blocks;i++) {
                for(size_t j=0;j<block_size;j++) {
                    buf[j] = *itr - prev;
                    if(type == dict_type::docs) prev = *itr;
                    ++itr;
                }
                auto written_codes = dint_block::encode(dict,buf.data(),block_size,output.data());
                stats.full_blocks.update(dict,block_size,output,written_codes);
            }


            if(left) {
                for(size_t j=0;j<left;j++) {
                    buf[j] = *itr - prev;
                    if(type == dict_type::docs) prev = *itr;
                    ++itr;
                }
                auto written_codes = dint_block::encode(dict,buf.data(),block_size,output.data());
                stats.nonfull_blocks.update(dict,block_size,output,written_codes);
            }
        }
        progress += n+1;
    }
    return stats;
}


int main(int argc, const char** argv) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << ":\n"
                  << "\t<collection basename>"
                  << std::endl;
        return 1;
    }

    const char* input_basename = argv[1];

    // define dict type
    const uint32_t encoding_block_size = 256;
    const uint32_t max_entry_width = 16;
    const uint32_t dict_entries = 65536;
    using block_stat_type = ds2i::block_stats_full_stride_geom<max_entry_width>;
    using dict_constructor_type = ds2i::dint_dict_builder_smc<block_stat_type,dict_entries, max_entry_width>;

    // first docs
    {
        auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::docs);

        auto dict = build_dict<block_stat_type,dict_constructor_type>(block_stats);

        auto enc_stats = encode_lists(dict,input_basename,dict_type::docs,encoding_block_size);

        std::cout << enc_stats << std::endl;
    }
    // second freqs
    {
        auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::freqs);

        auto dict = build_dict<block_stat_type,dict_constructor_type>(block_stats);

        auto enc_stats = encode_lists(dict,input_basename,dict_type::freqs,encoding_block_size);

        std::cout << enc_stats << std::endl;
    }


    return 0;
}
