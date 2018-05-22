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
ds2i::dictionary::builder build_dict(std::ostream& os,block_stat_type& block_stats)
{
    ds2i::dictionary::builder dict_builder;
    dict_constructor_type::build(os,dict_builder,block_stats);
    dict_builder.prepare_for_encoding();
    return dict_builder;
}

struct block_enc_stats {
    std::vector<int64_t> dict_usage_lens;
    std::vector<int64_t> dict_usage;
    std::vector<int64_t> dict_entry_lens;
    std::vector<int64_t> dict_entry_freqs;
    std::vector<uint32_t> codes_per_block;
    std::vector<uint32_t> exceptions_per_block;
    std::vector<int32_t> block_size;
    std::vector<std::string> dict_entries;
    size_t num_blocks = 0;
    size_t total_codes = 0;
    size_t total_exception_codes_u16 = 0;
    size_t total_exception_codes_u32 = 0;
    size_t postings_encoded = 0;

    void init(ds2i::dictionary::builder& dict) {
        size_t capacity = dict.capacity();
        for(size_t i=0;i<capacity;i++) {
            dict_entry_freqs.push_back(dict.freq(i));
        }
        dict_usage.resize(capacity);
        dict_entry_lens.resize(capacity);
        dict_usage_lens.resize(256+1);
        dict_entry_lens[0] = 0;
        dict_entries.push_back("[exception code 16bit]");
        dict_entry_lens[1] = 0;
        dict_entries.push_back("[exception code 32bit]");
        dict_entry_lens[2] = 256;
        dict_entries.push_back("[1]*256");
        dict_entry_lens[3] = 128;
        dict_entries.push_back("[1]*128");
        dict_entry_lens[4] = 64;
        dict_entries.push_back("[1]*64");
        dict_entry_lens[5] = 32;
        dict_entries.push_back("[1]*32");
        for(size_t i=6;i<capacity;i++) {
            dict_entry_lens[i] = dict.size(i);
            dict_entries.push_back(dict.entry_string(i));
        }

    }

    template<class t_dict>
    void update(t_dict& dict,size_t input_size,std::vector<uint32_t>& codes,size_t num_codes,size_t num_exceptions)
    {
        num_blocks++;
        postings_encoded += input_size;
        total_codes += num_codes;
        block_size.push_back(input_size);
        codes_per_block.push_back(num_codes+num_exceptions);
        exceptions_per_block.push_back(num_exceptions);
        for(size_t i=0;i<num_codes;i++) {
            dict_usage[codes[i]]++;
            if(codes[i] < dict.special_cases()) {
                switch(codes[i]) {
                    case 0:
                        dict_usage_lens[0]++;
                        total_exception_codes_u16++;
                        break;
                    case 1:
                        dict_usage_lens[1]++;
                        total_exception_codes_u32++;
                        break;
                    case 2:
                        dict_usage_lens[256]++;
                        break;
                    case 3:
                        dict_usage_lens[128]++;
                        break;
                    case 4:
                        dict_usage_lens[64]++;
                        break;
                    case 5:
                        dict_usage_lens[32]++;
                        break;
                }
            } else {
                dict_usage_lens[dict.size(codes[i])]++;
            }

        }
    }

    void print(std::ostream &os,ds2i::dictionary::builder& dict) {
        os << "DICT LEN USAGE:\n";
        for(size_t i=0;i<dict_usage_lens.size();i++) {
            if(dict_usage_lens[i] != 0) {
                size_t encoded_nums = i;
                if(i == 0) encoded_nums = 1;
                os  << "\tlen = " << std::setw(3) << i
                    << "\tnum_codes = "  << std::setw(10) << dict_usage_lens[i]
                    << "\tnum_postings = " << std::setw(10) << dict_usage_lens[i] * encoded_nums
                    << "\tpercent of codes = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(dict_usage_lens[i]) / double(total_codes) * 100
                    << "\tpercent of postings = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(dict_usage_lens[i]*encoded_nums) / double(postings_encoded) * 100 << "\n";
            }
        }
        os << "CODE WORD USAGE:\n";
        size_t total_bits = 0;
        size_t exceptions_bits = 0;
        for(size_t i=0;i<dict_usage.size();i++) {
            int64_t encoded_nums = dict.size(i);
            int64_t shorts_used = 1;
            if(i == 0 || i == 1) {
                encoded_nums = 1;
                shorts_used = 2;
                if(i == 1) shorts_used = 3;
                exceptions_bits += 16*(shorts_used-1);
            }
            size_t cur_bits = dict_usage[i] * shorts_used * 16;
            total_bits += cur_bits;
            os  << "\tcode = " << std::setw(5) << i
                << "\tcode_u16 = " << std::setw(1) << shorts_used
                << "\tentry_len = " << std::setw(3) << dict.size(i)
                << "\tfreq = " << std::setw(12) << dict_usage[i]
                << "\tpredicted_freq = " << std::setw(12) << dict_entry_freqs[i]
                << "\tpercent of codes = " << std::setw(6) << std::fixed << std::setprecision(2) << double(dict_usage[i]) / double(total_codes) * 100
                << "\tpercent of postings = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(dict_usage[i]*encoded_nums) / double(postings_encoded) * 100
                << "\t" << dict.entry_string(i) << "\n";
        }

        std::sort(codes_per_block.begin(),codes_per_block.end());

        os << "CODES PER BLOCK:\n";
        size_t cur = codes_per_block[0];
        size_t freq = 1;
        for(size_t i=1;i<codes_per_block.size();i++) {
            if(codes_per_block[i] != cur) {
                os  << "\t\tnum_codes_in_block = " << std::setw(5) << cur
                    << "\tbpi = " << std::setw(6) << std::fixed <<  std::setprecision(3) << double(cur*16) / double(block_size[0])
                    << "\tfreq = " << std::setw(10) << freq
                    << "\tpercent of blocks = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(freq) / double(num_blocks) * 100 << "\n";
                freq = 1;
            } else {
                freq++;
            }
            cur = codes_per_block[i];
        }

        std::sort(exceptions_per_block.begin(),exceptions_per_block.end());

        os << "EXCEPTIONS PER BLOCK:\n";
        cur = exceptions_per_block[0];
        freq = 1;
        for(size_t i=1;i<exceptions_per_block.size();i++) {
            if(exceptions_per_block[i] != cur) {
                os  << "\t\texceptions_per_block = " << std::setw(5) << cur
                    << "\tfreq = " << std::setw(10) << freq
                    << "\tpercent of blocks = " << std::setw(6) << std::fixed <<  std::setprecision(2) << double(freq) / double(num_blocks) * 100 << "\n";
                freq = 1;
            } else {
                freq++;
            }
            cur = exceptions_per_block[i];
        }

        os << "POSTINGS = " << postings_encoded  << std::endl;
        os << "BPI = " <<  std::setprecision(5) << double(total_bits) / double(postings_encoded) << std::endl;
        os << "EXCEPTIONS_BPI = " << double(exceptions_bits) / double(postings_encoded) << std::endl;
        os << "TOTAL_NON_EXCEPTION_U16 = " << total_codes  << std::endl;
        os << "TOTAL_EXCEPTION_U16 = " << total_exception_codes_u16  << std::endl;
        os << "TOTAL_EXCEPTION_U32 = " << total_exception_codes_u32  << std::endl;
        os << "TOTAL_U16 = " << total_codes + total_exception_codes_u16 + (2 * total_exception_codes_u32)  << std::endl;
        os << "EXCEPTIONS BPI = " << double((total_exception_codes_u16 + (2 * total_exception_codes_u32)) * 16) / double(postings_encoded) << std::endl;
    }
};



struct encoding_stats {
    encoding_stats(ds2i::dictionary::builder& dict) {
        small_lists.init(dict);
        small_lists.init(dict);
        full_blocks.init(dict);
        nonfull_blocks.init(dict);
        all_blocks.init(dict);
    }

    block_enc_stats small_lists;
    block_enc_stats full_blocks;
    block_enc_stats nonfull_blocks;
    block_enc_stats all_blocks;

    void print(std::ostream& os,ds2i::dictionary::builder& dict) {
        all_blocks.print(os,dict);
    }
};


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
            auto bs = dint_block::encode(dict,buf.data(),n,output.data());
            stats.small_lists.update(dict,n,output,bs.written_codes,bs.written_exceptions);
            stats.all_blocks.update(dict,n,output,bs.written_codes,bs.written_exceptions);
        } else {
            size_t full_blocks = n / block_size;
            size_t left = n % block_size;
            for(size_t i=0;i<full_blocks;i++) {
                for(size_t j=0;j<block_size;j++) {
                    buf[j] = *itr - prev;
                    if(type == dict_type::docs) prev = *itr;
                    ++itr;
                }
                auto bs = dint_block::encode(dict,buf.data(),block_size,output.data());
                stats.full_blocks.update(dict,block_size,output,bs.written_codes,bs.written_exceptions);
                stats.all_blocks.update(dict,block_size,output,bs.written_codes,bs.written_exceptions);
            }


            if(left) {
                for(size_t j=0;j<left;j++) {
                    buf[j] = *itr - prev;
                    if(type == dict_type::docs) prev = *itr;
                    ++itr;
                }
                auto bs = dint_block::encode(dict,buf.data(),left,output.data());
                stats.nonfull_blocks.update(dict,left,output,bs.written_codes,bs.written_exceptions);
                stats.all_blocks.update(dict,left,output,bs.written_codes,bs.written_exceptions);
            }
        }
        progress += n+1;
    }
    return stats;
}


template<class dict_constructor_type,class block_stat_type,uint32_t encoding_block_size>
void eval_dict(std::string input_basename,std::string log_prefix)
{
    {
        std::ofstream log_file(log_prefix + "-docs-" + dict_constructor_type::type());

        auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::docs);
        auto dict = build_dict<block_stat_type,dict_constructor_type>(log_file,block_stats);

        dict.print_stats(log_file);

        auto enc_stats = encode_lists(dict,input_basename,dict_type::docs,encoding_block_size);

        enc_stats.print(log_file,dict);
    }
    {
        std::ofstream log_file(log_prefix + "-freqs-" + dict_constructor_type::type());

        auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::freqs);
        auto dict = build_dict<block_stat_type,dict_constructor_type>(log_file,block_stats);

        dict.print_stats(log_file);

        auto enc_stats = encode_lists(dict,input_basename,dict_type::freqs,encoding_block_size);

        enc_stats.print(log_file,dict);
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
    using block_stat_type = ds2i::block_stats_full_stride_geom<max_entry_width>;

    // PDF
    {
        using dict_type = ds2i::dint_dict_builder_PDF<block_stat_type,dict_entries, max_entry_width>;
        eval_dict<dict_type,block_stat_type,encoding_block_size>(input_basename,log_prefix);
    }
    // DSF
    {
        using dict_type = ds2i::dint_dict_builder_DSF<block_stat_type,dict_entries, max_entry_width>;
        eval_dict<dict_type,block_stat_type,encoding_block_size>(input_basename,log_prefix);
    }
    

    return 0;
}
