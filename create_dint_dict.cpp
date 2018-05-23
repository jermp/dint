#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <numeric>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/format.hpp>

#include "dictionary_types.hpp"
#include "block_statistics.hpp"
#include "block_codecs.hpp"
#include "logging.hpp"

using namespace ds2i;


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

struct encoding_stats {
    encoding_stats(ds2i::dictionary::builder& dict_,size_t block_size_)
        : dict(dict_) , block_size(block_size_)
    {
        code_usage.resize(dict.capacity()+1);
        codes_per_block.resize(block_size*4);
        exceptions_per_block.resize(block_size*4);
    }

    void update(const std::vector<uint16_t>& codes,size_t num_codes,size_t num_postings) {
        total_blocks++;
        if(num_postings == block_size) {
            codes_per_block[num_codes]++;
            total_full_blocks++;
        }
        total_postings += num_postings;
        total_codes_u16 += num_codes;
        size_t exceptions = 0;
        for(size_t i=0;i<num_codes;i++) {
            uint16_t code_word = codes[i];
            code_usage[code_word]++;
            if(code_word == 0) {
                total_exceptions_u16++;
                exceptions++;
            }
            if(code_word == 1) {
                total_exceptions_u32++;
                exceptions++;
            }
        }
        if(num_postings == block_size) exceptions_per_block[exceptions]++;
    }

    void print() {
        DS2I_LOG << "(1) dictionary contents:";
        dict.print();

        DS2I_LOG << "(2) code lens distribution:";
        std::vector<uint64_t> len_stats(512);
        for(size_t i=0;i<code_usage.size();i++) {
            auto code_len = dict.size(i);
            len_stats[code_len] += code_usage[i];
        }
        boost::format fmtl("\t len = %1$3d #codes = %2$11d #postings = %3$11d \%%codes = %4$4.2f %%postings = %5$4.2f");
        for(size_t l=0;l<len_stats.size();l++) {
            size_t num_postings = l * code_usage[l];
            double percent_codes = double(len_stats[l])
                / double(total_codes_u16-total_exceptions_u32) * 100;
            if(l <= 1) num_postings = code_usage[l];
            if(len_stats[l] == 0) continue;
            double percent_postings = double(num_postings)
                / double(total_postings) * 100;
            DS2I_LOG << fmtl % l % len_stats[l] % num_postings % percent_codes % percent_postings;
        }


        DS2I_LOG << "(3) codes per block distribution:";
        boost::format fmtd("\t codes = %1$3d blocks = %2$12d bpi = %3$4.3f percent = %4$3.2f");
        for(size_t l=0;l<codes_per_block.size();l++) {
            if(codes_per_block[l] == 0) continue;
            double bpi = double(l*16) / double(block_size);
            double percentage = double(codes_per_block[l]) / double(total_full_blocks) * 100;
            DS2I_LOG << fmtd % l % codes_per_block[l] % bpi % percentage;
        }


        DS2I_LOG << "(4) exceptions per block distribution:";
        boost::format fmt("\t codes = %1$6d blocks = %2$12d  percent = %3$5.2f");
        for(size_t l=0;l<exceptions_per_block.size();l++) {
            if(exceptions_per_block[l] == 0) continue;
            double percentage = double(exceptions_per_block[l]) / double(total_full_blocks) * 100;
            DS2I_LOG << fmt %  l % exceptions_per_block[l] % percentage;
        }

        DS2I_LOG << "(5) overall stats:";

        DS2I_LOG << "\tblock_size = " << block_size;
        DS2I_LOG << "\tencoded blocks = " << total_blocks;
        DS2I_LOG << "\tencoded full blocks = " << total_full_blocks;
        double percent_full = double(total_full_blocks)/double(total_blocks)*100;
        DS2I_LOG << "\tpercent full blocks = " << boost::format("%1$/3.2f") % percent_full;
        DS2I_LOG << "\tencoded postings = " << total_postings;
        DS2I_LOG << "\ttotal u16 codes (inc. exceptions) = " << total_codes_u16;
        DS2I_LOG << "\ttotal u16 exceptions = " << total_exceptions_u16;
        DS2I_LOG << "\ttotal u32 exceptions = " << total_exceptions_u32;
        DS2I_LOG << "\tBPI = " << double(total_codes_u16*16)/double(total_postings);
        auto except_bits = total_exceptions_u16*16+total_exceptions_u32*32;
        DS2I_LOG << "\tEXCEPTION BPI = " << double(except_bits)/double(total_postings);
    }

    std::vector<uint32_t> codes_per_block;
    std::vector<uint32_t> exceptions_per_block;
    std::vector<uint32_t> code_usage;

    size_t total_blocks = 0;
    size_t total_postings = 0;
    size_t total_codes_u16 = 0;
    size_t total_exceptions_u16 = 0;
    size_t total_exceptions_u32 = 0;

    ds2i::dictionary::builder& dict;
    size_t block_size = 0;
};

void encode_lists(ds2i::dictionary::builder& dict,std::string input_basename,dict_type type,size_t block_size)
{
    DS2I_LOG << "encoding lists";

    encoding_stats stats(dict,block_size);

    std::string file_name = input_basename + ".docs";
    if(type == dict_type::freqs) file_name = input_basename + ".freqs";

    binary_collection input(file_name.c_str());

    boost::progress_display progress(input.data_size());
    std::vector<uint32_t> buf(block_size);
    std::vector<uint16_t> output(block_size*4);
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
            stats.update(output,written_codes,block_size);
        }
        if(left) {
            for(size_t j=0;j<left;j++) {
                buf[j] = *itr - prev;
                if(type == dict_type::docs) prev = *itr;
                ++itr;
            }
            size_t written_codes = dint_block::encode(dict,buf.data(),left,output.data());
            stats.update(output,written_codes,left);
        }
        progress += n+1;
    }
    stats.print();
}


template<class dict_constructor_type,class block_stat_type,uint32_t encoding_block_size>
void eval_dict(std::string input_basename,std::string log_prefix)
{
    {
        ds2i::start_logging_to_file(log_prefix + "-docs-" + dict_constructor_type::type());
	    auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::docs);
        auto dict = build_dict<block_stat_type,dict_constructor_type>(block_stats);
        encode_lists(dict,input_basename,dict_type::docs,encoding_block_size);
        stop_logging_to_file();
    }
    {
        ds2i::start_logging_to_file(log_prefix + "-freqs-" + dict_constructor_type::type());
        auto block_stats = create_block_stats<block_stat_type>(input_basename,dict_type::freqs);
        auto dict = build_dict<block_stat_type,dict_constructor_type>(block_stats);
        encode_lists(dict,input_basename,dict_type::freqs,encoding_block_size);
        stop_logging_to_file();
    }
}

int main(int argc, const char** argv) {
    ds2i::init_logging();

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
