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

    void print_usage_rbo(double p) {
        // (1) create lists to compare
        struct rbo_stat_type {
            int orig_id;
            size_t freq;
        };
        std::vector<rbo_stat_type> list;
        for(size_t i=dict.special_cases();i<dict.capacity();i++) {
            rbo_stat_type r;
            r.orig_id = i - dict.special_cases() + 1;
            r.freq = code_usage[i];
            list.push_back(r);
        }
        auto freq_cmp = [](const auto& a,const auto& b) {
            if(a.freq == b.freq) return a.orig_id < b.orig_id;
            return a.freq > b.freq;
        };
        std::sort(list.begin(),list.end(),freq_cmp);
        std::vector<int> A;
        std::vector<int> B;
        for(size_t i=0;i<list.size();i++) {
            A.push_back(list[i].orig_id);
            B.push_back(i+1);
        }
        // (2) perform rbo based on ammoffat code
        double weight = 1.0 - p;
        double rbo_min = 0.0;
        uint64_t overlap = 0;
        double contrib = 0;
        std::unordered_set<int> seen_A;
        std::unordered_set<int> seen_B;
        size_t n = 1;
        for(;n<=A.size();n++) {
            if (seen_A.count(A[n-1]) == 1 || seen_B.count(B[n-1]) == 1) {
		        DS2I_LOG << "Duplicate value at line " << n;
	        }

            std::cout << "A["<<n-1<<"] = " << A[n-1] << " B["<<n-1<<"] = " << B[n-1];

            if (A[n-1]==B[n-1]) overlap++;
            else {
                overlap += seen_A.count(A[n-1]);
                overlap += seen_B.count(B[n-1]);
            }
            seen_A.insert(A[n-1]);
            seen_B.insert(B[n-1]);
            contrib = weight*double(overlap)/double(n);
            rbo_min += contrib;

            DS2I_LOG << "n=" << n
                << " weight=" << weight
                << " overlap=" << overlap
                << " contrib=" << contrib
                << " rbo_min=" << rbo_min;

            weight *= p;
        }
        auto max_overlap = overlap;
        auto rbo_max = rbo_min;
        const double EPSILON = 1e-15;
        while (weight>EPSILON) {
            n++;
            contrib = weight*overlap/double(n);
            rbo_min += contrib;
            if (max_overlap==n-1) {
                // both new elements must be the same novel one
                max_overlap += 1;
            } else {
                // two new elements can be assumed, both ones that
                // appeared already in the other list
                max_overlap += 2;
            }
            rbo_max += weight*max_overlap/double(n);
            // prepare for the next pair of imaginary values
		    weight *= p;
        }
        boost::format rbofmt("\t rbo(p = %1$.5f) = %2$.6f + %3$.6f (n=%4$7d, nrows=%5$7d)");
        DS2I_LOG << rbofmt % p % rbo_min % (rbo_max-rbo_min) % A.size() % n;
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

        DS2I_LOG << "(5) RBO stats:";
        auto P = {0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 0.999, 0.999};
        for(auto p : P) {
            print_usage_rbo(p);
        }

        DS2I_LOG << "(6) overall stats:";

        DS2I_LOG << "\tblock_size = " << block_size;
        DS2I_LOG << "\tencoded blocks = " << total_blocks;
        DS2I_LOG << "\tencoded full blocks = " << total_full_blocks;
        double percent_full = double(total_full_blocks)/double(total_blocks)*100;
        DS2I_LOG << "\tpercent full blocks = " << boost::format("%1$3.2f") % percent_full;
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
    size_t total_full_blocks = 0;

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
