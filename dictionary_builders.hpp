#pragma once

#include <boost/progress.hpp>

#include "binary_blocks_collection.hpp"
#include "dictionary.hpp"
#include "heap.hpp"
#include "hash_utils.hpp"

#include <unordered_map>
#include <queue>

enum class dict_type: char
{
     docs='d',
     freqs='f',
};

namespace ds2i {

    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_builder_DSF
    {
        static std::string type() {
            return "DSF-" + std::to_string(num_entries) + "-" + std::to_string(entry_width);
        }

        template<class block_stat_type>
        static void build(dictionary::builder& builder,block_stat_type& block_stats)
        {
            // (1) init dictionary
            logger() << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width,type());

            // (2) find the top-K most covering blocks
            logger() << "(2) find the top-K most covering blocks" << std::endl;
            using btype = typename block_stats_type::block_type;
            auto coverage_cmp = [](const btype& left,const btype& right) {
                return (left.freq < right.freq);
            };
            std::priority_queue<btype,std::vector<btype>,decltype(coverage_cmp)> pq(coverage_cmp);
            {
                boost::progress_display progress(block_stats.blocks.size());
                for(const auto& block : block_stats.blocks) {
                    ++progress;
                    if(pq.size() < num_entries-1) {
                        pq.push(block);
                    } else {
                        if( coverage_cmp(pq.top(),block) ) {
                            pq.pop();
                            pq.push(block);
                        }
                    }
                }
            }

            logger() << "(3) add blocks to dict in decreasing freq order" << std::endl;
            std::vector<std::pair<int64_t,btype>> final_blocks;
            while(!pq.empty()) {
                auto& block = pq.top();
                final_blocks.emplace_back(-1 * int64_t(block.freq),block);
                pq.pop();
            }
            std::sort(final_blocks.begin(),final_blocks.end());
            for(auto& dict_entry : final_blocks) {
                auto& block = dict_entry.second;
                builder.append(block.entry,block.entry_len,block.freq);
            }
        }

    };

}