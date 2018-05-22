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



    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_builder_PDF
    {
        static std::string type() {
            return "PDF-" + std::to_string(num_entries) + "-" + std::to_string(entry_width);
        }

        template<class block_stat_type>
        static void build(dictionary::builder& builder,block_stat_type& block_stats)
        {
            logger() << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width,type());

            logger() << "(2) preparing initial estimates" << std::endl;
            std::vector<int64_t> freedom(block_stats.blocks.size());
            std::vector<uint8_t> dictionary(block_stats.blocks.size());
            using pqdata_t = std::pair<int64_t,size_t>;
            auto cmp = [](const pqdata_t& left,const pqdata_t& right) { return left.first < right.first;};
            std::priority_queue<pqdata_t, std::vector<pqdata_t>, decltype(cmp) > pq(cmp);
            for(size_t i=0;i<block_stats.blocks.size();i++) {
                const auto& block = block_stats.blocks[i];
                freedom[i] = block.freq;
                pq.emplace(freedom[i],i);
            }

            logger() << "(3) find the top-K most covering blocks" << std::endl;
            size_t needed = num_entries;
            while(needed != 0) {
                // (a) get top item
                auto item = pq.top(); pq.pop();
                auto cur_max_id = item.second;
                if(freedom[cur_max_id] != item.first) {
                    // is the item 'dirty?'
                    continue;
                }
                std::cout << "needed = " << needed << " - dequeue_and_add_to_dict(freedom=" << freedom[cur_max_id] << ",id=" 
                          << cur_max_id << ") - " << block_stats.block_string(cur_max_id) << std::endl;

                // (b) add to dict and adjust freedom of top item
                auto adjust = freedom[cur_max_id];
                dictionary[cur_max_id] = 1;
                auto& block = block_stats.blocks[cur_max_id];

                // (c) add freedom of prefixes
                for(size_t p = block.num_prefixes;p != 0; p--) {
                    auto p_id = block.prefix_ids[p-1];
                    adjust = adjust * 2;
                    auto padjust = freedom[p_id];
                    std::cout << "\tadjust_prefix_freedom(before_freedom=" << freedom[block_id] << ",prefix_id=" 
                          << p_id << ",after_freedom=" << freedom[p_id] - adjust 
                          << ") - " << block_stats.block_string(p_id) << std::endl;
                    freedom[p_id] = freedom[p_id] - adjust;
                    if(dictionary[p_id] == 1) {
                        std::cout << "\tprefix was in dict -> remove and re-add to queue" << std::endl;
                        dictionary[p_id] = 0;
                        adjust = adjust - padjust;
                        needed = needed + 1;
                    }
                    pq.emplace(freedom[p_id],p_id);
                }
                needed = needed - 1;
            }

            logger() << "(3) add blocks to dict in decreasing freq order" << std::endl;
            using btype = typename block_stat_type::block_type;
            std::vector<std::pair<int64_t,btype>> final_blocks;
            for(size_t i=0;i<dictionary.size();i++) {
                if(dictionary[i] == 1) {
                    auto& block = block_stats.blocks[i];
                    final_blocks.emplace_back(-1 * freedom[i],block);
                }
            }
            std::sort(final_blocks.begin(),final_blocks.end());
            for(auto& dict_entry : final_blocks) {
                auto& block = dict_entry.second;
                builder.append(block.entry,block.entry_len,uint64_t(dict_entry.first * -1));
            }
        }

    };


}