#pragma once

#include <boost/progress.hpp>

#include "binary_blocks_collection.hpp"
#include "dictionary.hpp"
#include "heap.hpp"
#include "hash_utils.hpp"
#include "util.hpp"

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
    struct dint_dict_type_DSF
    {
        static std::string type() {
            return "DSF-" + std::to_string(num_entries) + "-" + std::to_string(entry_width);
        }

        template<class block_stat_type>
        static void build(dictionary::builder& builder,block_stat_type& block_stats)
        {
            // (1) init dictionary
            DS2I_LOG << "(1) init dictionary";
            builder.init(num_entries, entry_width,type());

            // (2) find the top-K most covering blocks
            DS2I_LOG << "(2) find the top-K most covering blocks";
            using btype = typename block_stats_type::block_type;
            auto coverage_cmp = [](const btype& left,const btype& right) {
                return (left.freq > right.freq);
            };
            std::priority_queue<btype,std::vector<btype>,decltype(coverage_cmp)> pq(coverage_cmp);
            {
                boost::progress_display progress(block_stats.blocks.size());
                for(const auto& block : block_stats.blocks) {
                    ++progress;
                    if(pq.size() < num_entries) {
                        pq.push(block);
                    } else {
                        if( coverage_cmp(block,pq.top()) ) {
                            pq.pop();
                            pq.push(block);
                        }
                    }
                }
            }

            DS2I_LOG << "(3) add blocks to dict in decreasing freq order";
            std::vector<std::pair<int64_t,btype>> final_blocks;
            while(!pq.empty()) {
                auto& block = pq.top();
                final_blocks.emplace_back(block.freq,block);
                pq.pop();
            }
            // if estimated freq is the same we prefer longer entries
            auto cmp = [](const auto& a,const auto& b) {
                return a.first > b.first || ((a.first == b.first) && a.second.entry_len > b.second.entry_len);
            };
            std::sort(final_blocks.begin(),final_blocks.end(),cmp);
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
    struct dint_dict_type_PDF
    {
        using hash_map = std::unordered_map<uint64_t,uint64_t>;

        static std::string type() {
            return "PDF-" + std::to_string(num_entries) + "-" + std::to_string(entry_width);
        }

        template<class block_type>
        static size_t compute_prefix_ids(const hash_map& hash_id_map,const block_type& block,std::vector<int64_t>& prefix_ids)
        {
            size_t num_prefixes = 0;
    		for(size_t size_u32=1;size_u32<block.entry_len;size_u32*=2) {
                auto hash = hash_u32(block.entry,size_u32);
                auto itr = hash_id_map.find(hash);
                prefix_ids[num_prefixes++] = itr->second;
    		}
            return num_prefixes;
        }

        template<class block_stat_type>
        static void build(dictionary::builder& builder,block_stat_type& block_stats)
        {
            DS2I_LOG << "(1) init dictionary";
            builder.init(num_entries, entry_width,type());

            DS2I_LOG << "(2) preparing initial estimates";
            std::vector<int64_t> freedom(block_stats.blocks.size());
            std::vector<uint8_t> dictionary(block_stats.blocks.size());
            using pqdata_t = std::pair<int64_t,size_t>;
            auto cmp = [](const pqdata_t& left,const pqdata_t& right) { return left.first < right.first;};
            std::priority_queue<pqdata_t, std::vector<pqdata_t>, decltype(cmp) > pq(cmp);
            hash_map hash_id_map;
            hash_id_map.max_load_factor(0.05);
            hash_id_map.reserve(block_stats.blocks.size());
            for(size_t i=0;i<block_stats.blocks.size();i++) {
                const auto& block = block_stats.blocks[i];
                if(block.entry_len < entry_width) hash_id_map[block.hash] = i;
                freedom[i] = block.freq;
                pq.emplace(freedom[i],i);
            }

            DS2I_LOG << "(3) find the top-K most covering blocks";
            size_t needed = num_entries;
            std::vector<int64_t> prefix_ids(entry_width);
            while(needed != 0) {
                // (a) get top item
                auto item = pq.top(); pq.pop();
                auto cur_max_id = item.second;
                if(freedom[cur_max_id] != item.first || dictionary[cur_max_id] == 1) {
                    // is the item 'dirty?'
                    continue;
                }

                // (b) add to dict and adjust freedom of top item
                auto adjust = freedom[cur_max_id];
                dictionary[cur_max_id] = 1;
                auto& block = block_stats.blocks[cur_max_id];
                // DS2I_LOG << "\tADD TO DICT with freedom = " << adjust << " - " 
                //     << block_stats.block_string(cur_max_id);

                // (c) add freedom of prefixes
                auto num_prefixes = compute_prefix_ids(hash_id_map,block,prefix_ids);
                for(size_t p = num_prefixes;p != 0; p--) {
                    auto p_id = prefix_ids[p-1];
                    adjust = adjust * 2;
                    auto padjust = freedom[p_id];
                    // DS2I_LOG << "\t\tadjust freedom " << freedom[p_id] << " -> " 
                    //     << freedom[p_id] - adjust << " - "
                    //     << block_stats.block_string(p_id);
                    freedom[p_id] = freedom[p_id] - adjust;
                    if(dictionary[p_id] == 1) {
                        // DS2I_LOG << "\t\tremove from dict " << p_id;
                        dictionary[p_id] = 0;
                        adjust = adjust - padjust;
                        needed = needed + 1;
                    }
                    pq.emplace(freedom[p_id],p_id);
                }
                needed = needed - 1;
            }

            DS2I_LOG << "(4) add blocks to dict in decreasing freq order";
            using btype = typename block_stat_type::block_type;
            std::vector<std::pair<int64_t,btype>> final_blocks;
            for(size_t i=0;i<dictionary.size();i++) {
                if(dictionary[i] == 1) {
                    auto& block = block_stats.blocks[i];
                    final_blocks.emplace_back(freedom[i],block);
                }
            }
            auto final_cmp = [](const auto& a,const auto& b) {
                return a.first > b.first || ((a.first == b.first) && a.second.entry_len > b.second.entry_len);
            };
            std::sort(final_blocks.begin(),final_blocks.end(),final_cmp);
            for(auto& dict_entry : final_blocks) {
                auto& block = dict_entry.second;
                builder.append(block.entry,block.entry_len,uint64_t(dict_entry.first));
            }
        }
    };




    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_type_SDF
    {
        static std::string type() {
            return "SDF-" + std::to_string(num_entries) + "-" + std::to_string(entry_width);
        }

        template<class block_stat_type>
        static void build(std::ostream& os,dictionary::builder& builder,block_stat_type& block_stats)
        {
            // DS2I_LOG << "(1) init dictionary";
            // builder.init(num_entries, entry_width,type());

            // DS2I_LOG << "(2) preparing initial estimates";
            // std::vector<int64_t> freedom(block_stats.blocks.size());
            // std::vector<uint8_t> dictionary(block_stats.blocks.size());
            // using pqdata_t = std::pair<int64_t,size_t>;
            // auto cmp = [](const pqdata_t& left,const pqdata_t& right) { return left.first < right.first;};
            // std::priority_queue<pqdata_t, std::vector<pqdata_t>, decltype(cmp) > pq(cmp);
            // hash_map hash_id_map;
            // hash_id_map.max_load_factor(0.05);
            // hash_id_map.reserve(block_stats.blocks.size());
            // for(size_t i=0;i<block_stats.blocks.size();i++) {
            //     const auto& block = block_stats.blocks[i];
            //     if(block.entry_len < entry_width) hash_id_map[block.hash] = i;
            //     freedom[i] = block.freq;
            //     pq.emplace(freedom[i],i);
            // }

            // DS2I_LOG << "(3) find the top-K most covering blocks";
            // size_t needed = num_entries;
            // std::vector<int64_t> prefix_ids(entry_width);
            // while(needed != 0) {
            //     // (a) get top item
            //     auto item = pq.top(); pq.pop();
            //     auto cur_max_id = item.second;
            //     if(freedom[cur_max_id] != item.first || dictionary[cur_max_id] == 1) {
            //         // is the item 'dirty?'
            //         continue;
            //     }

            //     // (b) add to dict and adjust freedom of top item
            //     auto adjust = freedom[cur_max_id];
            //     dictionary[cur_max_id] = 1;
            //     auto& block = block_stats.blocks[cur_max_id];
            //     // DS2I_LOG << "\tADD TO DICT with freedom = " << adjust << " - " 
            //     //     << block_stats.block_string(cur_max_id);

            //     // (c) add freedom of prefixes
            //     auto num_prefixes = compute_prefix_ids(hash_id_map,block,prefix_ids);
            //     for(size_t p = num_prefixes;p != 0; p--) {
            //         auto p_id = prefix_ids[p-1];
            //         adjust = adjust * 2;
            //         auto padjust = freedom[p_id];
            //         // DS2I_LOG << "\t\tadjust freedom " << freedom[p_id] << " -> " 
            //         //     << freedom[p_id] - adjust << " - "
            //         //     << block_stats.block_string(p_id);
            //         freedom[p_id] = freedom[p_id] - adjust;
            //         if(dictionary[p_id] == 1) {
            //             // DS2I_LOG << "\t\tremove from dict " << p_id;
            //             dictionary[p_id] = 0;
            //             adjust = adjust - padjust;
            //             needed = needed + 1;
            //         }
            //         pq.emplace(freedom[p_id],p_id);
            //     }
            //     needed = needed - 1;
            // }

            // DS2I_LOG << "(4) add blocks to dict in decreasing freq order";
            // using btype = typename block_stat_type::block_type;
            // std::vector<std::pair<int64_t,btype>> final_blocks;
            // for(size_t i=0;i<dictionary.size();i++) {
            //     if(dictionary[i] == 1) {
            //         auto& block = block_stats.blocks[i];
            //         final_blocks.emplace_back(freedom[i],block);
            //     }
            // }
            // auto final_cmp = [](const auto& a,const auto& b) {
            //     return a.first > b.first || ((a.first == b.first) && a.second.entry_len > b.second.entry_len);
            // };
            // std::sort(final_blocks.begin(),final_blocks.end(),final_cmp);
            // for(auto& dict_entry : final_blocks) {
            //     auto& block = dict_entry.second;
            //     builder.append(block.entry,block.entry_len,uint64_t(dict_entry.first));
            // }
        }

    };



}