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
        static void build(std::ostream& os,dictionary::builder& builder,block_stat_type& block_stats)
        {
            // (1) init dictionary
            os << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width,type());

            // (2) find the top-K most covering blocks
            os << "(2) find the top-K most covering blocks" << std::endl;
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

            os << "(3) add blocks to dict in decreasing freq order" << std::endl;
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
        static void build(std::ostream& os,dictionary::builder& builder,block_stat_type& block_stats)
        {
            os << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width,type());

            os << "(2) preparing initial estimates" << std::endl;
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

            os << "(3) find the top-K most covering blocks" << std::endl;
            size_t needed = num_entries;
            while(needed != 0) {
                // (a) get top item
                auto item = pq.top(); pq.pop();
                auto cur_max_id = item.second;
                if(freedom[cur_max_id] != item.first || dictionary[cur_max_id] == 1) {
                    // is the item 'dirty?'
                    continue;
                }
                os << "needed = " << needed << " - dequeue_and_add_to_dict(freedom=" << freedom[cur_max_id] << ",id=" 
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
                    os << "\tadjust_prefix_freedom(before_freedom=" << freedom[p_id] << ",prefix_id=" 
                          << p_id << ",after_freedom=" << freedom[p_id] - adjust 
                          << ") - " << block_stats.block_string(p_id) << std::endl;
                    freedom[p_id] = freedom[p_id] - adjust;
                    if(dictionary[p_id] == 1) {
                        os << "\tprefix was in dict -> remove and re-add to queue" << std::endl;
                        dictionary[p_id] = 0;
                        adjust = adjust - padjust;
                        needed = needed + 1;
                    }
                    pq.emplace(freedom[p_id],p_id);
                }
                needed = needed - 1;
            }

            os << "(3) add blocks to dict in decreasing freq order" << std::endl;
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




    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_builder_SDF
    {
        static std::string type() {
            return "PDF-" + std::to_string(num_entries) + "-" + std::to_string(entry_width);
        }

        template<class block_type>
        static size_t compute_prefix_ids(const std::unordered_map<uint64_t,uint64_t>& hash_id_map,const block_type& block,std::vector<int64_t>& prefix_ids)
        {
            size_t num_prefixes = 0;
            auto entry_len = block.entry_len;
            for(size_t stride_len=1;stride_len<entry_len;stride_len*=2) {
                for(size_t pos=1;pos<=entry_len;pos+=stride_len) {
                    const uint64_t m = 0x5bd1e9955bd1e995;
                    const int r = 37;
                    uint64_t hash = 0xDEADBEEFDEADBEEF;
                    for(size_t j=0;j<stride_len;j++) {
                        uint64_t key = block.entry[pos-1+j];
                        key *= m;
                        key ^= key >> r;
                        key *= m;
                        hash *= m;
                        hash ^= key;
                    }
                    auto itr = hash_id_map.find(hash);
                    if(itr != hash_id_map.end()) {
                        prefix_ids[num_prefixes++] = itr->second;
                    }
                }
            }
            return num_prefixes;
        }

        template<class block_type>
        static size_t compute_superstring_ids(std::unordered_map<uint64_t,uint64_t>& hash_id_map,const block_type& block,
            std::vector<int64_t>& super_ids,std::vector<int64_t>& super_mult)
        {
            uint64_t buf[entry_width];
            size_t num_super = 0;
            auto entry_len = block.entry_len;
            for(size_t i=0;i<entry_len;i++) {
                buf[i] = block.entry[i%entry_len];
            }

            for(size_t stride_len=entry_len*2;stride_len<=entry_width;stride_len*=2) {
                for(size_t pos=1;pos<=entry_width;pos+=stride_len) {
                    const uint64_t m = 0x5bd1e9955bd1e995;
                    const int r = 37;
                    uint64_t hash = 0xDEADBEEFDEADBEEF;
                    for(size_t j=0;j<stride_len;j++) {
                        uint64_t key = buf[pos-1+j];
                        key *= m;
                        key ^= key >> r;
                        key *= m;
                        hash *= m;
                        hash ^= key;
                    }
                    auto itr = hash_id_map.find(hash);
                    if(itr != hash_id_map.end()) {
                        super_ids[num_super] = itr->second;
                        super_mult[num_super++] = stride_len/entry_len;
                    }
                }
            }
            return num_super;
        }


        template<class block_stat_type>
        static void build(std::ostream& os,dictionary::builder& builder,block_stat_type& block_stats)
        {
            os << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width,type());

            os << "(2) preparing initial estimates" << std::endl;
            std::vector<int64_t> freedom(block_stats.blocks.size());
            std::vector<uint8_t> dictionary(block_stats.blocks.size());
            std::vector<int64_t> prefix_ids(entry_width*2);
            std::vector<int64_t> super_ids(entry_width*2);
            std::vector<int64_t> super_mult(entry_width*2);
            std::unordered_map<uint64_t,uint64_t> hash_id_map;
            using pqdata_t = std::pair<int64_t,size_t>;
            auto cmp = [](const pqdata_t& left,const pqdata_t& right) { return left.first < right.first;};
            std::priority_queue<pqdata_t, std::vector<pqdata_t>, decltype(cmp) > pq(cmp);
            for(size_t i=0;i<block_stats.blocks.size();i++) {
                const auto& block = block_stats.blocks[i];
                hash_id_map[block.hash] = i;
                freedom[i] = block.freq;
                pq.emplace(freedom[i],i);
            }

            os << "(3) find the top-K most covering blocks" << std::endl;
            size_t needed = num_entries;
            while(needed != 0) {
                // (a) get top item
                auto item = pq.top(); pq.pop();
                auto cur_max_id = item.second;
                if(freedom[cur_max_id] != item.first || dictionary[cur_max_id] == 1) {
                    // is the item 'dirty?'
                    continue;
                }
                os << "needed = " << needed << " - dequeue_and_add_to_dict(freedom=" << freedom[cur_max_id] << ",id=" 
                          << cur_max_id << ") - " << block_stats.block_string(cur_max_id) << std::endl;

                // (b) add to dict and adjust freedom of top item
                auto adjust = freedom[cur_max_id];
                dictionary[cur_max_id] = 1;
                auto& block = block_stats.blocks[cur_max_id];

                auto num_prefix = compute_prefix_ids(hash_id_map,block,prefix_ids);

                // (c) adjust freedom of prefixes
                for(size_t p = 0;p < num_prefix;p++) {
                    auto p_id = prefix_ids[p];
                    freedom[p_id] -= adjust;
                    if(dictionary[p_id] == 1) {
                        dictionary[p_id] = 0;
                        needed = needed + 1;
                    }
                    pq.emplace(freedom[p_id],p_id);
                }

                // (c) adjust freedom of superstrings
                auto num_super = compute_superstring_ids(hash_id_map,block,super_ids,super_mult);
                for(size_t s = 0;s < num_super;s++) {
                    auto s_id = super_ids[s];
                    freedom[s_id] -= adjust/super_mult[s];
                    if(dictionary[s_id] == 1) {
                        dictionary[s_id] = 0;
                        needed = needed + 1;
                    }
                    pq.emplace(freedom[s_id],s_id);
                }

                needed = needed - 1;
            }

            os << "(3) add blocks to dict in decreasing freq order" << std::endl;
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