#pragma once

#include <boost/progress.hpp>
#include <boost/filesystem.hpp>

#include "binary_blocks_collection.hpp"
#include "dictionary.hpp"
#include "heap.hpp"
#include "hash_utils.hpp"
#include "util.hpp"

#include <unordered_map>
#include <queue>

namespace ds2i {

    template<typename Dictionary, typename Statistics>
    struct decreasing_static_frequencies
    {
        typedef Dictionary dictionary_type;
        typedef Statistics statistics_type;

        static std::string type() {
            return "DSF-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static void build(typename Dictionary::builder& builder,
                          Statistics& block_stats)
        {
            builder.init();

            using block_type = typename Statistics::block_type;
            auto coverage_cmp = [](block_type const& l, block_type const& r) {
                return l.freq > r.freq;
            };

            std::priority_queue<block_type,
                                std::vector<block_type>,
                                decltype(coverage_cmp)> pq(coverage_cmp);
            {
                boost::progress_display progress(block_stats.blocks.size());
                for (const auto& block : block_stats.blocks) {
                    ++progress;
                    if (pq.size() < dictionary_type::num_entries) {
                        pq.push(block);
                    } else {
                        if (coverage_cmp(block,pq.top())) {
                            pq.pop();
                            pq.push(block);
                        }
                    }
                }
            }

            logger() << "(3) add blocks to dict in decreasing freq order";
            std::vector<std::pair<int64_t, block_type>> final_blocks;
            while (!pq.empty()) {
                auto& block = pq.top();
                final_blocks.emplace_back(block.freq,block);
                pq.pop();
            }
            // if estimated freq is the same we prefer longer entries
            auto cmp = [](auto const& l, auto const& r) {
                return l.first > r.first || ((l.first == r.first) && l.second.entry_len > r.second.entry_len);
            };
            std::sort(final_blocks.begin() ,final_blocks.end(), cmp);
            for (auto& dict_entry : final_blocks) {
                auto& block = dict_entry.second;
                builder.append(block.entry,block.entry_len,block.freq);
            }
        }
    };

    template<typename Dictionary, typename Statistics>
    struct prefix_discounted_frequencies
    {
        typedef Dictionary dictionary_type;
        typedef Statistics statistics_type;
        using hash_map = std::unordered_map<uint64_t, uint64_t>;

        static std::string type() {
            return "PDF-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static void build(typename Dictionary::builder& builder,
                          Statistics& block_stats)
        {
            builder.init();

            logger() << "(2) preparing initial estimates";
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
                if(block.entry_len < dictionary_type::max_entry_size) hash_id_map[block.hash] = i;
                freedom[i] = block.freq;
                pq.emplace(freedom[i],i);
            }

            logger() << "(3) find the top-K most covering blocks";
            size_t needed = dictionary_type::num_entries;
            std::vector<int64_t> prefix_ids(dictionary_type::max_entry_size);
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
                // logger() << "\tADD TO DICT with freedom = " << adjust << " - "
                //     << block_stats.block_string(cur_max_id);

                // (c) add freedom of prefixes
                auto num_prefixes = compute_prefix_ids(hash_id_map,block,prefix_ids);
                for(size_t p = num_prefixes;p != 0; p--) {
                    auto p_id = prefix_ids[p-1];
                    adjust = adjust * 2;
                    auto padjust = freedom[p_id];
                    // logger() << "\t\tadjust freedom " << freedom[p_id] << " -> "
                    //     << freedom[p_id] - adjust << " - "
                    //     << block_stats.block_string(p_id);
                    freedom[p_id] = freedom[p_id] - adjust;
                    if(dictionary[p_id] == 1) {
                        // logger() << "\t\tremove from dict " << p_id;
                        dictionary[p_id] = 0;
                        adjust = adjust - padjust;
                        needed = needed + 1;
                    }
                    pq.emplace(freedom[p_id],p_id);
                }
                needed = needed - 1;
            }

            logger() << "(4) add blocks to dict in decreasing freq order";
            using block_type = typename statistics_type::block_type;
            std::vector<std::pair<int64_t, block_type>> final_blocks;
            for (size_t i=0;i<dictionary.size();i++) {
                if (dictionary[i] == 1) {
                    auto& block = block_stats.blocks[i];
                    final_blocks.emplace_back(freedom[i],block);
                }
            }
            std::sort(final_blocks.begin(), final_blocks.end(), [](auto const& l, auto const& r) {
                return l.first > r.first || ((l.first == r.first) && l.second.entry_len > r.second.entry_len);
            });
            for (auto& dict_entry : final_blocks) {
                auto& block = dict_entry.second;
                builder.append(block.entry, block.entry_len, uint64_t(dict_entry.first));
            }
        }

    private:
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
    };

    template<typename Dictionary, typename Statistics>
    struct longest_to_shortest_sweep
    {
        typedef Dictionary dictionary_type;
        typedef Statistics statistics_type;

        static std::string type() {
            return "LSS-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static void build(typename Dictionary::builder& builder,
                          Statistics& block_stats)
        {
            (void) builder;
            (void) block_stats;
        }

    };
}
