#pragma once

#include <boost/progress.hpp>
#include <boost/filesystem.hpp>

#include <unordered_map>
#include <queue>

#include "dint_configuration.hpp"
#include "statistics_collectors.hpp"
#include "binary_blocks_collection.hpp"
#include "dictionary.hpp"
#include "heap.hpp"
#include "hash_utils.hpp"
#include "util.hpp"
#include "model_build_utils.hpp"

namespace ds2i {

    template<typename Block>
    struct freq_sorter {
        bool operator()(Block const& l, Block const& r) {
            return l.freq > r.freq;
        }
    };

    template<typename Block>
    struct length_freq_sorter {
        bool operator()(Block const& l, Block const& r) {
            if (l.size == r.size) {
                return l.freq > r.freq;
            }
            return l.size > r.size;
        }
    };

    template<typename Dictionary, typename Statistics>
    struct decreasing_static_frequencies
    {
        typedef Dictionary dictionary_type;
        typedef Statistics statistics_type;
        typedef freq_sorter<typename statistics_type::block_type> sorter_type;

        static std::string type() {
            return "DSF-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static void build(typename dictionary_type::builder& builder,
                          statistics_type& stats)
        {
            builder.init(stats.total_integers);

            using block_type = typename statistics_type::block_type;
            auto coverage_cmp = [](block_type const& l, block_type const& r) {
                return l.freq > r.freq;
            };

            logger() << "(1) finding the most covering blocks";
            std::priority_queue<block_type,
                                std::vector<block_type>,
                                decltype(coverage_cmp)> pq(coverage_cmp);
            {
                boost::progress_display progress(stats.blocks.size());
                for (auto const& block : stats.blocks) {
                    ++progress;
                    if (pq.size() < dictionary_type::num_entries) {
                        pq.push(block);
                    } else {
                        if (coverage_cmp(block, pq.top())) {
                            pq.pop();
                            pq.push(block);
                        }
                    }
                }
            }

            logger() << "(2) adding entries to the dictionary";
            std::vector<std::pair<int64_t, block_type>> final_blocks;
            while (!pq.empty()) {
                auto& block = pq.top();
                final_blocks.emplace_back(block.freq,block);
                pq.pop();
            }

            std::sort(final_blocks.begin(), final_blocks.end(),
            [](auto const& l, auto const& r) { // if estimated freq is the same we prefer longer entries
                return (l.first  > r.first) or ((l.first == r.first) and l.second.size > r.second.size);
            });

            for (auto& dict_entry : final_blocks) {
                auto& block = dict_entry.second;
                builder.append(block.entry, block.size, block.freq);
            }
            logger() << "DONE" << std::endl;
        }
    };

    template<typename Dictionary, typename Statistics>
    struct prefix_discounted_frequencies
    {
        typedef Dictionary dictionary_type;
        typedef Statistics statistics_type;
        typedef freq_sorter<typename statistics_type::block_type> sorter_type;

        using hash_map = std::unordered_map<uint64_t, uint64_t>;

        static std::string type() {
            return "PDF-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static void build(typename dictionary_type::builder& builder,
                          statistics_type& stats)
        {
            builder.init(stats.total_integers);

            logger() << "(1) preparing initial estimates";
            std::vector<int64_t> freedom(stats.blocks.size());
            std::vector<uint8_t> dictionary(stats.blocks.size());
            using pqdata_t = std::pair<int64_t,size_t>;
            auto cmp = [](pqdata_t const& l, pqdata_t const& r) {
                return l.first < r.first;
            };
            std::priority_queue<pqdata_t, std::vector<pqdata_t>, decltype(cmp) > pq(cmp);
            hash_map hash_id_map;
            hash_id_map.max_load_factor(0.05);
            hash_id_map.reserve(stats.blocks.size());
            for (size_t i = 0; i < stats.blocks.size(); ++i) {
                auto const& block = stats.blocks[i];
                if (block.size < dictionary_type::max_entry_size) {
                    hash_id_map[block.hash] = i;
                }
                freedom[i] = block.freq;
                pq.emplace(freedom[i], i);
            }

            logger() << "(2) finding the most covering blocks";
            size_t needed = dictionary_type::num_entries;
            std::vector<int64_t> prefix_ids(dictionary_type::max_entry_size);
            while (needed != 0) {
                // (a) get top item
                auto item = pq.top(); pq.pop();
                auto cur_max_id = item.second;
                if (freedom[cur_max_id] != item.first or dictionary[cur_max_id] == 1) {
                    // is the item 'dirty?'
                    continue;
                }

                // (b) add to dict and adjust freedom of top item
                auto adjust = freedom[cur_max_id];
                dictionary[cur_max_id] = 1;
                auto& block = stats.blocks[cur_max_id];
                // logger() << "\tADD TO DICT with freedom = " << adjust << " - "
                //     << block_stats.block_string(cur_max_id);

                // (c) add freedom of prefixes
                auto num_prefixes = compute_prefix_ids(hash_id_map, block, prefix_ids);
                for (size_t p = num_prefixes; p != 0; --p) {
                    auto p_id = prefix_ids[p-1];
                    adjust = adjust * 2;
                    auto padjust = freedom[p_id];
                    // logger() << "\t\tadjust freedom " << freedom[p_id] << " -> "
                    //     << freedom[p_id] - adjust << " - "
                    //     << block_stats.block_string(p_id);
                    freedom[p_id] = freedom[p_id] - adjust;
                    if (dictionary[p_id] == 1) {
                        // logger() << "\t\tremove from dict " << p_id;
                        dictionary[p_id] = 0;
                        adjust = adjust - padjust;
                        needed = needed + 1;
                    }
                    pq.emplace(freedom[p_id],p_id);
                }
                needed = needed - 1;
            }

            logger() << "(3) adding entries to the dictionary";
            using block_type = typename statistics_type::block_type;
            std::vector<std::pair<int64_t, block_type>> final_blocks;
            for (size_t i = 0; i < dictionary.size(); ++i) {
                if (dictionary[i] == 1) {
                    auto& block = stats.blocks[i];
                    final_blocks.emplace_back(freedom[i], block);
                }
            }
            std::sort(final_blocks.begin(), final_blocks.end(), [](auto const& l, auto const& r) {
                return (l.first  > r.first) or
                      ((l.first == r.first) and
                        l.second.size > r.second.size);
            });

            for (auto& dict_entry : final_blocks) {
                auto& block = dict_entry.second;
                builder.append(block.entry, block.size, uint64_t(dict_entry.first));
            }
        }

    private:
        template<typename block_type>
        static size_t compute_prefix_ids(hash_map const& hash_id_map,
                                         block_type const& block,
                                         std::vector<int64_t>& prefix_ids)
        {
            size_t num_prefixes = 0;
            for (size_t size_u32=1; size_u32 < block.size; size_u32 *= 2) {
                auto hash = hash_bytes64(block.entry, size_u32);
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
        typedef length_freq_sorter<typename statistics_type::block_type> sorter_type;

        static std::string type() {
            return "LSS-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static void build(typename dictionary_type::builder& builder,
                          statistics_type& stats)
        {
            builder.init(stats.total_integers);

            // <hash of block, unique_id>
            std::unordered_map<uint64_t,
                               uint32_t // assume less than 2^32 candidates
                               > map;

            std::vector<uint64_t> id_lowerbounds(constants::max_fractal_steps + 1, 0);
            uint32_t id = 0;
            uint32_t curr_block_size = constants::max_block_length;
            uint32_t k = 0;
            for (auto const& block: stats.blocks)
            {
                if (block.size == curr_block_size / 2) {
                    id_lowerbounds[++k] = id;
                    curr_block_size /= 2;
                }

                map[block.hash] = id;
                ++id;
            }

            assert(id == stats.blocks.size());
            id_lowerbounds[constants::max_fractal_steps] = id;

            curr_block_size = constants::max_block_length;
            k = 0;
            for (uint32_t i = 0; i < stats.blocks.size(); ++i)
            {
                if (builder.full()) break;

                auto const& block = stats.blocks[i];
                builder.append(block.entry, block.size, block.freq);

                // decrease frequencies of smaller blocks (if any)
                for (uint32_t block_size = block.size / 2; block_size != 0; block_size /= 2) {
                    for (uint32_t begin = 0; begin < block.size; begin += block_size) {
                        uint8_t const* b = reinterpret_cast<uint8_t const*>(&(block.entry[begin]));
                        uint8_t const* e = b + std::min<uint64_t>(block_size, block.size - begin) * sizeof(uint32_t);
                        assert(std::min<uint64_t>(block_size, block.size - begin) == block_size);
                        uint64_t hash = hash_bytes64(byte_range(b, e));
                        auto it = map.find(hash);
                        if (it != map.end()) {
                            uint32_t id = map[hash];
                            if (stats.blocks[id].freq >= block.freq) {
                                stats.blocks[id].freq -= block.freq;
                            }
                        }
                    }
                }

                if (block.size == curr_block_size / 2) {
                    logger() << "covering " << builder.coverage() << "% of integers "
                             << "with entries of size " << block.size << std::endl;
                    curr_block_size /= 2;
                    ++k;
                    // sort sub-blocks after decreasing of frequencies
                    sorter_type sorter;
                    std::sort(
                        stats.blocks.begin() + id_lowerbounds[k],
                        stats.blocks.begin() + id_lowerbounds[k + 1],
                        sorter
                    );
                }
            }
        }
    };
}
