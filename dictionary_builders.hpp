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

// namespace ds2i {

//     template<uint32_t num_entries = 65536,
//              uint32_t entry_width = 16,
//              target_len_seq target_lens>
//     struct dint_dict_builder_bbc
//     {
//         static double cost(uint32_t block_size, uint32_t block_frequency) {
//             return block_frequency * (32.0 * block_size - 16.0);
//         }

//         // Giulio: we can exploit the property that the distribution
//         // is highly skewed to speed up the computation.
//         // Instead of selecting top_k entries from EACH statistics file,
//         // we can use an adaptive distribution, i.e., select more entries
//         // if they appear a lot.
//         static const uint32_t top_k = 50000;

//         // <block, index into the frequency array>
//         typedef std::pair<std::vector<uint32_t>, uint32_t> entry_type;

//         static double bpi(uint32_t block_size, uint32_t block_frequency, uint64_t total_integers) {
//             return cost(block_size, block_frequency) / total_integers;
//         };

//         struct bpi_comparator {
//             bpi_comparator()
//                 : m_total_integers(0)
//                 , m_frequencies(nullptr)
//             {}

//             bpi_comparator(std::vector<uint64_t> const* frequencies, uint64_t total_integers)
//                 : m_total_integers(total_integers)
//                 , m_frequencies(frequencies)
//             {}

//             bool operator()(entry_type const& x, entry_type const& y) {
//                 return bpi(x.first.size(), m_frequencies->operator[](x.second), m_total_integers)
//                      < bpi(y.first.size(), m_frequencies->operator[](y.second), m_total_integers);
//             }

//         private:
//             uint64_t m_total_integers;
//             std::vector<uint64_t> const* m_frequencies;
//         };

//         template<typename dict_type>
//         static void build(dictionary::builder& builder, double const* percentages,
//                           uint64_t total_integers, std::string prefix_name)
//         {
//             builder.init(num_entries, entry_width);

//             // assume that everything NOT covered by the dictionary is left uncompressed
//             double final_bpi = 32.0;

//             // 1. select the top-k entries from each statistics file
//             // and put:
//             //   - the blocks into a max heap;
//             //   - the frequencies into a vector;
//             // 2. create an hash table mapping blocks to indices of such vector
//             // (in order to access efficiently the frequencies)

//             std::vector<uint64_t> blocks_frequencies;
//             blocks_frequencies.reserve(top_k * 5); // top-k entries for each different block_size

//             bpi_comparator bc(&blocks_frequencies, total_integers);

//             typedef heap<entry_type, bpi_comparator> max_heap_type;
//             std::vector<max_heap_type> max_heaps(5, max_heap_type(top_k, bc));

//             // <hash of block, index into frequency array>
//             std::unordered_map<uint64_t, uint32_t> map;

//             for (uint32_t block_size = 16, i = 0; block_size != 0; block_size /= 2, ++i)
//             {
//                 // logger() << "loading " << top_k << " most frequent blocks of size " << block_size << "..." << std::endl;
//                 std::string collection_name("./" + prefix_name + ".blocks_stats." + std::to_string(block_size) + ".bin");
//                 binary_blocks_collection input(collection_name.c_str());
//                 auto begin = input.begin();
//                 for (uint32_t k = 0; k < top_k; ++k, ++begin)
//                 {
//                     uint32_t index = blocks_frequencies.size();
//                     entry_type entry(std::vector<uint32_t>(), index);
//                     auto const& block = *begin;

//                     // if (k < 10) std::cout << "block.freq() = " << block.freq() << std::endl;

//                     entry.first.reserve(block.size());
//                     for (uint32_t x: block) {
//                         entry.first.push_back(x);
//                     }

//                     max_heaps[i].push_back(std::move(entry));
//                     blocks_frequencies.push_back(block.freq());

//                     uint8_t const* b = reinterpret_cast<uint8_t const*>(entry.first.data());
//                     uint8_t const* e = b + block.size() * sizeof(uint32_t);
//                     uint64_t hash = hash_bytes64(byte_range(b, e));
//                     map[hash] = index;
//                 }
//                 max_heaps[i].make();
//             }


//             // for (int i = 0; i < 10; ++i) {
//             //     auto const& best = max_heap.top();
//             //     auto const& best_block = best.first;
//             //     std::cout << "block.size() " << best_block.size() << std::endl;
//             //     std::cout << "block.freq() " << blocks_frequencies[best.second] << std::endl;
//             //     max_heap.pop();
//             // }

//             // 3. keep adding entries to the dictionary until it fills up:
//             //   - at each step we remove the entry from the heap (in O(1))
//             //     that gives the highest reduction in used bits per integer (bpi);
//             //   - we decrease the frequency of all its sub-entries (if found) and
//             //     call heapify() to restore the heap condition
//             // uint64_t covered_integers = 0;

//             double total_coverage = 0.0;
//             uint32_t added_entries = 0;

//             for (int i = 0; i < 5 and not builder.full(); ++i) {
//                 double p = 0.0;
//                 while (p < percentages[i])
//                 {
//                     if (max_heaps[i].empty()) break;
//                     auto const& best = max_heaps[i].top();
//                     auto const& best_block = best.first;

//                     if (not builder.append(best_block.data(), best_block.size())) {
//                         break;
//                     }

//                     ++added_entries;
//                     uint64_t best_block_freq = blocks_frequencies[best.second];
//                     p += best_block_freq * best_block.size() * 100.0 / total_integers;
//                     double cost_saving = bpi(best_block.size(), best_block_freq, total_integers);
//                     final_bpi -= cost_saving;

//                     if (added_entries % 200 == 0) {
//                         std::cout << "added_entries " << added_entries << "/65536" << std::endl;
//                         std::cout << "p = " << p << "/" << percentages[i] << std::endl;
//                         // std::cout << "saving " << cost_saving << " bits x int" << std::endl;

//                         // std::cout << "adding a target of size " << best.first.size() << " to the dictionary" << std::endl;

//                         // for (auto x: best_block) {
//                         //     std::cout << x << " ";
//                         // }
//                         // std::cout << std::endl;

//                         std::cout << "current bits x integer: " << final_bpi << std::endl;
//                         std::cout << "covering " << total_coverage + p << "% of integers" << std::endl;
//                     }

//                     // logger() << "decreasing freq. of sub-blocks..." << std::endl;
//                     for (uint32_t block_size = best_block.size() / 2, j = i + 1; block_size != 0; block_size /= 2, ++j) {
//                         // std::cout << "decreasing freqs of sub-blocks of size " << block_size << std::endl;
//                         for (uint32_t begin = 0; begin < best_block.size(); begin += block_size) {

//                             // uint32_t end = std::min<uint32_t>(begin + block_size, best_block.size());
//                             // std::cout << "block[" << begin << ", " << end << ") = ";
//                             // for (uint32_t kk = begin; kk != end; ++kk) {
//                             //     std::cout << best_block[kk] << " ";
//                             // }
//                             // std::cout << std::endl;

//                             uint8_t const* b = reinterpret_cast<uint8_t const*>(&best_block[begin]);
//                             uint8_t const* e = b + std::min<uint64_t>(block_size, best_block.size() - begin) * sizeof(uint32_t);
//                             uint64_t hash = hash_bytes64(byte_range(b, e));
//                             uint32_t index = map[hash];

//                             // std::cout << "decreasing " << blocks_frequencies[index] << " by " << best_block_freq << std::endl;
//                             blocks_frequencies[index] -= best_block_freq;
//                             assert(blocks_frequencies[index] >= 0);
//                         }
//                         // Giulio: if we keep the positions of the elements in the heap,
//                         // we can use sink(position): O(log n) vs O(n)
//                         max_heaps[j].make();
//                     }

//                     max_heaps[i].pop();
//                 }

//                 total_coverage += p;
//             }

//             logger() << "using " << final_bpi << " bits x integer" << std::endl;
//             logger() << "covering " << total_coverage << "% of integers" << std::endl;
//         }

//     };
// }

namespace ds2i {

    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_builder_sbc
    {
        static const size_t MIN_SIZE_THRES = 4096;

        template<class block_stat_type>
        static void build(dictionary::builder& builder,block_stat_type& block_stats)
        {
            // (1) init dictionary
            logger() << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width);

            // (2) find the top-K most covering blocks
            logger() << "(2) compute initial savings" << std::endl;
            std::vector<uint64_t> block_ids;
            std::vector<uint64_t> savings;
            std::vector<uint64_t> F;
            std::unordered_map<uint64_t,uint64_t> dictionary;
            std::unordered_map<uint64_t,uint64_t> bid_map;
            std::unordered_map<uint64_t,uint64_t> rbid_map;
            for(size_t i=0;i<block_stats.blocks.size();i++) {
                const auto& b = block_stats.blocks[i];
                auto max_savings = (3*b.entry_len -1) * b.freq;
                if(max_savings > 32) {
                    F.push_back(b.freq);
                    savings.push_back(max_savings);
                    bid_map[i] = block_ids.size();
                    rbid_map[block_ids.size()] = i;
                        .push_back(i);
                }
            }
            logger() << "Considering " << block_ids.size() << " blocks out of " << block_stats.blocks.size() << std::endl;

            logger() << "(3) performing coverage with adjusted counts" << std::endl;
            {
                boost::progress_display progress(num_entries);
                size_t needed = num_entries;
                size_t next = num_entries-1;
                while(needed != 0) {
                    // (1) find next best block to add
                    auto max_savings_mapped_bid = std::distance(savings.begin(),std::max_element(savings.begin(),savings.end()));
                    dictionary[max_savings_mapped_bid] = savings[max_savings_mapped_bid];
                    savings[max_savings_mapped_bid] = 0;
                    auto max_freq = F[max_savings_mapped_bid];
                    // (2) find the prefixes and adjust freqs
                    auto orig_max_id = rbid_map[max_savings_mapped_bid];
                    auto& max_block = block_stats.blocks[orig_max_id];
                    for(size_t i=0;i<max_block.num_prefixes;i++) {
                        auto prefix_id = max_block.prefix_ids[i];
                        auto map_itr = bid_map.find(prefix_id);
                        if(map_itr != bid_map.end()) {
                            auto mapped_prefix_id = map_itr->second;
                            auto& pb = block_stats.blocks[prefix_id];
                            F[mapped_prefix_id] -= max_freq;
                            savings[mapped_prefix_id] = (3*pb.entry_len -1) * F[mapped_prefix_id];

                            // (3) savings decreased. this guy has to try again!
                            auto ditr = dictionary.find(prefix_id);
                            if(ditr != dictionary.end()) {
                                dictionary.erase(ditr);
                                needed++;
                            }
                        }

                    }
                    needed--;
                    if(next == needed) {
                        ++progress;
                        next--;
                    }
                }
            }

            logger() << "(3) add blocks to dict" << std::endl;
            for(auto& dict_entry : dictionary) {
                auto mapped_block_id = dict_entry.first;
                block_id = rbid_map[mapped_block_id];
                auto savings = dict_entry.second;
                auto freq = F[mapped_block_id];
                auto& block = block_stats.blocks[block_id];
                builder.append(block.entry,block.entry_len,freq,savings);
            }
        }


    };

    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_builder_smc
    {
        static const size_t MIN_SIZE_THRES = 4096;

        template<class block_stat_type>
        static void build(dictionary::builder& builder,block_stat_type& block_stats)
        {
            // (1) init dictionary
            logger() << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width);

            // (2) find the top-K most covering blocks
            logger() << "(2) find the top-K most covering blocks" << std::endl;
            using btype = typename block_stats_type::block_type;
            auto coverage_cmp = [](const btype& left,const btype& right) { return left.coverage < right.coverage;};
            std::priority_queue<btype,std::vector<btype>,decltype(coverage_cmp)> pq(coverage_cmp);
            {
                boost::progress_display progress(block_stats.blocks.size());
                for(const auto& block : block_stats.blocks) {
                    ++progress;
                    if(pq.size() < num_entries-1) {
                        pq.push(block);
                    } else {
                        if( pq.top().coverage < block.coverage ) {
                            pq.pop();
                            pq.push(block);
                        }
                    }
                }
            }

            logger() << "(3) add blocks to dict" << std::endl;
            while(!pq.empty()) {
                auto block = pq.top();
                builder.append(block.entry,block.entry_len,block.freq,block.coverage);
                pq.pop();
            }
        }


    };



    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_builder_smc3l
    {
        static const size_t MIN_SIZE_THRES = 4096;

        template<class block_stat_type>
        static void build(dictionary::builder& builder,block_stat_type& block_stats)
        {
            // (1) init dictionary
            logger() << "(1) init dictionary" << std::endl;
            builder.init(num_entries, entry_width);

            // (2) find the top-K most covering blocks
            logger() << "(2) find the top-K most covering blocks" << std::endl;
            using btype = typename block_stats_type::block_type;
            auto coverage_cmp = [](const btype& left,const btype& right) {
                return (left.freq * (3*left.entry_len-1)) < (right.freq * (3*right.entry_len-1));
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

            logger() << "(3) add blocks to dict" << std::endl;
            while(!pq.empty()) {
                auto block = pq.top();
                builder.append(block.entry,block.entry_len,block.freq,block.coverage);
                pq.pop();
            }
        }


    };

}