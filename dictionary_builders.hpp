#pragma once

#include "binary_blocks_collection.hpp"
#include "dictionary.hpp"
#include "heap.hpp"
#include "hash_utils.hpp"

#include <unordered_map>

enum class dict_type: char
{
     docs='d',
     freqs='f',
};

enum class target_len_seq : char
{
     geometric='g',
     linear='l',
     fibonacci='f',
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

    template<uint32_t max_entry_width>
    struct block_stats_full_stride_geom {

        template<uint32_t width>
        struct block_info {
            size_t freq;
            size_t coverage;
            size_t entry_len;
            uint32_t entry[width];
        };

        using block_type = block_info<max_entry_width>;

        template<class t_list>
        void process_list(t_list& list,dict_type type) {
            thread_local std::vector<uint32_t> buf(max_entry_width);
            size_t n = list.docs.size();
            size_t full_strides = n / max_entry_width;
            logger() << "Procces list of size " << n << " full_strides " << full_strides;

            auto lst_itr = list.docs.begin();
            if(type != dict_type::docs) lst_itr = list.freqs.begin();

            uint32_t prev = 0;
            for(size_t i=0;i<full_strides;i++) {
                // fill buf
                std::string bufstr = "[";
                for(size_t j=0;j<max_entry_width;j++) {
                    buf[j] = *lst_itr - prev;
                    bufstr += "," + std::to_string(buf[j]);
                    if(type == dict_type::docs) prev = *lst_itr;
                    ++lst_itr;
                }

                for()
                logger() << "buf = " << bufstr << "]";
                // compute hashes (modified from murmur! hopefully still good...)
                const uint32_t m = 0x5bd1e995;
                const int r = 24;
                uint32_t hash = 0xDEADBEEF;
                size_t cur_len = 1;
                for(size_t j=1;j<=max_entry_width;j++) {
                    uint32_t key = buf[j-1];
                    key *= m;
                    key ^= key >> r;
                    key *= m;
                    hash *= m;
                    hash ^= key;
                    if(j == cur_len) {
                        logger() << "hash("<<j<<") = " << hash;
                        auto itr = block_map.find(hash);
                        if(itr != block_map.end()) {
                            auto block_idx = itr->second;
                            blocks[block_idx].freq += (max_entry_width >> (j-1));
                            blocks[block_idx].coverage = blocks[block_idx].freq * j;
                        } else {
                            // create a new block
                            block_type new_block;
                            new_block.freq = (max_entry_width >> (j-1));
                            new_block.coverage = new_block.freq * j;
                            new_block.entry_len = j;
                            for(size_t k=0;k<j;k++) {
                                new_block.entry[k] = buf[k];
                            }
                            block_map[hash] = blocks.size();
                            blocks.emplace_back(std::move(new_block));
                        }
                        cur_len *= 2;
                    }
                }
            }
        }

        std::vector<uint32_t> valid_target_lens;
        std::unordered_map<uint32_t,uint64_t> block_map;
        std::vector<block_type> blocks;
    };
}

namespace ds2i {

    template<typename block_stats_type,
             uint32_t num_entries = 65536,
             uint32_t entry_width = 16
             >
    struct dint_dict_builder_smc
    {
        static const size_t MIN_SIZE_THRES = 4096;

        template<typename InputCollection>
        static void build(dictionary::builder& builder,InputCollection const& input,enum dict_type type)
        {
            // (1) create statistics
            block_stats_type stats;
            for (auto const& plist: input) {
                size_t n = plist.docs.size();
                if (n < MIN_SIZE_THRES) continue;
                stats.process_list(plist,type);
            }

            // (2) init dictionary
            builder.init(num_entries, entry_width);

            // (3) find the top-K most covering blocks
            using btype = typename block_stats_type::block_type;
            auto coverage_cmp = [](const btype& left,const btype& right) { return left.coverage < right.coverage;};
            std::priority_queue<btype,std::vector<btype>,decltype(coverage_cmp)> pq(coverage_cmp);
            for(const auto& block : stats.blocks) {
                if(pq.size() < num_entries-1) {
                    pq.push(block);
                } else {
                    if( pq.top().coverage < block.coverage ) {
                        pq.pop();
                        pq.push(block);
                    }
                }
            }

            // (4) add them to the dict builder
            while(!pq.empty()) {
                auto block = pq.top();
                builder.append(block.entry,block.entry_len);
                pq.pop();
            }
        }


    };

}