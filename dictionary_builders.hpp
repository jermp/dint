#pragma once

#include "dictionary.hpp"
#include "heap.hpp"

namespace ds2i {

    namespace util {
        double cost(uint32_t block_size, uint32_t block_frequency) {
            return block_frequency * (32.0 * block_size - 16.0);
        }
    }

    template<uint32_t num_entries = 65536,
             uint32_t entry_width = 16>
    struct dictionary_builder
    {
        // Giulio: we can exploit the property that the distribution
        // is highly skewed to speed up the computation.
        // Instead of selecting top_k entries from EACH statistics file,
        // we can use an adaptive distribution, i.e., select more entries
        // if they appear a lot.
        static const uint32_t top_k = 500000;

        static void build(dictionary& dict, uint64_t total_num_integers, std::string prefix_name)
        {
            auto bpi = [&](uint32_t block_size, uint32_t block_frequency) {
                return util::cost(block_size, block_frequency) / total_num_integers;
            };

            dictionary::builder builder(num_entries, entry_width);

            // assume that everything NOT covered by the dictionary is left uncompressed
            double final_bpi = 32.0;

            // 1. select the top-1M entries from each statistics file
            // and put:
            //   - the blocks into a max heap;
            //   - the frequencies into a vector;
            // 2. create an hash table mapping blocks to indices of such vector
            // (in order to access efficiently the frequencies)

            std::vector<uint64_t> blocks_frequencies;
            blocks_frequencies.reserve(top_k * 5); // top_k entries for each different block_size

            // <block, index into the frequency array>
            typedef std::pair<std::vector<uint32_t>, uint32_t> entry_type;

            auto compare_bpi = [&](entry_type const& x, entry_type const& y) {
                return bpi(x.first.size(), frequencies[x.second])
                     < bpi(y.first.size(), frequencies[y.second]);
            };

            // std::priority_queue<entry_type,
            //                     std::vector<entry_type>,
            //                     decltype(compare_bpi)> max_heap(compare_bpi);

            heap<entry, decltype(compare_bpi)> max_heap(top_k * 5, compare_bpi);

            for (uint32_t block_size = 16; block_size != 0; block_size /= 2)
            {
                logger() << "loading " << top_k << " most frequent blocks of size " << block_size << "..." << std::endl;
                std::string collection_name("./" + prefix_name + ".blocks_stats." + std::to_string(block_size) + ".bin");
                binary_blocks_collection input(collection_name.c_str());
                auto begin = input.begin();
                for (uint32_t i = 0; i < top_k; ++i, ++begin)
                {
                    uint32_t index = frequencies.size();
                    entry_type entry(std::vector<uint32_t>(), index);
                    auto const& block = *begin;
                    entry.first.reserve(block.size());
                    for (uint32_t x: block) {
                        entry.first.push_back(x);
                    }

                    max_heap.push(std::move(entry));
                    frequencies.push_back(block.freq());

                    uint8_t const* b = reinterpret_cast<uint8_t const*>(entry.first.data());
                    uint8_t const* e = b + block.size() * sizeof(uint32_t);
                    uint64_t hash = hash_bytes64(byte_range(b, e));
                    map[hash] = index;
                }
            }

            // 3. keep adding entries to the dictionary until it fills up:
            //   - at each step we remove the entry from the heap (in O(1))
            //     that gives the highest reduction in used bits per integer (bpi);
            //   - we decrease the frequency of all its sub-entries (if found) and
            //     call heapify() to restore the heap condition
            while (not builder.full() or not max_heap.empty()) {
                auto const& best = max_heap.top();
                auto const& best_block = best.first;
                uint64_t best_block_freq = frequencies[best.second];
                double cost_saving = bpi(best_block.size(), best_block_freq);
                bpi -= cost_saving;
                logger() << "added a target of size " << best.first.size() << " to the dictionary" << std::endl;
                logger() << "current bits x integer: " << bpi << std::endl;

                logger() << "heapifing..." << std::endl;
                for (uint32_t block_size = best_block.size() / 2; block_size != 0; block_size /= 2) {
                    for (uint32_t begin = 0; begin != end; begin += block_size) {
                        uint8_t const* b = reinterpret_cast<uint8_t const*>(&best_block[begin]);
                        uint8_t const* e = reinterpret_cast<uint8_t const*>(&best_block[end]);
                        uint64_t hash = hash_bytes64(byte_range(b, e));
                        uint32_t index = map[hash];
                        frequencies[index] -= best_block_freq;
                        assert(frequencies[index] >= 0);
                        heap.heapify();
                    }
                }
            }

            builder.build(dict);
            logger() << "using " << bpi << " bits x integer" << std::endl;
        }

    };
}
