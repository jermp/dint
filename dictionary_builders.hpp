#pragma once

#include "dictionary.hpp"

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

        static void build(dictionary& dict, uint64_t total_num_integers)
        {
            auto bpi = [&](uint32_t block_size, uint32_t block_frequency) {
                return util::cost(block_size, block_frequency) / total_num_integers;
            };

            dictionary::builder builder(num_entries, entry_width);

            // assume that everything NOT covered by the dictionary
            // is left uncompressed
            double final_bpi = 32.0;

            // 1. select the top-1M entries from each statistics file
            // and put:
            //   - the blocks into a max heap;
            //   - the frequencies into a vector;
            // 2. create an hash table mapping blocks to indices of such vector
            // (in order to access efficiently the frequencies)
            for (uint32_t block_size = 16; block_size != 0; block_size /= 2)
            {

            }

            // 3. keep adding entries to the dictionary until it fills up:
            //   - at each step we remove the entry from the heap (in O(1))
            //     that gives the highest reduction in used bits per integer (bpi);
            //   - we decrease the frequency of all its sub-entries (if found);
            //   - we re-establish the max-heap condition (O(log n) for each
            //     decrease-freq operation)
            while (not builder.full() or not heap.empty()) {

            }

            builder.build(dict);
            logger() << "using " << bpi << " bits x integer" << std::endl;
        }

    };
}
