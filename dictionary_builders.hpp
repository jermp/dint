#pragma once

#include "binary_blocks_collection.hpp"
#include "dictionary.hpp"
#include "heap.hpp"
#include "hash_utils.hpp"

#include <unordered_map>

namespace ds2i {

    namespace util {
        double cost(uint32_t block_size, uint32_t block_frequency) {
            return block_frequency * (48.0 * block_size - 16.0);
        }
    }

    template<uint32_t num_entries, uint32_t entry_width>
    struct dint_dictionary_builder
    {
        // Giulio: we can exploit the property that the distribution
        // is highly skewed to speed up the computation.
        // Instead of selecting top_k entries from EACH statistics file,
        // we can use an adaptive distribution, i.e., select more entries
        // if they appear a lot.
        static const uint32_t top_k = 50000;

        // <block, index into the frequency array>
        typedef std::pair<std::vector<uint32_t>, uint32_t> entry_type;

        static double bpi(uint32_t block_size, uint32_t block_frequency, uint64_t total_integers) {
            return util::cost(block_size, block_frequency) / total_integers;
        };

        struct bpi_comparator {
            bpi_comparator()
                : m_total_integers(0)
                , m_frequencies(nullptr)
            {}

            bpi_comparator(std::vector<uint64_t> const* frequencies, uint64_t total_integers)
                : m_total_integers(total_integers)
                , m_frequencies(frequencies)
            {}

            bool operator()(entry_type const& x, entry_type const& y) {
                return bpi(x.first.size(), m_frequencies->operator[](x.second), m_total_integers)
                     < bpi(y.first.size(), m_frequencies->operator[](y.second), m_total_integers);
            }

        private:
            uint64_t m_total_integers;
            std::vector<uint64_t> const* m_frequencies;
        };


        typedef std::pair<std::vector<uint32_t>, uint64_t> entry_type2;

        struct bpi_comparator2 {
            bpi_comparator2()
                : m_total_integers(0)
                , m_data(nullptr)
            {}

            bpi_comparator2(std::vector<entry_type2> const* data, uint64_t total_integers)
                : m_total_integers(total_integers)
                , m_data(data)
            {}

            int operator()(uint32_t id_x, uint32_t id_y) const {
                auto const& entry_x = m_data->operator[](id_x);
                auto const& entry_y = m_data->operator[](id_y);

                double bpi_x = bpi(entry_x.first.size(), entry_x.second, m_total_integers);
                double bpi_y = bpi(entry_y.first.size(), entry_y.second, m_total_integers);

                if (bpi_x == bpi_y) return 0;
                if (bpi_x  > bpi_y) return 1;
                return -1;
            }

            bool operator()(entry_type2 const& x, entry_type2 const& y) const
            {
                return bpi(x.first.size(), x.second, m_total_integers)
                     > bpi(y.first.size(), y.second, m_total_integers);
            }

        private:
            uint64_t m_total_integers;
            std::vector<entry_type2> const* m_data;
        };

        template<typename C>
        struct set_heap {
            set_heap()
                : m_size(0)
            {}

            set_heap(C const& comparator)
                : m_size(0)
                , m_comparator(comparator)
            {
                m_heap.push_back(-1); // dummy id
            }

            void push_back(uint32_t id) {
                m_heap.push_back(id);
                ++m_size;
            }

            void make() {
                m_size = m_heap.size() - 1;
                m_positions.reserve(m_size);
                for (size_t pos = 1; pos <= m_size; ++pos) {
                    m_positions.push_back(pos);
                }
                for (uint32_t pos = m_size / 2; pos > 0; --pos) {
                    max_heapify(pos);
                }
            }

            uint32_t pop()
            {
                uint32_t max = m_heap[1];
                std::swap(m_heap[1], m_heap[m_size]);
                std::swap(m_positions[m_heap[1]],
                          m_positions[m_heap[m_size]]);
                --m_size;
                max_heapify(1);
                return max;
            }

            void remove(uint32_t id)
            {
                uint32_t pos = m_positions[id];

                std::swap(m_heap[pos],
                          m_heap[m_size]);
                std::swap(m_positions[m_heap[pos]],
                          m_positions[m_heap[m_size]]);

                --m_size;
                swim(pos);
                sink(pos);
            }

            void increase(uint32_t id) {
                swim(m_positions[id]);
            }

            void decrease(uint32_t id) {
                sink(m_positions[id]);
            }

            bool empty() {
                std::cout << "m_size " << m_size << std::endl;
                return m_size == 0;
            }

        private:
            size_t m_size;
            std::vector<uint32_t> m_positions; // maps an id to its current position in the heap
            std::vector<uint32_t> m_heap;
            C m_comparator;

            bool greater(uint32_t pos_x, uint32_t pos_y) const {
                return m_comparator(m_heap[pos_x], m_heap[pos_y]) > 0;
            }

            bool less(uint32_t pos_x, uint32_t pos_y) const {
                return m_comparator(m_heap[pos_x], m_heap[pos_y]) < 0;
            }

            void max_heapify(uint32_t pos)
            {
                assert(pos);
                uint32_t l = 2 * pos;
                uint32_t r = l + 1;
                uint32_t max_pos = pos;

                if (l <= m_size and greater(l, pos)) {
                    max_pos = l;
                }

                if (r <= m_size and greater(r, max_pos)) {
                    max_pos = r;
                }

                if (max_pos != pos) {
                    std::swap(m_heap[pos],
                              m_heap[max_pos]);
                    std::swap(m_positions[m_heap[pos]],
                              m_positions[m_heap[max_pos]]);
                    max_heapify(max_pos);
                }
            }

            void sink(uint32_t pos) {
                while (2 * pos <= m_size) {
                    uint32_t i = 2 * pos;
                    if (i < m_size and less(i, i + 1)) {
                        ++i;
                    }
                    if (not less(pos, i)) break;

                    std::swap(m_heap[pos],
                              m_heap[i]);
                    std::swap(m_positions[m_heap[pos]],
                              m_positions[m_heap[i]]);
                    pos = i;
                }
            }

            void swim(uint32_t pos) {
                while (pos > 1 and less(pos / 2, pos)) {
                    std::swap(m_heap[pos],
                              m_heap[pos / 2]);
                    std::swap(m_positions[m_heap[pos]],
                              m_positions[m_heap[pos / 2]]);
                    pos = pos / 2;
                }
            }
        };

        static void build(dictionary::builder& builder, double const* percentages,
                          uint64_t total_integers, std::string prefix_name)
        {
            builder.init(num_entries, entry_width);

            // at the beginning, everything is exception(al)
            double final_bpi = 48.0;

            // 1. select the top-k entries from each statistics file
            // and put:
            //   - the blocks into a max heap;
            //   - the frequencies into a vector;
            // 2. create an hash table mapping blocks to indices of such vector
            // (in order to access efficiently the frequencies)

            std::vector<uint64_t> blocks_frequencies;
            blocks_frequencies.reserve(top_k * 5); // top-k entries for each different block_size

            bpi_comparator bc(&blocks_frequencies, total_integers);

            typedef heap<entry_type, bpi_comparator> max_heap_type;
            std::vector<max_heap_type> max_heaps(5, max_heap_type(top_k, bc));

            // <hash of block, index into frequency array>
            std::unordered_map<uint64_t, uint32_t> map;

            for (uint32_t block_size = 16, i = 0; block_size != 0; block_size /= 2, ++i)
            {
                // logger() << "loading " << top_k << " most frequent blocks of size " << block_size << "..." << std::endl;
                std::string collection_name("./" + prefix_name + ".blocks_stats." + std::to_string(block_size) + ".bin");
                binary_blocks_collection input(collection_name.c_str());
                auto begin = input.begin();
                for (uint32_t k = 0; k < top_k; ++k, ++begin)
                {
                    uint32_t index = blocks_frequencies.size();
                    entry_type entry(std::vector<uint32_t>(), index);
                    auto const& block = *begin;

                    // if (k < 10) std::cout << "block.freq() = " << block.freq() << std::endl;

                    entry.first.reserve(block.size());
                    for (uint32_t x: block) {
                        entry.first.push_back(x);
                    }

                    max_heaps[i].push_back(std::move(entry));
                    blocks_frequencies.push_back(block.freq());

                    uint8_t const* b = reinterpret_cast<uint8_t const*>(entry.first.data());
                    uint8_t const* e = b + block.size() * sizeof(uint32_t);
                    uint64_t hash = hash_bytes64(byte_range(b, e));
                    map[hash] = index;
                }
                max_heaps[i].make();
            }


            // for (int i = 0; i < 10; ++i) {
            //     auto const& best = max_heap.top();
            //     auto const& best_block = best.first;
            //     std::cout << "block.size() " << best_block.size() << std::endl;
            //     std::cout << "block.freq() " << blocks_frequencies[best.second] << std::endl;
            //     max_heap.pop();
            // }

            // 3. keep adding entries to the dictionary until it fills up:
            //   - at each step we remove the entry from the heap (in O(1))
            //     that gives the highest reduction in used bits per integer (bpi);
            //   - we decrease the frequency of all its sub-entries (if found) and
            //     call heapify() to restore the heap condition
            // uint64_t covered_integers = 0;

            double total_coverage = 0.0;
            uint32_t added_entries = 0;

            for (int i = 0; i < 5 and not builder.full(); ++i) {
                double p = 0.0;
                while (p < percentages[i])
                {
                    if (max_heaps[i].empty()) break;
                    auto const& best = max_heaps[i].top();
                    auto const& best_block = best.first;

                    if (not builder.append(best_block.data(), best_block.size())) {
                        break;
                    }

                    ++added_entries;
                    uint64_t best_block_freq = blocks_frequencies[best.second];
                    p += best_block_freq * best_block.size() * 100.0 / total_integers;
                    double cost_saving = bpi(best_block.size(), best_block_freq, total_integers);
                    final_bpi -= cost_saving;

                    if (added_entries % 500 == 0) {
                        std::cout << "added_entries " << added_entries << "/65536" << std::endl;
                        std::cout << "p = " << p << "/" << percentages[i] << std::endl;
                        // std::cout << "saving " << cost_saving << " bits x int" << std::endl;

                        // std::cout << "adding a target of size " << best.first.size() << " to the dictionary" << std::endl;

                        // for (auto x: best_block) {
                        //     std::cout << x << " ";
                        // }
                        // std::cout << std::endl;

                        std::cout << "current bits x integer: " << final_bpi << std::endl;
                        std::cout << "covering " << total_coverage + p << "% of integers" << std::endl;
                    }

                    // logger() << "decreasing freq. of sub-blocks..." << std::endl;
                    for (uint32_t block_size = best_block.size() / 2, j = i + 1; block_size != 0; block_size /= 2, ++j) {
                        // std::cout << "decreasing freqs of sub-blocks of size " << block_size << std::endl;
                        for (uint32_t begin = 0; begin < best_block.size(); begin += block_size) {

                            // uint32_t end = std::min<uint32_t>(begin + block_size, best_block.size());
                            // std::cout << "block[" << begin << ", " << end << ") = ";
                            // for (uint32_t kk = begin; kk != end; ++kk) {
                            //     std::cout << best_block[kk] << " ";
                            // }
                            // std::cout << std::endl;

                            uint8_t const* b = reinterpret_cast<uint8_t const*>(&best_block[begin]);
                            uint8_t const* e = b + std::min<uint64_t>(block_size, best_block.size() - begin) * sizeof(uint32_t);
                            uint64_t hash = hash_bytes64(byte_range(b, e));
                            uint32_t index = map[hash];

                            // std::cout << "decreasing " << blocks_frequencies[index] << " by " << best_block_freq << std::endl;
                            blocks_frequencies[index] -= best_block_freq;
                            assert(blocks_frequencies[index] >= 0);
                        }
                        // Giulio: if we keep the positions of the elements in the heap,
                        // we can use sink(position): O(log n) vs O(n)
                        max_heaps[j].make();
                    }

                    max_heaps[i].pop();
                }

                total_coverage += p;
            }

            logger() << "using " << final_bpi << " bits x integer" << std::endl;
            logger() << "covering " << total_coverage << "% of integers" << std::endl;
        }

        static void build1(dictionary::builder& builder, uint64_t total_integers,
                          std::string prefix_name, double eps = 0.0001)
        {
            builder.init(num_entries, entry_width);

            // at the beginning, everything is exception(al)
            double final_bpi = 48.0;


            std::vector<entry_type2> candidates;
            bpi_comparator2 bc2(&candidates, total_integers);
            set_heap<bpi_comparator2> ids(bc2);

            // <hash of block, unique_id>
            std::unordered_map<uint64_t, uint32_t> map;

            // push to the set_heap the candidates
            for (uint32_t block_size = 16; block_size != 0; block_size /= 2)
            {
                std::string collection_name("./" + prefix_name + ".blocks_stats." + std::to_string(block_size) + ".bin");
                binary_blocks_collection input(collection_name.c_str());
                uint32_t entries = 0;

                for (auto begin = input.begin(); begin != input.end(); ++begin)
                {
                    auto const& block = *begin;
                    double saving = bpi(block.size(), block.freq(), total_integers);
                    // std::cout << saving << std::endl;
                    if (saving > eps) {
                        uint32_t id = candidates.size();
                        entry_type2 entry(std::vector<uint32_t>(), block.freq());
                        entry.first.reserve(block.size());
                        for (uint32_t x: block) {
                            entry.first.push_back(x);
                        }
                        ids.push_back(id);

                        uint8_t const* b = reinterpret_cast<uint8_t const*>(entry.first.data());
                        uint8_t const* e = b + block.size() * sizeof(uint32_t);
                        uint64_t hash = hash_bytes64(byte_range(b, e));
                        map[hash] = id;
                        ++entries;

                        candidates.push_back(std::move(entry));

                    } else {
                        break; // bpi cost is decreasing, since blocks are sorted in decreasing frequency
                    }
                }

                logger() << "added " << entries <<  " " << block_size << "-int entries out of "
                         << input.num_blocks() << " (" << entries * 100.0 / input.num_blocks() << "%)"
                         << std::endl;
            }

            ids.make();

            double total_coverage = 0.0;
            uint32_t added_entries = 0;
            std::vector<uint32_t> larger_block;
            larger_block.reserve(16);

            while (not builder.full())
            {
                if (ids.empty()) break;

                uint32_t id = ids.pop();
                auto const& block_data = candidates[id];
                auto const& block = block_data.first;
                uint64_t freq = block_data.second;
                std::cout << "frequency = " << freq << std::endl;

                logger() << "selecting entry of size " << block.size() << ": ";
                for (auto x: block) {
                    std::cout << x << " ";
                }
                std::cout << std::endl;

                builder.append(block.data(), block.size());
                ++added_entries;


                total_coverage += freq * block.size() * 100.0 / total_integers;
                double cost_saving = bpi(block.size(), freq, total_integers);
                final_bpi -= cost_saving;

                // if (added_entries % 500 == 0) {
                    std::cout << "added_entries " << added_entries << "/65536" << std::endl;
                    std::cout << "bits x integer: " << final_bpi << std::endl;
                    std::cout << "covering " << total_coverage << "% of integers" << std::endl;
                // }

                // logger() << "decreasing freqs of small blocks..." << std::endl;
                // decrease frequencies of smaller blocks (if any)
                for (uint32_t block_size = block.size() / 2; block_size != 0; block_size /= 2) {
                    for (uint32_t begin = 0; begin < block.size(); begin += block_size) {
                        uint8_t const* b = reinterpret_cast<uint8_t const*>(&block[begin]);
                        uint8_t const* e = b + std::min<uint64_t>(block_size, block.size() - begin) * sizeof(uint32_t);
                        uint64_t hash = hash_bytes64(byte_range(b, e));
                        uint32_t id = map[hash];
                        std::cout << "decreasing " << candidates[id].second << " by " << freq << std::endl;
                        candidates[id].second -= freq;
                        assert(candidates[id].second >= 0);
                        ids.decrease(id);
                    }
                }

                logger() << "decreasing freq. of big blocks..." << std::endl;

                // delete all larger blocks (if any)
                for (uint32_t block_size = block.size() * 2; block_size != 32; block_size *= 2)
                {
                    larger_block.clear();
                    for (uint32_t i = 0; i < block_size / block.size(); ++i) {
                        for (auto x: block) {
                            larger_block.push_back(x);
                            std::cout << x << " ";
                        }
                        std::cout << " - ";
                    }
                    std::cout << std::endl;

                    uint8_t const* b = reinterpret_cast<uint8_t const*>(larger_block.data());
                    uint8_t const* e = b + larger_block.size() * sizeof(uint32_t);
                    uint64_t hash = hash_bytes64(byte_range(b, e));
                    auto it = map.find(hash);
                    if (it != map.end()) {
                        std::cout << "found" << std::endl;
                        uint32_t id = (*it).second;
                        ids.remove(id);
                    } else {
                        std::cout << "not found" << std::endl;
                    }
                }
            }

            logger() << "using " << final_bpi << " bits x integer" << std::endl;
            logger() << "covering " << total_coverage << "% of integers" << std::endl;
        }

        static const uint32_t MAX_BLOCK_LEN = 8;     // 16
        static const uint32_t MAX_FRACTAL_STEPS = 4; // 5

        static void build2(dictionary::builder& builder, uint64_t total_integers,
                           std::string prefix_name, double eps = 0.0001)
        {
            builder.init(num_entries, entry_width);

            // at the beginning, everything is exception(al)
            double final_bpi = 48.0;

            std::vector<std::vector<entry_type2>> candidates(5, std::vector<entry_type2>());
            // typedef set_heap<bpi_comparator2> set_heap_type;
            // std::vector<set_heap_type> ids(5, set_heap_type(bc2));

            // <hash of block, unique_id>
            std::unordered_map<uint64_t, uint32_t> map;

            std::vector<uint64_t> id_lowerbounds(MAX_FRACTAL_STEPS + 1, 0);
            // push to the set_heap the candidates
            uint32_t i = 0;
            uint32_t id = 0;
            for (uint32_t block_size = MAX_BLOCK_LEN; block_size != 0; block_size /= 2)
            {
                std::string collection_name("./" + prefix_name + ".blocks_stats." + std::to_string(block_size) + ".bin");
                binary_blocks_collection input(collection_name.c_str());
                uint32_t entries = 0;
                auto& c = candidates[i];

                for (auto begin = input.begin(); begin != input.end(); ++begin)
                {
                    auto const& block = *begin;
                    double saving = bpi(block.size(), block.freq(), total_integers);
                    // std::cout << saving << std::endl;
                    if (saving > eps or (i == MAX_FRACTAL_STEPS - 1 and (id < num_entries))) {
                        entry_type2 entry(std::vector<uint32_t>(), block.freq());
                        entry.first.reserve(block.size());
                        for (uint32_t x: block) {
                            entry.first.push_back(x);
                        }
                        uint8_t const* b = reinterpret_cast<uint8_t const*>(entry.first.data());
                        uint8_t const* e = b + block.size() * sizeof(uint32_t);
                        uint64_t hash = hash_bytes64(byte_range(b, e));
                        map[hash] = id;
                        ++id;
                        ++entries;
                        c.push_back(std::move(entry));
                    } else {
                        break; // bpi cost is decreasing, since blocks are sorted in decreasing frequency
                    }
                }

                ++i;
                id_lowerbounds[i] = id;

                logger() << "added " << entries <<  " " << block_size << "-int entries out of "
                         << input.num_blocks() << " (" << entries * 100.0 / input.num_blocks() << "%)"
                         << std::endl;
            }

            double total_coverage = 0.0;
            uint32_t added_entries = 0;

            i = 0;
            while (not builder.full() and i != MAX_FRACTAL_STEPS)
            {
                for (auto const& x: candidates[i])
                {
                    auto const& block = x.first;
                    uint64_t freq = x.second;
                    double cost_saving = bpi(block.size(), freq, total_integers);
                    if (i != 4 and cost_saving < eps) break;

                    // std::cout << "frequency = " << freq << std::endl;
                    // std::cout << "selected entry of size " << block.size() << ": ";
                    // for (auto x: block) {
                    //     std::cout << x << " ";
                    // }
                    // std::cout << std::endl;

                    builder.append(block.data(), block.size());
                    ++added_entries;
                    total_coverage += freq * block.size() * 100.0 / total_integers;
                    final_bpi -= cost_saving;

                    if (added_entries % 500 == 0) {
                        logger() << "entries in dictionary " << added_entries << "/65536" << std::endl;
                        logger() << "current bits x integer: " << final_bpi << std::endl;
                        logger() << "covering " << total_coverage << "% of integers" << std::endl;
                    }

                    // decrease frequencies of smaller blocks (if any)
                    for (uint32_t block_size = block.size() / 2, j = i + 1; block_size != 0; block_size /= 2, ++j) {
                        for (uint32_t begin = 0; begin < block.size(); begin += block_size) {
                            uint8_t const* b = reinterpret_cast<uint8_t const*>(&block[begin]);
                            uint8_t const* e = b + std::min<uint64_t>(block_size, block.size() - begin) * sizeof(uint32_t);
                            uint64_t hash = hash_bytes64(byte_range(b, e));
                            uint32_t id = map[hash];
                            // std::cout << "decreasing " << candidates[j][id - id_lowerbounds[j]].second << " by " << freq << std::endl;
                            candidates[j][id - id_lowerbounds[j]].second -= freq;
                            assert(candidates[j][id - id_lowerbounds[j]].second >= 0);
                        }
                    }
                }

                logger() << "covering " << total_coverage << "% of integers "
                         << "with entries of size " << (MAX_BLOCK_LEN >> i) << std::endl;

                if (i != MAX_FRACTAL_STEPS - 1) {
                    auto& c = candidates[i + 1];
                    bpi_comparator2 comp(&c, total_integers);
                    std::sort(c.begin(), c.end(), comp);
                }
                ++i;
            }

            logger() << "using " << final_bpi << " bits x integer" << std::endl;
            logger() << "covering " << total_coverage << "% of integers" << std::endl;
        }
    };
}
