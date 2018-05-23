#pragma once

#include "hash_utils.hpp"
#include "model_build_utils.hpp"

#include <fstream>
#include <unordered_map>

namespace ds2i {

    struct blocks_statistics {

        // <frequency of block, block>
        typedef std::pair<uint64_t, std::vector<uint32_t>> value_type;

        blocks_statistics(uint32_t block_size)
            : m_block_size(block_size)
        {
            logger() << "initialize blocks_statistics with block_size = " << block_size << std::endl;
            assert(is_power_of_two(m_block_size));
        }

        void process(uint32_t const* begin, uint64_t n)
        {
            // uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
            // uint8_t const* e = b + n * sizeof(uint32_t);
            // while (b != e) {
            //     size_t num_bytes = std::min<uint64_t>(bytes_per_block(), e - b);
            //     uint64_t hash = hash_bytes64(byte_range(b, b + num_bytes));

            //     auto it = m_map.find(hash);
            //     if (it == m_map.end()) {

            //         value_type pair(1, std::vector<uint32_t>());
            //         auto& block = pair.second;
            //         uint64_t m = std::min<uint64_t>(block_size(), num_bytes / sizeof(uint32_t));
            //         block.reserve(m);
            //         auto ptr = reinterpret_cast<uint32_t const*>(b);
            //         while (block.size() < m) {
            //             block.push_back(*ptr);
            //             ++ptr;
            //         }
            //         m_map[hash] = std::move(pair);

            //         if (m_map.size() % 1000000 == 0) {
            //             logger() << m_map.size() << " distinct blocks of size " << block_size() << std::endl;
            //         }

            //     } else {
            //         (*it).second.first += 1;
            //     }

            //     b += num_bytes;
            // }

            uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
            uint8_t const* e = b + n * sizeof(uint32_t);
            uint64_t bytes = bytes_per_block();
            while (b < e) { // ignore block tail
                uint64_t hash = hash_bytes64(byte_range(b, b + bytes));
                auto it = m_map.find(hash);
                if (it == m_map.end())
                {
                    value_type pair(1, std::vector<uint32_t>());
                    auto& block = pair.second;
                    block.reserve(block_size());
                    auto ptr = reinterpret_cast<uint32_t const*>(b);
                    while (block.size() < block_size()) {
                        block.push_back(*ptr++);
                    }
                    m_map[hash] = std::move(pair);

                    if (m_map.size() % 1000000 == 0) {
                        logger() << m_map.size() << " distinct blocks of size " << block_size() << std::endl;
                    }

                } else {
                    (*it).second.first += 1;
                }

                b += bytes;
            }
        }

        void sort_and_write(std::string const& output_filename,
                            uint64_t total_integers, double eps = 0.0001)
        {
            std::vector<value_type> freq_blocks;
            freq_blocks.reserve(m_map.size());
            for (auto& pair: m_map) {
                if (bpi(block_size(), pair.second.first, total_integers) > eps) { // cost pruning
                    value_type p(pair.second.first, std::vector<uint32_t>());
                    p.second.swap(pair.second.second);
                    freq_blocks.push_back(std::move(p));
                }
            }

            logger() << "sorting " << freq_blocks.size() << " blocks by decreasing freq..." << std::endl;
            std::sort(freq_blocks.begin(), freq_blocks.end(),
                      [](value_type const& x, value_type const& y) {
                            return x.first > y.first;
                      });

            std::ofstream out(output_filename.c_str());

            logger() << "writing binary file..." << std::endl;
            // write header
            uint32_t num_blocks = freq_blocks.size();
            std::streamsize bytes = sizeof(uint32_t);
            out.write(reinterpret_cast<char const*>(&num_blocks), bytes);
            for (uint64_t i = 0; i < freq_blocks.size(); ++i) {
                auto const& block = freq_blocks[i].second;
                uint32_t size = block.size();
                uint32_t freq = freq_blocks[i].first;
                out.write(reinterpret_cast<char const*>(&size), bytes);
                out.write(reinterpret_cast<char const*>(&freq), bytes);
                out.write(reinterpret_cast<char const*>(block.data()), size * bytes);
            }
            out.close();

            // free memory
            logger() << "releasing memory..." << std::endl;
            std::unordered_map<uint64_t, value_type>().swap(m_map);
            logger() << "DONE" << std::endl;
        }

        uint64_t bytes_per_block() const {
            return block_size() * sizeof(uint32_t);
        }

        uint32_t block_size() const {
            return m_block_size;
        }

    private:
        uint32_t m_block_size;
        std::unordered_map<uint64_t, value_type> m_map;
    };

}
