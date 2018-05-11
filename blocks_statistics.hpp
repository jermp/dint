#pragma once

#include "hash_utils.hpp"

#include <fstream>
#include <unordered_map>

namespace ds2i {

    namespace util {
        bool is_power_of_two(uint64_t x) {
            return x & (x - 1) == 0;
        }
    }

    struct blocks_statistics {

        // <frequency of block, block>
        typedef std::pair<uint64_t, std::vector<uint32_t>> value_type;

        blocks_statistics(uint32_t block_size)
            : m_block_size(block_size)
            , m_bytes_per_block(block_size * sizeof(uint32_t))
        {
            assert(util::is_power_of_two(m_block_size));
        }

        void process(uint32_t const* begin, uint64_t n)
        {
            uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
            uint8_t const* e = b + n * sizeof(uint32_t);
            while (b != e) {
                size_t num_bytes = std::min<uint64_t>(bytes_per_block(), e - b);
                uint64_t hash = hash_bytes64(byte_range(b, b + num_bytes));

                auto it = m_map.find(hash);
                if (it == m_map.end()) {

                    value_type pair(1, std::vector<uint32_t>());
                    auto& block = pair.second;
                    uint64_t m = std::min<uint64_t>(block_size(), num_bytes / sizeof(uint32_t));
                    block.reserve(m);
                    auto ptr = reinterpret_cast<uint32_t const*>(b);
                    while (block.size() < m) {
                        block.push_back(*ptr);
                        ++ptr;
                    }
                    m_map[hash] = std::move(pair);

                    if (m_map.size() % 1000000 == 0) {
                        logger() << m_map.size() << " distinct blocks of size " << block_size() << std::endl;
                    }

                } else {
                    (*it).second.first += 1;
                }

                begin += num_bytes;
            }
        }

        void sort_and_write(std::string const& output_filename)
        {
            std::vector<value_type> freq_blocks;
            freq_blocks.reserve(occs.size());
            for (auto& pair: occs) {
                value_type p(pair.second.first, gaps_type());
                p.second.swap(pair.second.second);
                freq_blocks.push_back(std::move(p));
            }

            logger() << "sorting by decreasing freq..." << std::endl;
            std::sort(freq_blocks.begin(), freq_blocks.end(),
                      [](auto const& x, auto const& y) {
                            return x.first > y.first;
                      });

            std::ofstream out(output_filename.c_str());

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
                for (auto x: block) {
                    out.write(reinterpret_cast<char const*>(&x), bytes);
                }
            }

            out.close();
        }

        uint64_t bytes_per_block() const {
            return block_size() * sizeof(uint32_t);
        }

        uint32_t block_size() const {
            return m_block_size;
        }

    private:
        uint32_t m_block_size;
        uint32_t m_bytes_per_block;
        std::unordered_map<uint64_t, value_type> m_map;
    };

}
