#pragma once

#include "dint_configuration.hpp"

namespace ds2i {

    struct block_type {
        block_type()
            : freq(1)
        {}

        uint64_t hash() const {
            return hash_bytes64(data.data(),
                                data.size());
        }

        uint64_t freq;
        std::vector<uint32_t> data;
    };

    typedef std::unordered_map<uint64_t, block_type> map_type;

    void increase_frequency(uint32_t const* entry, size_t n,
                            map_type& bmap, uint32_t amount = 1)
    {
        auto hash = hash_bytes64(entry, n);
        auto it = bmap.find(hash);
        if (it != bmap.end()) {
            (*it).second.freq += amount;
        } else {
            block_type block;
            block.data.reserve(n);
            while (block.data.size() < n) {
                block.data.push_back(*entry++);
            }
            bmap[hash] = std::move(block);
        }
    }

    // struct geometric // Giulio: what is the reason for this?
    // {
    //     static std::string type() {
    //         return "geometric";
    //     }

    //     static void collect(const std::vector<uint32_t>& buf, map_type& block_map)
    //     {
    //         auto b = buf.data();
    //         for (size_t block_size = 1; block_size <= buf.size(); block_size *= 2) {
    //             increase_frequency(b, size_u32, block_map);
    //         }
    //     }
    // };

    template<uint32_t t_max_block_size>
    struct adjusted
    {
        static const uint32_t max_block_size = t_max_block_size;

        static std::string type() {
            return "adjusted";
        }

        static void collect(std::vector<uint32_t>& buf, map_type& bmap) {
            auto b = buf.data();
            for (uint32_t block_size = 1; block_size <= max_block_size; block_size *= 2) {
                uint32_t blocks = buf.size() / block_size;
                for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += block_size) {
                    increase_frequency(b + pos, block_size, bmap);
                }
            }
        }
    };

    template<uint32_t t_max_block_size>
    struct full
    {
        static const uint32_t max_block_size = t_max_block_size;

        static std::string type() {
            return "full";
        }

        static void collect(std::vector<uint32_t>& buf, map_type& bmap) {
            auto b = buf.data();
            uint32_t blocks = buf.size() / max_block_size;
            for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += max_block_size) {
                for (uint32_t block_size = 1; block_size <= max_block_size; block_size *= 2) {
                    uint32_t amount = max_block_size / block_size;
                    increase_frequency(b + pos, block_size, bmap, amount);
                }
            }
        }
    };

};
