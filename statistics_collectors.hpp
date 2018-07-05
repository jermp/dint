#pragma once

#include "hash_utils.hpp"
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

    struct freq_sorter {
        bool operator()(block_type const& l, block_type const& r) {
            return l.freq > r.freq;
        }
    };

    struct length_freq_sorter {
        bool operator()(block_type const& l, block_type const& r) {
            if (l.data.size() == r.data.size()) {
                return l.freq > r.freq;
            }
            return l.data.size() > r.data.size();
        }
    };

    struct freq_length_sorter {
        bool operator()(block_type const& l, block_type const& r) {
            if (l.freq == r.freq) {
                return l.data.size() > r.data.size();
            }
            return l.freq > r.freq;
        }
    };

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

            for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
                uint32_t block_size = constants::target_sizes[s];
                uint32_t blocks = buf.size() / block_size;
                for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += block_size) {
                    increase_frequency(b + pos, block_size, bmap);
                }
            }

            // for (uint32_t block_size = max_block_size; block_size != 0; block_size /= 2) {
            //     uint32_t blocks = buf.size() / block_size;
            //     for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += block_size) {
            //         increase_frequency(b + pos, block_size, bmap);
            //     }
            // }
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
                for (uint32_t block_size = max_block_size; block_size != 0; block_size /= 2) {
                    uint32_t amount = max_block_size / block_size;
                    increase_frequency(b + pos, block_size, bmap, amount);
                }
            }
        }
    };

    template<uint32_t t_max_block_size>
    struct fixed
    {
        static const uint32_t max_block_size = t_max_block_size;

        static std::string type() {
            return "fixed";
        }

        static void collect(std::vector<uint32_t>& buf, map_type& bmap) {
            auto b = buf.data();
            uint32_t blocks = buf.size() / max_block_size;
            for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += max_block_size) {
                increase_frequency(b + pos, max_block_size, bmap);
            }
        }
    };

};
