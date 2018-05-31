#pragma once

#include "dint_configuration.hpp"

namespace ds2i {

    #pragma pack(push, 1)
    template<uint32_t N>
    struct block {
        uint64_t hash;
        uint64_t freq;
        uint8_t size;
        uint32_t entry[N];
    };
    #pragma pack(pop)

    typedef std::pair<uint64_t, block> hash_block_type;
    typedef std::unordered_map<uint64_t, block> blocks_map;


    // Giulio: why not an abstract class with the increase_frequency method
    // and a blocks_map object ??



    void increase_frequency(const uint32_t* entry, size_t n,
                            blocks_map& bmap, uint32_t amount = 1)
    {
        auto hash = hash_bytes64(entry, n);
        auto it = bmap.find(hash);
        if (it != bmap.end()) {
            it->second.freq += amount;
        } else {
            block b;
            b.hash = hash;
            b.freq = 1;
            b.size = n;
            std::copy(entry, entry + n, b.entry);
            bmap[hash] = std::move(b);
        }
    }

    // struct geometric // Giulio: what is the reason for this??
    // {
    //     static std::string type() {
    //         return "geometric";
    //     }

    //     template<class Map>
    //     static void collect(const std::vector<uint32_t>& buf, Map& block_map)
    //     {
    //         auto b = buf.data();
    //         for (size_t block_size = 1; block_size <= buf.size(); block_size *= 2) {
    //             increase_frequency(b, size_u32, block_map);
    //         }
    //     }
    // };

    struct adjusted
    {
        static std::string type() {
            return "adjusted";
        }

        static void collect(std::vector<uint32_t>& buf, blocks_map& bmap) {
            auto b = buf.data();
            for (uint32_t block_size = 1; block_size <= constants::max_block_length; block_size *= 2) {
                for (size_t pos = 0; pos < buf.size(); pos += block_size) {
                    increase_frequency(b + pos, block_size, bmap);
                }
            }
        }
    };

    struct full
    {
        static std::string type() {
            return "full";
        }

        static void collect(std::vector<uint32_t>& buf, blocks_map& bmap) {
            auto b = buf.data();
            for (size_t pos = 0; pos < buf.size(); pos += constants::max_block_length) {
                for (uint32_t block_size = 1; block_size <= constants::max_block_length; block_size *= 2) {
                    uint32_t amount = constants::max_block_length / block_size;
                    increase_frequency(b + pos, block_size, bmap, amount);
                }
            }
        }
    };

};
