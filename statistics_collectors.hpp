#pragma once

namespace ds2i {

    static const size_t MAX_STRIDE_SIZE = 64;

    #pragma pack(push, 1)
    template<uint32_t width>
    struct block_info {
        uint64_t hash;
        uint64_t freq;
        uint8_t entry_len;
        uint32_t entry[width];
    };
    #pragma pack(pop)

    // Giulio: why not an abstract class with the increase_frequency method
    // and the block_map ??

    template<typename bm_type>
    void increase_frequency(const uint32_t* entry, size_t n, bm_type& block_map,
                            uint32_t amount = 1)
    {
        using b_type = typename bm_type::value_type::second_type;
        auto hash = hash_u32(entry,n);
        auto itr = block_map.find(hash);
        if (itr != block_map.end()) {
            itr->second.freq += amount;
        } else {
            b_type new_block;
            new_block.hash = hash;
            new_block.freq = 1;
            new_block.entry_len = n;
            std::copy(entry,entry+n,new_block.entry);
            block_map[hash] = new_block;
        }
    }

    // struct geometric // Giulio: what is the reason for this??
    // {
    //     static std::string type() {
    //         return "geometric";
    //     }

    //     template<class bm_type>
    //     static void collect(const std::vector<uint32_t>& buf, bm_type& block_map)
    //     {
    //         auto b = buf.data();
    //         for (size_t stride_size = 1; stride_size <= buf.size(); stride_size *= 2) {
    //             increase_frequency(b, size_u32, block_map);
    //         }
    //     }
    // };

    struct adjusted
    {
        static std::string type() {
            return "adjusted";
        }

        template<class bm_type>
        static void collect(std::vector<uint32_t>& buf, bm_type& block_map)
        {
            auto b = buf.data();
            for (size_t stride_size = 1; stride_size <= MAX_STRIDE_SIZE; stride_size *= 2) {
                for (size_t pos = 0; pos < buf.size(); pos += stride_size) {
                    increase_frequency(b + pos, stride_size, block_map);
                }
            }
        }
    };

    struct full
    {
        static std::string type() {
            return "full";
        }

        template<class bm_type>
        static void collect(std::vector<uint32_t>& buf, bm_type& block_map) {
            auto b = buf.data();
            for (size_t pos = 0; pos < buf.size(); pos += MAX_STRIDE_SIZE) {
                for (size_t stride_size = 1; stride_size <= MAX_STRIDE_SIZE; stride_size *= 2)
                {
                    uint32_t amount = MAX_STRIDE_SIZE / stride_size;
                    increase_frequency(b + pos, stride_size, block_map, amount);
                }
            }
        }
    };

};
