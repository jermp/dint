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

    struct selector {
        // selector(uint32_t block_size)
        //     : m_block_size(block_size)
        // {
        //     // assert(m_block_size % 2 == 0);
        //     // if (constants::context == constants::block_selector::median) {
        //     //     m_buf.reserve(m_block_size);
        //     // }
        //     // if (constants::context == constants::block_selector::mode) {
        //     //     m_occs.reserve(m_block_size);
        //     // }
        // }

        uint32_t get(uint32_t const* entry, size_t n) {
            uint32_t x = 0;
            if (constants::context == constants::block_selector::max) {
                for (auto ptr = entry; ptr != entry + n; ++ptr) {
                    if (*ptr > x) {
                        x = *ptr;
                    }
                }
                // std::cout << "block max " << x << "; ";
            }
            // else
            // if (constants::context == constants::block_selector::median) {
            //     auto ptr = entry;
            //     for (uint32_t s = 0; s != m_block_size; ++s) {
            //         m_buf[s] = *ptr;
            //         ++ptr;
            //     }
            //     std::sort(m_buf.begin(), m_buf.end());
            //     x = (m_buf[m_block_size / 2 - 1] + m_buf[m_block_size / 2]) / 2;
            // } else
            // if (constants::context == constants::block_selector::mode)
            // {
            //     m_occs.clear();
            //     auto ptr = entry;
            //     uint32_t best_occ = 0;
            //     for (uint32_t s = 0; s != m_block_size; ++s)
            //     {
            //         auto it = m_occs.find(*ptr);
            //         if (it == m_occs.end()) {
            //             m_occs[*ptr] = 1;
            //         } else {
            //             uint32_t occ = (*it).second;
            //             if (occ > best_occ) {
            //                 best_occ = occ;
            //                 x = *ptr;
            //             }
            //         }
            //         ++ptr;
            //     }
            // }
            else {
                throw std::runtime_error("Unsupported context");
            }

            // uint32_t selector_code = (x == 0 ? 0 : ceil_log2(x)) + 1;
            // uint32_t selector_code = 0;
            // if (x > 1) {
            //     selector_code = ceil_log2(ceil_log2(x));
            // }

            uint32_t selector_code = 0;
            if (x > 1) {
                selector_code = ceil_log2(ceil_log2(x + 1));
            }

            return selector_code;
            // std::cout << "selector_code " << selector_code << "; ";
            // uint32_t index = selector_code;
            // uint32_t index = 0;

            // while (selector_code > constants::selector_codes[index]
            //        and index != constants::num_selectors) {
            //     ++index;
            // }

            // // std::cout << "index " << index << "\n";
            // assert(index < constants::num_selectors);

            // return index;
        }

    // private:
    //     uint32_t m_block_size;
        // std::vector<uint32_t> m_buf;
        // std::unordered_map<uint32_t, uint32_t> m_occs;
    };

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

    template<uint32_t t_max_block_size>
    struct adjusted
    {
        static const uint32_t max_block_size = t_max_block_size;

        static std::string type() {
            return "adjusted";
        }

        static void collect(std::vector<uint32_t>& buf, std::vector<map_type>& block_maps) {
            auto b = buf.data();
            uint32_t blocks = buf.size() / constants::block_size;
            selector sct;
            for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += constants::block_size)
            {
                uint32_t index = sct.get(b + pos, constants::block_size);
                for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
                    uint32_t jump_size = constants::target_sizes[s];
                    uint32_t jumps = constants::block_size / jump_size;
                    for (uint32_t j = 0, p = 0; j < jumps; ++j, p += jump_size) {
                        increase_frequency(b + pos + p, jump_size, block_maps[index]);
                    }
                }
            }
        }

        static void collect(std::vector<uint32_t>& buf, map_type& block_map) {
            auto b = buf.data();
            for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
                uint32_t block_size = constants::target_sizes[s];
                uint32_t blocks = buf.size() / block_size;
                for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += block_size) {
                    increase_frequency(b + pos, block_size, block_map);
                }
            }
        }
    };

    // template<uint32_t t_max_block_size>
    // struct full
    // {
    //     static const uint32_t max_block_size = t_max_block_size;

    //     static std::string type() {
    //         return "full";
    //     }

    //     static void collect(std::vector<uint32_t>& buf, map_type& bmap) {
    //         auto b = buf.data();
    //         uint32_t blocks = buf.size() / max_block_size;
    //         for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += max_block_size) {
    //             for (uint32_t block_size = max_block_size; block_size != 0; block_size /= 2) {
    //                 uint32_t amount = max_block_size / block_size;
    //                 increase_frequency(b + pos, block_size, bmap, amount);
    //             }
    //         }
    //     }
    // };

    // template<uint32_t t_max_block_size>
    // struct fixed
    // {
    //     static const uint32_t max_block_size = t_max_block_size;

    //     static std::string type() {
    //         return "fixed";
    //     }

    //     static void collect(std::vector<uint32_t>& buf, map_type& bmap) {
    //         auto b = buf.data();
    //         uint32_t blocks = buf.size() / max_block_size;
    //         for (uint32_t i = 0, pos = 0; i < blocks; ++i, pos += max_block_size) {
    //             increase_frequency(b + pos, max_block_size, bmap);
    //         }
    //     }
    // };

};
