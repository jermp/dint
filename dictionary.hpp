#pragma once

#include <unordered_map>

#include <succinct/mappable_vector.hpp>
#include <fstream>

#include "hash_utils.hpp"
#include "util.hpp"

namespace ds2i {

    struct dictionary {

        static const uint32_t invalid_index = uint32_t(-1);
        static const uint32_t reserved = 5; // 1 for exceptions
                                            // 4 for runs

        struct builder {

            builder()
                : m_pos(0)
                , m_size(reserved)
                , m_reserved(reserved)
                , m_capacity(0)
                , m_entry_size(0)
                , m_table(0, 0)
            {}

            void init(uint32_t capacity, uint32_t entry_size) {
                m_pos = reserved * (entry_size + 1);
                m_size = reserved;
                m_capacity = capacity;
                m_entry_size = entry_size;
                m_table.resize(capacity * (entry_size + 1), 0);
                m_freq.resize(capacity);
                m_metric.resize(capacity);

                m_freq.push_back(0);    // exceptions
                m_metric.push_back(0);  // exceptions

                m_freq.push_back(0);    // 256
                m_metric.push_back(0);  // 256

                m_freq.push_back(0);    // 128
                m_metric.push_back(0);  // 128

                m_freq.push_back(0);    // 64
                m_metric.push_back(0);  // 64

                m_freq.push_back(0);    // 32
                m_metric.push_back(0);  // 32

            }

            builder(uint32_t capacity, uint32_t entry_size) {
                init(capacity, entry_size);
            }

            bool full() {
                return m_size == m_capacity;
            }

            bool append(uint32_t const* entry, uint32_t entry_size,uint64_t freq,uint64_t metric) {
                if (full()) {
                    return false;
                }
                assert(m_pos % (m_entry_size + 1) == 0);
                std::copy(entry, entry + entry_size, &m_table[m_pos]);
                m_pos += m_entry_size + 1;
                m_table[m_pos - 1] = entry_size;
                m_freq.push_back(freq);
                m_metric.push_back(metric);
                ++m_size;
                return true;
            }

            void prepare_for_encoding() {
                logger() << "building mapping for encoding " << std::endl;
                std::vector<uint32_t> run(256, 1);
                uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
                uint32_t i = 0;
                for (uint32_t n = 256; n != 16; n /= 2, ++i) {
                    uint64_t hash = hash_bytes64(byte_range(ptr, ptr + n * sizeof(uint32_t)));
                    m_map[hash] = i;
                }
                for (; i < capacity(); ++i) {
                    uint8_t const* ptr = reinterpret_cast<uint8_t const*>(get(i));
                    uint32_t entry_size = size(i);
                    uint64_t hash = hash_bytes64(byte_range(ptr, ptr + entry_size * sizeof(uint32_t)));
                    m_map[hash] = i;
                }
            }

            // Giulio: return the index of the pattern if found in the table,
            // otherwise return the invalid_index value
            uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
            {
                uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
                uint8_t const* e = b + entry_size * sizeof(uint32_t);
                uint64_t hash = hash_bytes64(byte_range(b, e));
                auto it = m_map.find(hash);
                if (it != m_map.end()) {
                    assert((*it).second <= capacity());
                    return (*it).second;
                }
                return invalid_index;
            }

            size_t capacity() const {
                return m_capacity;
            }

            size_t entry_size() const {
                return m_entry_size;
            }

            void build(dictionary& dict) {
                std::swap(m_capacity, dict.m_capacity);
                std::swap(m_entry_size, dict.m_entry_size);
                dict.m_table.steal(m_table);
                builder().swap(*this);
            }

            void swap(builder& other) {
                std::swap(m_pos, other.m_pos);
                std::swap(m_size, other.m_size);
                std::swap(m_capacity, other.m_capacity);
                std::swap(m_entry_size, other.m_entry_size);
                m_table.swap(other.m_table);
                m_map.swap(other.m_map);
            }

            uint32_t size(uint32_t i) const {
                uint32_t begin = i * (m_entry_size + 1);
                uint32_t end = begin + m_entry_size;
                return m_table[end];
            }

            uint64_t freq(uint32_t i) const {
                returm m_freq[i];
            }

            uint64_t metric(uint32_t i) const {
                returm m_metric[i];
            }

            uint32_t special_cases() const {
                return m_reserved;
            }

            uint32_t const* get(uint32_t i) const {
                uint32_t begin = i * (m_entry_size + 1);
                return &m_table[begin];
            }

            std::string entry_string(uint32_t i) const {
                uint32_t begin = i * (m_entry_size + 1);
                uint32_t end = begin + m_entry_size;
                uint32_t const* entry = &m_table[begin];
                auto esize =  m_table[end];
                std::string estr = "[";
                for(size_t i=0;i<esize-1;i++) {
                    estr += std::to_string(entry[i]) + ",";
                }
                return estr + std::to_string(entry[esize-1]) + "]";
            }

        private:
            uint32_t m_pos;
            uint32_t m_size;
            uint32_t m_reserved;
            uint32_t m_capacity;
            uint32_t m_entry_size;
            std::vector<uint32_t> m_table;
            std::vector<uint64_t> m_freq;
            std::vector<uint64_t> m_metric;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;


        };

        dictionary()
            : m_capacity(0)
            , m_entry_size(0)
        {}

        // void optimize() {
        //     m_T.resize(65536, std::vector<uint32_t>(17, 0));
        //     int sum = 0;
        //     for (int i = 0; i < 65536; ++i) {
        //         for (int j = 0; j < 17; ++j) {
        //             m_T[i][j] = m_table[sum + j];
        //         }
        //         sum += 17;
        //     }
        // }

        uint32_t copy(uint32_t i, uint32_t* out) const {

            assert(i < 65536);

            // uint32_t begin = i * (m_entry_size + 1);
            // uint32_t end = begin + m_entry_size;
            // uint32_t size = m_table[end];
            // uint32_t const* ptr = &m_table[begin];
            // // std::copy(ptr, ptr + size, out);
            // std::copy(ptr, &m_table[end], out);
            // return size;

            // // std::copy(&m_T[i][0], &m_T[i][16], out);
            // // return m_T[i][16];


            // 1.78 ns x int
            uint32_t begin = i * 17;
            uint32_t end = begin + 16;
            uint32_t size = m_table[end];
            uint32_t const* ptr = &m_table[begin];
            memcpy(out, ptr, 64);

            // // 1.83 ns x int
            // uint32_t begin = i * (m_entry_size + 1); // i * 17
            // uint32_t end = begin + m_entry_size; // begin + 16
            // uint32_t size = m_table[end];
            // uint32_t const* ptr = &m_table[begin];
            // memcpy(out, ptr, 64);

            // // 1.95 ns x int
            // uint32_t begin = i * (m_entry_size + 1); // i * 17
            // uint32_t end = begin + m_entry_size; // begin + 16
            // uint32_t size = m_table[end];
            // uint32_t const* ptr = &m_table[begin];
            // std::copy(ptr, ptr + 16, out);

            return size;
        }

        size_t capacity() const {
            return m_capacity;
        }

        size_t entry_size() const {
            return m_entry_size;
        }

        void swap(dictionary& other) {
            std::swap(m_capacity, other.m_capacity);
            std::swap(m_entry_size, other.m_entry_size);
            m_table.swap(other.m_table);
        }

        template<typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_capacity, "m_capacity")
                (m_entry_size, "m_entry_size")
                (m_table, "m_table")
                ;
        }

    private:
        uint32_t m_capacity;
        uint32_t m_entry_size;
        succinct::mapper::mappable_vector<uint32_t> m_table;

        // std::vector<std::vector<uint32_t>> m_T;
    };

}
