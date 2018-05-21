#pragma once

#include <succinct/mappable_vector.hpp>
#include <fstream>

#include "hash_utils.hpp"
#include "util.hpp"

namespace ds2i {

    struct dictionary {

        static const uint32_t invalid_index = uint32_t(-1);
        static const uint32_t reserved = 6; // 1 for exceptions
                                            // 5 for runs

        struct builder {

            builder()
                : m_pos(0)
                , m_size(reserved)
                , m_capacity(0)
                , m_entry_size(0)
                , m_table(0, 1)
            {}

            void init(uint32_t capacity, uint32_t entry_size) {
                if (!is_power_of_two(entry_size)) {
                    throw std::runtime_error("entry size must be a power of 2");
                }
                m_pos = reserved * (entry_size + 1);
                m_size = reserved;
                m_capacity = capacity;
                m_entry_size = entry_size;

                // XXX: Giulio
                /*
                    We prefill the table with 1s to let the
                    decoder skip whenever we have a run.
                    In fact, whenever we have to copy en entry
                    from the dictionary of size t, we copy to the output:
                    x ..... x 1 .... 1
                    --------- --------
                        t      unused
                    if a run will follow, the decoder will have to just skip
                    without performing any copy.

                    If we encode deltas by subtracting an additional 1,
                    we must change the default value to 0.
                */
                m_table.resize(capacity * (entry_size + 1), 1);
            }

            builder(uint32_t capacity, uint32_t entry_size) {
                init(capacity, entry_size);
            }

            void write(std::ofstream& dictionary_file) const {
                dictionary_file.write(reinterpret_cast<char const*>(&m_capacity), sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(&m_entry_size), sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), m_table.size() * sizeof(uint32_t));
            }

            void load(std::ifstream& dictionary_file) {
                uint32_t capacity, entry_size;
                dictionary_file.read(reinterpret_cast<char*>(&capacity), sizeof(uint32_t));
                dictionary_file.read(reinterpret_cast<char*>(&entry_size), sizeof(uint32_t));
                init(capacity, entry_size);
                dictionary_file.read(reinterpret_cast<char*>(m_table.data()), m_table.size() * sizeof(uint32_t));
                m_size = m_capacity;
            }

            bool full() {
                return m_size == m_capacity;
            }

            bool append(uint32_t const* entry, uint32_t entry_size) {
                if (full()) {
                    return false;
                }
                assert(m_pos % (m_entry_size + 1) == 0);
                std::copy(entry, entry + entry_size, &m_table[m_pos]);
                m_pos += m_entry_size + 1;
                m_table[m_pos - 1] = entry_size;
                ++m_size;
                return true;
            }

            void prepare_for_encoding() {
                logger() << "building mapping for encoding " << std::endl;
                std::vector<uint32_t> run(256, 1);
                uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
                uint32_t i = 0;
                for (uint32_t n = 256; n != 8; n /= 2, ++i) {
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

        private:
            uint32_t m_pos;
            uint32_t m_size;
            uint32_t m_capacity;
            uint32_t m_entry_size;
            std::vector<uint32_t> m_table;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;

            uint32_t const* get(uint32_t i) const {
                uint32_t begin = i * (m_entry_size + 1);
                return &m_table[begin];
            }

            uint32_t size(uint32_t i) const {
                uint32_t begin = i * (m_entry_size + 1);
                uint32_t end = begin + m_entry_size;
                return m_table[end];
            }
        };

        dictionary()
            : m_capacity(0)
            , m_entry_size(0)
        {}

        // void print() {
        //     uint64_t sum = 0;
        //     for (int i = 0; i < 65536; ++i) {
        //         for (int j = 0; j < 9; ++j) {
        //             std::cout << m_table[sum + j] << " ";
        //         }
        //         std::cout << std::endl;
        //         sum += 9;
        //     }
        // }

        uint32_t copy(uint32_t i, uint32_t* out) const {

            assert(i < 65536);


            // 1.78 ns x int
            // uint32_t begin = i * 17;
            // uint32_t end = begin + 16;
            // uint32_t size = m_table[end];
            // uint32_t const* ptr = &m_table[begin];
            // memcpy(out, ptr, 64);

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

            uint32_t begin = i * 17;
            uint32_t end = begin + 16;
            uint32_t size = m_table[end];
            uint32_t const* ptr = &m_table[begin];
            memcpy(out, ptr, 64);

            // uint32_t begin = i * 9;
            // uint32_t end = begin + 8;
            // uint32_t size = m_table[end];
            // uint32_t const* ptr = &m_table[begin];
            // memcpy(out, ptr, 32);

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
    };

}
