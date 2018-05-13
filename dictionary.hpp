#pragma once

#include <succinct/mappable_vector.hpp>
#include <fstream>

#include "util.hpp"

namespace ds2i {

    struct dictionary {

        struct builder {

            static const uint32_t reserved_entries = 4;
            static const uint32_t invalid_index = uint32_t(-1);

            builder()
                : m_pos(0)
                , m_size(reserved_entries)
                , m_capacity(0)
                , m_entry_size(0)
                , m_table(0, 0)
            {}

            builder(uint32_t capacity, uint32_t entry_size)
                // Giulio: take into account reserved entries for runs
                : m_pos(reserved_entries * (entry_size + 1))
                , m_size(reserved_entries)
                , m_capacity(capacity)
                , m_entry_size(entry_size)
                , m_table(capacity * (entry_size + 1), 0)
            {
                std::cout << "m_capacity = " << m_capacity << std::endl;
            }

            bool full() {
                return m_size == m_capacity;
            }

            // Giulio: return true if there is still space left for the entry;
            // false otherwise
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
            }

        private:
            uint32_t m_pos;
            uint32_t m_size;
            uint32_t m_capacity;
            uint32_t m_entry_size;
            std::vector<uint32_t> m_table;
        };

        dictionary()
            : m_capacity(0)
            , m_entry_size(0)
        {}

        void build_mapping() {
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

        // Giulio: return a pointer to the beginning of the i-th entry
        uint32_t const* get(uint32_t i) const {
            uint32_t begin = i * (m_entry_size + 1);
            return &m_table[begin];
        }

        // Giulio: return the index of the pattern if found in the dictionary,
        // invalid_index otherwise
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

        // Giulio: return the size of the i-th entry
        uint32_t size(uint32_t i) const {
            uint32_t begin = i * (m_entry_size + 1);
            uint32_t end = begin + m_entry_size;
            return m_table[end];
        }

        // Giulio: copy i-th entry to the output stream
        // and return the entry size as back-pointer
        uint32_t copy(uint32_t i, uint32_t* out) {
            uint32_t begin = i * (m_entry_size + 1);
            uint32_t end = begin + m_entry_size;
            std::copy(&m_table[begin], &m_table[end - 1], out);
            return m_table[end];
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
            m_map.swap(other.m_map);
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

        // Giulio: this table should be on the stack, that is
        // uint32_t m_table[65635][16 + 1]
        // make this so later
        succinct::mapper::mappable_vector<uint32_t> m_table;

        // map from hash codes to table indexes, used during encoding
        std::unordered_map<uint64_t, uint32_t> m_map;
    };

}
