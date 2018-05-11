#pragma once

#include <fstream>

#include "util.hpp"

namespace ds2i {

    struct dictionary {

        struct builder {

            static const uint32_t reserved_entries = 4;

            builder()
                : m_pos(0)
                , m_num_entries(0)
                , m_entry_width(0)
                , m_table(0, 0)
            {}

            builder(uint32_t num_entries, uint32_t entry_width)
                // Giulio: take into account reserved entries for runs
                : m_pos(reserved_entries * (entry_width + 1))
                , m_num_entries(num_entries)
                , m_entry_width(entry_width)
                , m_table(num_entries * (entry_width + 1), 0)
            {}

            // Giulio: return true if there is still space left for the entry;
            // false otherwise
            bool append(uint32_t const* entry, uint32_t entry_width) {
                if (m_pos == m_table.size()) {
                    return false;
                }
                assert(m_pos % (m_entry_width + 1) == 0);
                std::copy(entry, entry + entry_width, &m_table[m_pos]);
                m_pos += m_entry_width + 1;
                m_table[m_pos - 1] = entry_width;
                return true;
            }

            void build(dictionary& dict) {
                std::swap(m_num_entries, dict.num_entries);
                std::swap(m_entry_width, dict.m_entry_width);
                dict.m_table.steal(m_table);
                builder().swap(*this);
            }

            void swap(builder& other) {
                std::swap(m_pos, other.m_pos);
                std::swap(m_num_entries, other.m_num_entries);
                std::swap(m_entry_width, other.m_entry_width);
                m_table.swap(other.m_table);
            }

        private:
            uint32_t m_pos;
            uint32_t m_num_entries;
            uint32_t m_entry_width;
            std::vector<uint32_t> m_table;
        };

        dictionary()
            : m_num_entries(0)
            , m_entry_width(0)
        {}

        // Giulio: return a pointer to the beginning of the i-th entry
        uint32_t const* operator[](uint32_t i) const {
            uint32_t begin = i * (entry_width() + 1);
            return &m_table[begin];
        }

        // Giulio: return the size of the i-th entry
        uint32_t size(uint32_t i) const {
            uint32_t begin = i * (entry_width() + 1);
            uint32_t end = begin + entry_width();
            return m_table[end];
        }

        // Giulio: copy i-th entry to the output stream
        // and return the entry size as back-pointer
        uint32_t copy(uint32_t i, uint32_t* out) {
            uint32_t begin = i * (entry_width() + 1);
            uint32_t end = begin + entry_width();
            std::copy(&m_table[begin], &m_table[end - 1], out);
            return m_table[end];
        }

        size_t num_entries() const {
            return m_num_entries;
        }

        size_t entry_width() const {
            return m_entry_width;
        }

        void swap(dictionary& other) {
            std::swap(m_num_entries, other.m_num_entries);
            std::swap(m_entry_width, other.m_entry_width);
            m_table.steal(other.m_table);
        }

        template<typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_num_entries, "m_num_entries")
                (m_entry_width, "m_entry_width")
                (m_table, "m_table")
                ;
        }

    private:
        uint32_t m_num_entries;
        uint32_t m_entry_width;

        // Giulio: this table should be on the stack, that is
        // uint32_t m_table[65635][16 + 1]
        // make this so later
        succinct::mapper::mappable_vector<uint32_t> m_table;
    };

}
