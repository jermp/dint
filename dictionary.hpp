#pragma once

#include <unordered_map>

#include <succinct/mappable_vector.hpp>
#include <fstream>

#include "hash_utils.hpp"
#include "util.hpp"

#include <boost/format.hpp>

namespace ds2i {

    template<uint32_t t_max_entry_size>
    struct dictionary {
        static_assert(is_power_of_two(t_max_entry_size));

        static const uint32_t invalid_index = uint32_t(-1);
        static const uint32_t reserved = 6; // 2 for exceptions
                                            // 4 for runs

        struct builder {
            static const uint32_t invalid_index = dictionary::invalid_index;
            static const uint32_t reserved = dictionary::reserved;

            builder()
                : m_pos(0)
                , m_size(reserved)
                , m_reserved(reserved)
                , m_capacity(0)
                , m_table(0, 0)
            {}

            void init(uint32_t capacity,std::string type = "NULL") {
                m_type = type;
                m_pos = reserved * (t_max_entry_size + 1);
                m_size = reserved;
                m_capacity = capacity;

                m_table.resize(capacity * (t_max_entry_size + 1), 0);
                m_freq.resize(0);

                m_freq.push_back(0);    // exceptions u16
                m_freq.push_back(0);    // exceptions u32

                m_freq.push_back(0);    // 256
                m_freq.push_back(0);    // 128
                m_freq.push_back(0);    // 64
                m_freq.push_back(0);    // 32
            }

            builder(uint32_t capacity) {
                init(capacity);
            }

            void load_from_file(std::string dict_file) {
                std::ifstream ifs(dict_file);
                load(ifs);
            }

            void write(std::ofstream& dictionary_file) const {
                dictionary_file.write(reinterpret_cast<char const*>(&m_capacity), sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), m_table.size() * sizeof(uint32_t));
            }

            bool try_store_to_file(std::string dict_file) const {
                std::ofstream ofs(dict_file);
                if(ofs) {
                    write(ofs);
                    return true;
                }
                return false;
            }

            void load(std::ifstream& dictionary_file) {
                dictionary_file.read(reinterpret_cast<char*>(&m_capacity), sizeof(uint32_t));
                init(m_capacity);
                dictionary_file.read(reinterpret_cast<char*>(m_table.data()), m_table.size() * sizeof(uint32_t));
                m_size = m_capacity;
            }

            bool full() {
                return m_size == m_capacity;
            }

            bool append(uint32_t const* entry, uint32_t entry_size,uint64_t freq) {
                if (full()) {
                    return false;
                }
                assert(m_pos % (t_max_entry_size + 1) == 0);
                std::copy(entry, entry + entry_size, &m_table[m_pos]);
                m_pos += t_max_entry_size + 1;
                m_table[m_pos - 1] = entry_size;
                m_freq.push_back(freq);
                ++m_size;
                return true;
            }

            void prepare_for_encoding() {
                DS2I_LOG << "building mapping for encoding ";
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

            size_t max_entry_size() const {
                return t_max_entry_size;
            }

            void build(dictionary& dict) {
                std::swap(m_capacity, dict.m_capacity);
                dict.m_table.steal(m_table);
                builder().swap(*this);
            }

            void swap(builder& other) {
                std::swap(m_pos, other.m_pos);
                std::swap(m_size, other.m_size);
                std::swap(m_capacity, other.m_capacity);
                m_table.swap(other.m_table);
                m_freq.swap(other.m_freq);
                m_map.swap(other.m_map);
            }

            uint32_t size(uint32_t i) const {
                // special cases
                if(i == 0 || i == 1) {
                    return 0;
                }
                if(i == 2) return 256;
                if(i == 3) return 128;
                if(i == 4) return 64;
                if(i == 5) return 32;
                uint32_t begin = i * (t_max_entry_size + 1);
                uint32_t end = begin + t_max_entry_size;
                return m_table[end];
            }

            uint64_t freq(uint32_t i) const {
                if(i <= 5) return 0;
                return m_freq[i];
            }

            uint32_t special_cases() const {
                return m_reserved;
            }

            uint32_t const* get(uint32_t i) const {
                uint32_t begin = i * (t_max_entry_size + 1);
                return &m_table[begin];
            }

            std::string entry_string(uint32_t i) const {
                if(i == 0) return "[exception 16bit]";
                if(i == 1) return "[exception 32bit]";
                if(i == 2) return "[1]*256";
                if(i == 3) return "[1]*128";
                if(i == 4) return "[1]*64";
                if(i == 5) return "[1]*32";

                uint32_t begin = i * (t_max_entry_size + 1);
                uint32_t end = begin + t_max_entry_size;
                uint32_t const* entry = &m_table[begin];
                auto esize =  m_table[end];
                std::string estr = "[";
                for(size_t i=0;i<esize-1;i++) {
                    estr += std::to_string(entry[i]) + ",";
                }
                return estr + std::to_string(entry[esize-1]) + "]";
            }

            std::string type() const {
                return m_type;
            }

            void print() {
                DS2I_LOG << "type = " << type();
                DS2I_LOG << "     size = " << m_size;
                DS2I_LOG << "     special_cases = " << m_reserved;
                DS2I_LOG << "     max_entry_size = " << t_max_entry_size;

                std::vector<uint32_t> len_stats(257);
                for(size_t i=0;i<m_size;i++) {
                    len_stats[size(i)]++;
                }
                DS2I_LOG << " LEN DIST = ";
                boost::format fmt("\t   len = %1$3d count = %2$5d percent = %3$4.2f");
                for(size_t i=0;i<len_stats.size();i++) {
                    if(len_stats[i] != 0) {
                        DS2I_LOG << fmt % i % len_stats[i] % (double(len_stats[i]) / double(m_size) * 100);
                    }
                }
                DS2I_LOG << " CONTENT = ";
                boost::format fmt2("\t   code = %1$6d freq = %2$10d entry = %3%");
                for(size_t i=0;i<m_size;i++) {
                    DS2I_LOG << fmt2 % i % freq(i) %  entry_string(i);
                }
            }

        private:
            std::string m_type;
            uint32_t m_pos;
            uint32_t m_size;
            uint32_t m_reserved;
            uint32_t m_capacity;
            std::vector<uint32_t> m_table;
            std::vector<uint64_t> m_freq;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;

        };

        dictionary()
            : m_capacity(0)
        {}


        uint32_t copy(uint32_t i, uint32_t* out) const
        {
            assert(i < capacity());
            uint32_t begin = i * (t_max_entry_size+1);
            uint32_t end = begin + t_max_entry_size;
            uint32_t size = m_table[end];
            uint32_t const* ptr = &m_table[begin];
            memcpy(out, ptr, t_max_entry_size*sizeof(uint32_t));
            return size;
        }

        size_t capacity() const {
            return m_capacity;
        }

        void swap(dictionary& other) {
            std::swap(m_capacity, other.m_capacity);
            m_table.swap(other.m_table);
        }

        template<typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_capacity, "m_capacity")
                (m_table, "m_table")
                ;
        }

    private:
        uint32_t m_capacity;
        succinct::mapper::mappable_vector<uint32_t> m_table;
    };

}
