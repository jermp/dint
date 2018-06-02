#pragma once

#include <unordered_map>
#include <fstream>

#include <boost/format.hpp>
#include <succinct/mappable_vector.hpp>

#include "dint_configuration.hpp"
#include "model_build_utils.hpp"
#include "hash_utils.hpp"
#include "util.hpp"

namespace ds2i {

    template<uint32_t t_num_entries,
             uint32_t t_max_entry_size>
    struct dictionary
    {
        static_assert(is_power_of_two(t_max_entry_size));
        static const uint32_t num_entries = t_num_entries;
        static const uint32_t max_entry_size = t_max_entry_size;
        static const uint32_t invalid_index = uint32_t(-1);
        static const uint32_t reserved = 6; // 1 for exceptions
                                            // 5 for runs

        struct builder {
            static const uint32_t invalid_index = dictionary::invalid_index;
            static const uint32_t reserved = dictionary::reserved;

            builder()
                : m_pos(0)
                , m_size(reserved)

                , m_final_bpi(0.0)
                , m_total_coverage(0.0)
                , m_total_integers(0)

                , m_reserved(reserved)
                , m_capacity(0)
                , m_table(0, 0)
            {}

            void init(uint64_t total_integers = 0, std::string type = "NULL") {
                m_type = type;
                m_pos = reserved * (t_max_entry_size + 1);
                m_size = reserved;
                m_capacity = t_num_entries;

                // at the beginning, everything is exception(al)
                m_final_bpi = constants::initial_bpi;
                m_total_coverage = 0.0;
                m_total_integers = total_integers;

                m_table.resize(t_num_entries * (t_max_entry_size + 1), 0);
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
                if (ofs) {
                    write(ofs);
                    return true;
                }
                return false;
            }

            void load(std::ifstream& dictionary_file) {
                dictionary_file.read(reinterpret_cast<char*>(&m_capacity), sizeof(uint32_t));
                init();
                dictionary_file.read(reinterpret_cast<char*>(m_table.data()), m_table.size() * sizeof(uint32_t));
                m_size = m_capacity;
            }

            bool full() {
                return m_size == m_capacity;
            }

            bool append(uint32_t const* entry, uint32_t entry_size, uint64_t freq)
            {
                if (full()) return false;

                assert(m_pos % (t_max_entry_size + 1) == 0);
                std::copy(entry, entry + entry_size, &m_table[m_pos]);
                m_pos += t_max_entry_size + 1;
                m_table[m_pos - 1] = entry_size;
                ++m_size;

                // logging
                double cost_saving = compute_saving(entry_size, freq, m_total_integers);
                m_total_coverage += freq * entry_size * 100.0 / m_total_integers;
                m_final_bpi -= cost_saving;
                if (m_size % 500 == 0) {
                    logger() << "entries in dictionary " << m_size << "/" << num_entries << std::endl;
                    logger() << "current bits x integer: " << m_final_bpi << std::endl;
                    logger() << "covering " << m_total_coverage << "% of integers" << std::endl;
                }

                return true;
            }

            void prepare_for_encoding() {
                std::vector<uint32_t> run(256, 1);
                uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
                uint32_t i = 1;
                for (uint32_t n = 256; n != max_entry_size(); n /= 2, ++i) {
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
                m_map.swap(other.m_map);
            }

            uint32_t size(uint32_t i) const {
                uint32_t begin = i * (t_max_entry_size + 1);
                uint32_t end = begin + t_max_entry_size;
                return m_table[end];
            }

            double bpi() const {
                return m_final_bpi;
            }

            double coverage() const {
                return m_total_coverage;
            }

            uint32_t const* get(uint32_t i) const {
                uint32_t begin = i * (t_max_entry_size + 1);
                return &m_table[begin];
            }

            // std::string entry_string(uint32_t i) const {
            //     if(i == 0) return "[exception 16bit]";
            //     if(i == 1) return "[exception 32bit]";
            //     if(i == 2) return "[1]*256";
            //     if(i == 3) return "[1]*128";
            //     if(i == 4) return "[1]*64";
            //     if(i == 5) return "[1]*32";

            //     uint32_t begin = i * (t_max_entry_size + 1);
            //     uint32_t end = begin + t_max_entry_size;
            //     uint32_t const* entry = &m_table[begin];
            //     auto esize =  m_table[end];
            //     std::string estr = "[";
            //     for(size_t i=0;i<esize-1;i++) {
            //         estr += std::to_string(entry[i]) + ",";
            //     }
            //     return estr + std::to_string(entry[esize-1]) + "]";
            // }

            std::string const& type() const {
                return m_type;
            }

            // void print() {
            //     logger() << "type = " << type();
            //     logger() << "     size = " << m_size;
            //     logger() << "     special_cases = " << m_reserved;
            //     logger() << "     max_entry_size = " << t_max_entry_size;

            //     std::vector<uint32_t> len_stats(257);
            //     for(size_t i=0;i<m_size;i++) {
            //         len_stats[size(i)]++;
            //     }
            //     logger() << " LEN DIST = ";
            //     boost::format fmt("\t   len = %1$3d count = %2$5d percent = %3$4.2f");
            //     for(size_t i=0;i<len_stats.size();i++) {
            //         if(len_stats[i] != 0) {
            //             logger() << fmt % i % len_stats[i] % (double(len_stats[i]) / double(m_size) * 100);
            //         }
            //     }
            //     logger() << " CONTENT = ";
            //     boost::format fmt2("\t   code = %1$6d freq = %2$10d entry = %3%");
            //     for(size_t i=0;i<m_size;i++) {
            //         logger() << fmt2 % i % freq(i) %  entry_string(i);
            //     }
            // }

        private:
            std::string m_type;
            uint32_t m_pos;
            uint32_t m_size;

            double m_final_bpi;
            double m_total_coverage;
            uint64_t m_total_integers;

            uint32_t m_reserved;
            uint32_t m_capacity;
            std::vector<uint32_t> m_table;

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
            memcpy(out, ptr, t_max_entry_size * sizeof(uint32_t));
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

// #pragma once

// #include <succinct/mappable_vector.hpp>
// #include <fstream>

// #include "hash_utils.hpp"
// #include "util.hpp"

// namespace ds2i {

//     template<uint32_t t_num_entries,
//              uint32_t t_max_entry_size>
//     struct dictionary {

//         static_assert(is_power_of_two(t_max_entry_size));
//         static const uint32_t num_entries = t_num_entries;
//         static const uint32_t max_entry_size = t_max_entry_size;
//         static const uint32_t invalid_index = uint32_t(-1);
//         static const uint32_t reserved = 6; // 1 for exceptions
//                                             // 5 for runs

//         struct builder {

//             builder()
//                 : m_pos(0)
//                 , m_size(reserved)
//                 , m_capacity(0)
//                 , m_max_entry_size(0)
//                 , m_table(0, 0)
//             {}

//             void init(uint32_t capacity, uint32_t entry_size) {
//                 if (!is_power_of_two(entry_size)) {
//                     throw std::runtime_error("entry size must be a power of 2");
//                 }
//                 m_pos = reserved * (entry_size + 1);
//                 m_size = reserved;
//                 m_capacity = capacity;
//                 m_max_entry_size = entry_size;

//                 // XXX: Giulio
//                 /*
//                     We prefill the table with 1s to let the
//                     decoder skip whenever we have a run.
//                     In fact, whenever we have to copy en entry
//                     from the dictionary of size t, we copy to the output:
//                     x ..... x 1 .... 1
//                     --------- --------
//                         t      unused
//                     if a run will follow, the decoder will have to just skip
//                     without performing any copy.

//                     If we encode deltas by subtracting an additional 1,
//                     we must change the default value to 0.
//                 */
//                 m_table.resize(capacity * (entry_size + 1), 0);
//             }

//             builder(uint32_t capacity, uint32_t entry_size) {
//                 init(capacity, entry_size);
//             }

//             void write(std::ofstream& dictionary_file) const {
//                 dictionary_file.write(reinterpret_cast<char const*>(&m_capacity), sizeof(uint32_t));
//                 dictionary_file.write(reinterpret_cast<char const*>(&m_max_entry_size), sizeof(uint32_t));
//                 dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), m_table.size() * sizeof(uint32_t));
//             }

//             void load(std::ifstream& dictionary_file) {
//                 uint32_t capacity, entry_size;
//                 dictionary_file.read(reinterpret_cast<char*>(&capacity), sizeof(uint32_t));
//                 dictionary_file.read(reinterpret_cast<char*>(&entry_size), sizeof(uint32_t));
//                 init(capacity, entry_size);
//                 dictionary_file.read(reinterpret_cast<char*>(m_table.data()), m_table.size() * sizeof(uint32_t));
//                 m_size = m_capacity;
//             }

//             bool full() {
//                 return m_size == m_capacity;
//             }

//             bool append(uint32_t const* entry, uint32_t entry_size) {
//                 if (full()) {
//                     return false;
//                 }
//                 assert(m_pos % (m_max_entry_size + 1) == 0);
//                 std::copy(entry, entry + entry_size, &m_table[m_pos]);
//                 m_pos += m_max_entry_size + 1;
//                 m_table[m_pos - 1] = entry_size;
//                 ++m_size;
//                 return true;
//             }

//             void prepare_for_encoding() {
//                 logger() << "building mapping for encoding " << std::endl;
//                 std::vector<uint32_t> run(256, 0);
//                 uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
//                 uint32_t i = 1;
//                 for (uint32_t n = 256; n != 8; n /= 2, ++i) {
//                     uint64_t hash = hash_bytes64(byte_range(ptr, ptr + n * sizeof(uint32_t)));
//                     m_map[hash] = i;
//                 }
//                 for (; i < capacity(); ++i) {
//                     uint8_t const* ptr = reinterpret_cast<uint8_t const*>(get(i));
//                     uint32_t entry_size = size(i);
//                     uint64_t hash = hash_bytes64(byte_range(ptr, ptr + entry_size * sizeof(uint32_t)));
//                     m_map[hash] = i;
//                 }
//             }

//             // Giulio: return the index of the pattern if found in the table,
//             // otherwise return the invalid_index value
//             uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
//             {
//                 uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
//                 uint8_t const* e = b + entry_size * sizeof(uint32_t);
//                 uint64_t hash = hash_bytes64(byte_range(b, e));
//                 auto it = m_map.find(hash);
//                 if (it != m_map.end()) {
//                     assert((*it).second <= capacity());
//                     return (*it).second;
//                 }
//                 return invalid_index;
//             }

//             size_t capacity() const {
//                 return m_capacity;
//             }

//             size_t max_entry_size() const {
//                 return m_max_entry_size;
//             }

//             void build(dictionary& dict) {
//                 std::swap(m_capacity, dict.m_capacity);
//                 dict.m_table.steal(m_table);
//                 builder().swap(*this);
//             }

//             void swap(builder& other) {
//                 std::swap(m_pos, other.m_pos);
//                 std::swap(m_size, other.m_size);
//                 std::swap(m_capacity, other.m_capacity);
//                 std::swap(m_max_entry_size, other.m_max_entry_size);
//                 m_table.swap(other.m_table);
//                 m_map.swap(other.m_map);
//             }

//         private:
//             uint32_t m_pos;
//             uint32_t m_size;
//             uint32_t m_capacity;
//             uint32_t m_max_entry_size;
//             std::vector<uint32_t> m_table;

//             // map from hash codes to table indexes, used during encoding
//             std::unordered_map<uint64_t, uint32_t> m_map;

//             uint32_t const* get(uint32_t i) const {
//                 uint32_t begin = i * (m_max_entry_size + 1);
//                 return &m_table[begin];
//             }

//             uint32_t size(uint32_t i) const {
//                 uint32_t begin = i * (m_max_entry_size + 1);
//                 uint32_t end = begin + m_max_entry_size;
//                 return m_table[end];
//             }
//         };

//         dictionary()
//             : m_capacity(0)
//         {}

//         // void print() {
//         //     uint64_t sum = 0;
//         //     for (int i = 0; i < 65536; ++i) {
//         //         for (int j = 0; j < 9; ++j) {
//         //             std::cout << m_table[sum + j] << " ";
//         //         }
//         //         std::cout << std::endl;
//         //         sum += 9;
//         //     }
//         // }

//         uint32_t copy(uint32_t i, uint32_t* out) const
//         {
//             assert(i < capacity());
//             uint32_t begin = i * 17;
//             uint32_t end = begin + 16;
//             uint32_t size = m_table[end];
//             uint32_t const* ptr = &m_table[begin];
//             memcpy(out, ptr, 64);

//             // uint32_t begin = i * 9;
//             // uint32_t end = begin + 8;
//             // uint32_t size = m_table[end];
//             // uint32_t const* ptr = &m_table[begin];
//             // memcpy(out, ptr, 32);

//             return size;
//         }

//         size_t capacity() const {
//             return m_capacity;
//         }

//         void swap(dictionary& other) {
//             std::swap(m_capacity, other.m_capacity);
//             m_table.swap(other.m_table);
//         }

//         template<typename Visitor>
//         void map(Visitor& visit)
//         {
//             visit
//                 (m_capacity, "m_capacity")
//                 (m_table, "m_table")
//                 ;
//         }

//     private:
//         uint32_t m_capacity;
//         succinct::mapper::mappable_vector<uint32_t> m_table;
//     };

//     // struct dictionary {

//     //     static const uint32_t invalid_index = uint32_t(-1);
//     //     static const uint32_t reserved = 6; // 1 for exceptions
//     //                                         // 5 for runs

//     //     struct builder {

//     //         builder()
//     //         {}

//     //         void init(uint32_t capacity, uint32_t max_entry_size) {
//     //             m_size = reserved;
//     //             m_capacity = capacity;
//     //             m_max_entry_size = max_entry_size;
//     //             m_offsets.reserve(capacity);
//     //             for (uint32_t i = 0; i < reserved; ++i) { // unused offsets
//     //                 m_offsets.push_back(0);
//     //             }
//     //         }

//     //         builder(uint32_t capacity, uint32_t max_entry_size) {
//     //             init(capacity, max_entry_size);
//     //         }

//     //         void write(std::ofstream& dictionary_file) const {
//     //             uint32_t offsets_size = m_offsets.size();
//     //             uint32_t table_size = m_table.size();
//     //             dictionary_file.write(reinterpret_cast<char const*>(&m_size), sizeof(uint32_t));
//     //             dictionary_file.write(reinterpret_cast<char const*>(&m_max_entry_size), sizeof(uint32_t));
//     //             dictionary_file.write(reinterpret_cast<char const*>(&offsets_size), sizeof(uint32_t));
//     //             dictionary_file.write(reinterpret_cast<char const*>(&table_size), sizeof(uint32_t));
//     //             dictionary_file.write(reinterpret_cast<char const*>(m_offsets.data()), m_offsets.size() * sizeof(uint32_t));
//     //             dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), m_table.size() * sizeof(uint32_t));
//     //         }

//     //         void load(std::ifstream& dictionary_file) {
//     //             uint32_t num_entries, max_entry_size, offsets_size, table_size;
//     //             dictionary_file.read(reinterpret_cast<char*>(&num_entries), sizeof(uint32_t));
//     //             dictionary_file.read(reinterpret_cast<char*>(&max_entry_size), sizeof(uint32_t));
//     //             dictionary_file.read(reinterpret_cast<char*>(&offsets_size), sizeof(uint32_t));
//     //             dictionary_file.read(reinterpret_cast<char*>(&table_size), sizeof(uint32_t));
//     //             m_size = num_entries;
//     //             m_capacity = num_entries;
//     //             m_max_entry_size = max_entry_size;
//     //             m_table.resize(table_size);
//     //             m_offsets.resize(offsets_size);
//     //             dictionary_file.read(reinterpret_cast<char*>(m_offsets.data()), offsets_size * sizeof(uint32_t));
//     //             dictionary_file.read(reinterpret_cast<char*>(m_table.data()), table_size * sizeof(uint32_t));
//     //         }

//     //         bool full() {
//     //             return m_size == m_capacity;
//     //         }

//     //         bool append(uint32_t const* entry, uint32_t entry_size) {

//     //             if (entry_size > m_max_entry_size) {
//     //                 throw std::runtime_error("entry not allowed");
//     //             }

//     //             if (full()) {
//     //                 return false;
//     //             }

//     //             m_offsets.push_back(m_table.size());
//     //             // std::cout << "m_offsets.size() = " << m_offsets.size() << std::endl;

//     //             // [size, entry]
//     //             // std::cout << "pusing an entry of size " << entry_size << std::endl;
//     //             m_table.push_back(entry_size);
//     //             for (uint32_t i = 0; i < entry_size; ++i, ++entry) {
//     //                 m_table.push_back(*entry);
//     //                 // std::cout << *entry << " ";
//     //             }
//     //             // std::cout << "table size = " << m_table.size() << std::endl;

//     //             ++m_size;
//     //             return true;
//     //         }

//     //         void prepare_for_encoding() {
//     //             logger() << m_size << " entries in the dictionary" << std::endl;
//     //             // logger() << "building mapping for encoding " << std::endl;
//     //             std::vector<uint32_t> run(256, 0);
//     //             uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
//     //             uint32_t i = 1;
//     //             for (uint32_t n = 256; n != 8; n /= 2, ++i) {
//     //                 uint64_t hash = hash_bytes64(byte_range(ptr, ptr + n * sizeof(uint32_t)));
//     //                 m_map[hash] = i;
//     //             }
//     //             for (; i < m_size; ++i) {
//     //                 uint8_t const* ptr = reinterpret_cast<uint8_t const*>(get(i));
//     //                 uint32_t entry_size = size(i);

//     //                 // std::cout << "entry_size = " << entry_size << std::endl;
//     //                 // uint32_t const* p = get(i);
//     //                 // for (uint32_t k = 0; k < entry_size; ++k) {
//     //                 //     std::cout << *p << " ";
//     //                 //     ++p;
//     //                 // }
//     //                 // std::cout << std::endl;

//     //                 uint64_t hash = hash_bytes64(byte_range(ptr, ptr + entry_size * sizeof(uint32_t)));
//     //                 m_map[hash] = i;
//     //             }
//     //         }

//     //         // Giulio: return the index of the pattern if found in the table,
//     //         // otherwise return the invalid_index value
//     //         uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
//     //         {
//     //             uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
//     //             uint8_t const* e = b + entry_size * sizeof(uint32_t);
//     //             uint64_t hash = hash_bytes64(byte_range(b, e));
//     //             auto it = m_map.find(hash);
//     //             if (it != m_map.end()) {
//     //                 assert((*it).second <= capacity());
//     //                 return (*it).second;
//     //             }
//     //             return invalid_index;
//     //         }

//     //         uint32_t capacity() const {
//     //             return m_capacity;
//     //         }

//     //         uint32_t size() const {
//     //             return m_size;
//     //         }

//     //         uint32_t max_entry_size() const {
//     //             return m_max_entry_size;
//     //         }

//     //         void build(dictionary& dict) {
//     //             std::swap(m_size, dict.m_size);
//     //             dict.m_offsets.steal(m_offsets);
//     //             dict.m_table.steal(m_table);
//     //             builder().swap(*this);
//     //         }

//     //         void swap(builder& other) {
//     //             std::swap(m_size, other.m_size);
//     //             std::swap(m_capacity, other.m_capacity);
//     //             std::swap(m_max_entry_size, other.m_max_entry_size);
//     //             m_offsets.swap(other.m_offsets);
//     //             m_table.swap(other.m_table);
//     //             m_map.swap(other.m_map);
//     //         }

//     //         // void print() {
//     //         //     for (uint32_t i = reserved; i < capacity(); ++i) {
//     //         //         uint32_t entry_size = size(i);
//     //         //         std::cout << "entry_size = " << entry_size << std::endl;
//     //         //         uint32_t const* p = get(i);
//     //         //         for (uint32_t k = 0; k < entry_size; ++k) {
//     //         //             std::cout << *p << " ";
//     //         //             ++p;
//     //         //         }
//     //         //         std::cout << std::endl;
//     //         //     }
//     //         // }

//     //     private:
//     //         uint32_t m_size;
//     //         uint32_t m_capacity;
//     //         uint32_t m_max_entry_size;
//     //         std::vector<uint32_t> m_offsets;
//     //         std::vector<uint32_t> m_table;

//     //         // map from hash codes to table indexes, used during encoding
//     //         std::unordered_map<uint64_t, uint32_t> m_map;

//     //         uint32_t const* get(uint32_t i) const {
//     //             assert(i < size());
//     //             uint32_t offset = m_offsets[i];
//     //             return &m_table[offset + 1]; // skip size
//     //         }

//     //         uint32_t size(uint32_t i) const {
//     //             assert(i < size());
//     //             return m_table[m_offsets[i]];
//     //         }
//     //     };

//     //     dictionary()
//     //         : m_size(0)
//     //     {}

//     //     uint32_t copy(uint32_t i, uint32_t* out) const
//     //     {
//     //         assert(i < size());
//     //         uint32_t offset = m_offsets[i];
//     //         uint32_t size = m_table[offset];
//     //         uint32_t const* ptr = &m_table[offset + 1];
//     //         // memcpy(out, ptr, size * 4);
//     //         memcpy(out, ptr, 64);
//     //         return size;
//     //     }

//     //     size_t size() const {
//     //         return m_size;
//     //     }

//     //     void swap(dictionary& other) {
//     //         std::swap(m_size, other.m_size);
//     //         m_offsets.swap(other.m_offsets);
//     //         m_table.swap(other.m_table);
//     //     }

//     //     template<typename Visitor>
//     //     void map(Visitor& visit)
//     //     {
//     //         visit
//     //             (m_size, "m_size")
//     //             (m_offsets, "m_offsets")
//     //             (m_table, "m_table")
//     //             ;
//     //     }

//     // private:
//     //     uint32_t m_size;
//     //     succinct::mapper::mappable_vector<uint32_t> m_offsets;
//     //     succinct::mapper::mappable_vector<uint32_t> m_table;
//     // };
// }
