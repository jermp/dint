#pragma once

#include <unordered_map>
#include <unordered_set>
#include <fstream>

#include <boost/format.hpp>
#include <succinct/mappable_vector.hpp>

#include "dint_configuration.hpp"
#include "model_build_utils.hpp"
#include "hash_utils.hpp"
#include "util.hpp"

#include <x86intrin.h>

namespace ds2i {

    // template<uint32_t t_num_entries,
    //          uint32_t t_max_entry_size>
    // struct dictionary
    // {
    //     // static_assert(is_power_of_two(t_max_entry_size));
    //     static const uint32_t num_entries = t_num_entries;
    //     static const uint32_t max_entry_size = t_max_entry_size;
    //     static const uint32_t invalid_index = uint32_t(-1);
    //     static const uint32_t reserved = 6;

    //                                         // 6 reserved entries:
    //                                         // 1 for exception header
    //                                         // 5 for runs

    //                                         // 10 reserved entries:
    //                                         // 1 for exception header
    //                                         // 1 for exception 1 byte
    //                                         // 1 for exception 2 bytes
    //                                         // 1 for exception 3 bytes
    //                                         // 1 for exception 4 bytes
    //                                         // 5 for runs

    //     struct builder
    //     {
    //         static const uint32_t num_entries = dictionary::num_entries;
    //         static const uint32_t max_entry_size = dictionary::max_entry_size;
    //         static const uint32_t invalid_index = dictionary::invalid_index;
    //         static const uint32_t reserved = dictionary::reserved;

    //         uint64_t codewords = 0;
    //         uint64_t exceptions = 0;

    //         builder()
    //             : m_pos(0)
    //             , m_size(reserved)
    //             , m_final_bpi(0.0)
    //             , m_total_coverage(0.0)
    //             , m_total_integers(0)
    //             , m_table(0, 0)
    //         {}

    //         void init(uint64_t total_integers = 0) {
    //             m_pos = reserved * (max_entry_size + 1);
    //             m_size = reserved;

    //             // at the beginning, everything is exception(al)
    //             m_final_bpi = constants::initial_bpi;
    //             m_total_coverage = 0.0;
    //             m_total_integers = total_integers;

    //             m_table.resize(num_entries * (max_entry_size + 1), 0);
    //             m_freqs.reserve(num_entries);
    //         }

    //         void load_from_file(std::string dict_file) {
    //             std::ifstream ifs(dict_file);
    //             load(ifs);
    //         }

    //         bool try_store_to_file(std::string dict_file) const {
    //             std::ofstream ofs(dict_file);
    //             if (ofs) {
    //                 write(ofs);
    //                 return true;
    //             }
    //             return false;
    //         }

    //         void write(std::ofstream& dictionary_file) const {
    //             dictionary_file.write(reinterpret_cast<char const*>(&m_size), sizeof(uint32_t));
    //             dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), m_table.size() * sizeof(uint32_t));
    //         }

    //         void load(std::ifstream& dictionary_file) {
    //             uint32_t size = 0;
    //             dictionary_file.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    //             init();
    //             dictionary_file.read(reinterpret_cast<char*>(m_table.data()), m_table.size() * sizeof(uint32_t));
    //             m_size = size;
    //         }

    //         bool full() {
    //             return m_size == num_entries;
    //         }

    //         bool append(uint32_t const* entry, uint32_t entry_size, uint64_t freq)
    //         {
    //             if (full()) return false;

    //             assert(m_pos % (max_entry_size + 1) == 0);
    //             std::copy(entry, entry + entry_size, &m_table[m_pos]);
    //             m_pos += max_entry_size + 1;
    //             m_table[m_pos - 1] = entry_size;

    //             m_freqs.emplace_back(m_size, freq);

    //             ++m_size;

    //             // logging
    //             double cost_saving = compute_saving(entry_size, freq, m_total_integers);
    //             m_total_coverage += freq * entry_size * 100.0 / m_total_integers;
    //             m_final_bpi -= cost_saving;
    //             if (m_size % 5000 == 0) {
    //                 logger() << "entries in dictionary " << m_size << "/" << num_entries << std::endl;
    //                 logger() << "current bits x integer: " << m_final_bpi << std::endl;
    //                 logger() << "covering " << m_total_coverage << "% of integers" << std::endl;
    //             }

    //             return true;
    //         }

    //         void prepare_for_encoding() {
    //             logger() << size() << " entries in the dictionary" << std::endl;
    //             std::vector<uint32_t> run(256, 0);
    //             uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
    //             uint32_t i = 1;
    //             for (uint32_t n = 256; n >= max_entry_size; n /= 2, ++i) {
    //                 uint64_t hash = hash_bytes64(byte_range(ptr, ptr + n * sizeof(uint32_t)));
    //                 m_map[hash] = i;
    //             }
    //             for (; i < size(); ++i) {
    //                 uint8_t const* ptr = reinterpret_cast<uint8_t const*>(get(i));
    //                 uint32_t entry_size = size(i);
    //                 uint64_t hash = hash_bytes64(byte_range(ptr, ptr + entry_size * sizeof(uint32_t)));
    //                 m_map[hash] = i;
    //             }
    //         }

    //         // Giulio: return the index of the pattern if found in the table,
    //         // otherwise return the invalid_index value
    //         uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
    //         {
    //             uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
    //             uint8_t const* e = b + entry_size * sizeof(uint32_t);
    //             uint64_t hash = hash_bytes64(byte_range(b, e));
    //             auto it = m_map.find(hash);
    //             if (it != m_map.end()) {
    //                 assert((*it).second <= num_entries);
    //                 return (*it).second;
    //             }
    //             return invalid_index;
    //         }

    //         void build(dictionary& dict) {
    //             dict.m_table.steal(m_table);
    //             builder().swap(*this);
    //         }

    //         void swap(builder& other) {
    //             std::swap(m_pos, other.m_pos);
    //             std::swap(m_size, other.m_size);
    //             std::swap(m_final_bpi, other.m_final_bpi);
    //             std::swap(m_total_coverage, other.m_total_coverage);
    //             std::swap(m_total_integers, other.m_total_integers);
    //             m_table.swap(other.m_table);
    //             m_map.swap(other.m_map);
    //         }

    //         uint32_t size() const {
    //             return m_size;
    //         }

    //         double bpi() const {
    //             return m_final_bpi;
    //         }

    //         double coverage() const {
    //             return m_total_coverage;
    //         }

    //         static std::string type() {
    //             return "rectangular";
    //         }

    //         void print_indexes(std::string filename) {
    //             std::ofstream out(filename.c_str());
    //             std::sort(m_freqs.begin(), m_freqs.end(),
    //                 [](auto const& x, auto const& y) {
    //                     return x.second > y.second;
    //                 }
    //             );
    //             for (auto const& p: m_freqs) {
    //                 out << p.first << "\n";
    //             }
    //             out.close();
    //         }

    //         void print_entries(std::string filename)
    //         {
    //             std::ofstream out(filename.c_str());

    //             uint32_t sizes[5] = {0, 0, 0, 0, 0};

    //             std::sort(m_freqs.begin(), m_freqs.end(),
    //                 [](auto const& x, auto const& y) {
    //                     return x.second > y.second;
    //                 }
    //             );

    //             for (auto const& p: m_freqs)
    //             {
    //                 uint32_t i = p.first;
    //                 uint32_t const* ptr = get(i);
    //                 uint32_t n = size(i);
    //                 out << i << " freq: ";
    //                 out << p.second << " entry: ";
    //                 out << "[";
    //                 // std::cout << std::setw( 6) << i << " freq: ";
    //                 // std::cout << std::setw(16) << p.second << " entry: ";
    //                 // std::cout << "[";
    //                 for (uint32_t k = 0; k < n; ++k) {
    //                     // std::cout << *ptr;
    //                     out << *ptr;
    //                     if (k != n - 1) //std::cout << "|";
    //                         out << "|";
    //                     ++ptr;
    //                 }
    //                 // std::cout << "]" << std::endl;
    //                 out << "]\n";
    //                 sizes[ceil_log2(n)] += 1;
    //             }

    //             out << "{";
    //             for (uint32_t i = 0; i < 5; ++i) {
    //                 out << "\"" << (uint32_t(1) << i) << "\": \""
    //                     << sizes[i] << "\"";
    //                 if (i != 4) {
    //                     out << ", ";
    //                 }
    //             }
    //             out << "}";
    //             out.close();
    //         }

    //     private:
    //         std::string m_type;
    //         uint32_t m_pos;
    //         uint32_t m_size;
    //         double m_final_bpi;
    //         double m_total_coverage;
    //         uint64_t m_total_integers;
    //         std::vector<uint32_t> m_table;

    //         std::vector<std::pair<uint32_t, uint64_t>> m_freqs;

    //         // map from hash codes to table indexes, used during encoding
    //         std::unordered_map<uint64_t, uint32_t> m_map;

    //         uint32_t size(uint32_t i) const {
    //             uint32_t begin = i * (t_max_entry_size + 1);
    //             uint32_t end = begin + t_max_entry_size;
    //             return m_table[end];
    //         }

    //         uint32_t const* get(uint32_t i) const {
    //             uint32_t begin = i * (t_max_entry_size + 1);
    //             return &m_table[begin];
    //         }
    //     };

    //     dictionary()
    //     {}

    //     uint32_t copy(uint32_t i, uint32_t* out) const
    //     {
    //         assert(i < num_entries);
    //         uint32_t begin = i * (max_entry_size + 1);
    //         uint32_t end = begin + max_entry_size;
    //         uint32_t size = m_table[end];
    //         uint32_t const* ptr = &m_table[begin];
    //         memcpy(out, ptr, max_entry_size * sizeof(uint32_t));


    //         // TODO: try SIMD here
    //         // __m256i tmp = _mm256_load_si256((const __m256i *) &m_table[begin]);
    //         // _mm256_storeu_si256((__m256i *) out, tmp);

    //         return size;
    //     }

    //     void swap(dictionary& other) {
    //         m_table.swap(other.m_table);
    //     }

    //     template<typename Visitor>
    //     void map(Visitor& visit)
    //     {
    //         visit
    //             (m_table, "m_table")
    //             ;
    //     }

    //     void print(uint32_t i) {
    //         uint32_t begin = i * (max_entry_size + 1);
    //         uint32_t end = begin + max_entry_size;
    //         uint32_t size = m_table[end];
    //         uint32_t const* ptr = &m_table[begin];
    //         std::cout << "[";
    //         for (uint32_t k = 0; k < size; ++k) {
    //             std::cout << *ptr;
    //             ++ptr;
    //             if (k != size - 1) std::cout << "|";
    //         }
    //         std::cout << "]" << std::endl;
    //     }

    //     void print() {

    //         uint32_t sizes[5] = {0, 0, 0, 0, 0};

    //         for (uint64_t i = reserved * (max_entry_size + 1);
    //                       i < m_table.size(); i += max_entry_size + 1) {
    //             uint32_t size = m_table[i + max_entry_size];
    //             // std::cout << "[";
    //             // for (uint32_t k = 0; k < size; ++k) {
    //             //     std::cout << m_table[i + k];
    //             //     if (k != size - 1) std::cout << "|";
    //             // }
    //             // std::cout << "]" << std::endl;

    //             // std::cout << size << std::endl;
    //             sizes[ceil_log2(size)] += 1;
    //         }

    //         for (uint32_t i = 0; i < 5; ++i) {
    //             std::cout << "entries of size " << (uint32_t(1) << i)
    //                       << ": " << sizes[i] << std::endl;
    //         }

    //     }

    // private:
    //     succinct::mapper::mappable_vector<uint32_t> m_table;
    // };

    // template<uint32_t t_num_entries,
    //          uint32_t t_max_entry_size>
    // struct dictionary
    // {
    //     static_assert(is_power_of_two(t_max_entry_size));
    //     static const uint32_t num_entries = t_num_entries;
    //     static const uint32_t max_entry_size = t_max_entry_size;
    //     static const uint32_t invalid_index = uint32_t(-1);
    //     static const uint32_t reserved = 6; // 1 for exceptions
    //                                         // 5 for runs

    //     struct builder
    //     {
    //         static const uint32_t num_entries = dictionary::num_entries;
    //         static const uint32_t max_entry_size = dictionary::max_entry_size;
    //         static const uint32_t invalid_index = dictionary::invalid_index;
    //         static const uint32_t reserved = dictionary::reserved;

    //         builder()
    //             : m_size(reserved)
    //             , m_final_bpi(0.0)
    //             , m_total_coverage(0.0)
    //             , m_total_integers(0)
    //         {}

    //         void init(uint64_t total_integers = 0) {
    //             m_size = reserved;

    //             // at the beginning, everything is exception(al)
    //             m_final_bpi = constants::initial_bpi;
    //             m_total_coverage = 0.0;
    //             m_total_integers = total_integers;

    //             m_offsets.reserve(num_entries);
    //             for (uint32_t i = 0; i < reserved; ++i) { // unused offsets
    //                 m_offsets.push_back(0);
    //             }
    //         }

    //         void load_from_file(std::string dict_file) {
    //             std::ifstream ifs(dict_file);
    //             load(ifs);
    //         }

    //         bool try_store_to_file(std::string dict_file) const {
    //             std::ofstream ofs(dict_file);
    //             if (ofs) {
    //                 write(ofs);
    //                 return true;
    //             }
    //             return false;
    //         }

    //         void write(std::ofstream& dictionary_file) const {
    //             uint32_t offsets_size = m_offsets.size();
    //             uint32_t table_size = m_table.size();
    //             dictionary_file.write(reinterpret_cast<char const*>(&m_size), sizeof(uint32_t));
    //             dictionary_file.write(reinterpret_cast<char const*>(&offsets_size), sizeof(uint32_t));
    //             dictionary_file.write(reinterpret_cast<char const*>(&table_size), sizeof(uint32_t));
    //             dictionary_file.write(reinterpret_cast<char const*>(m_offsets.data()), offsets_size * sizeof(uint32_t));
    //             dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), table_size * sizeof(uint32_t));
    //         }

    //         void load(std::ifstream& dictionary_file) {
    //             uint32_t offsets_size = 0;
    //             uint32_t table_size = 0;
    //             dictionary_file.read(reinterpret_cast<char*>(&m_size), sizeof(uint32_t));
    //             dictionary_file.read(reinterpret_cast<char*>(&offsets_size), sizeof(uint32_t));
    //             dictionary_file.read(reinterpret_cast<char*>(&table_size), sizeof(uint32_t));
    //             m_table.resize(table_size);
    //             m_offsets.resize(offsets_size);
    //             dictionary_file.read(reinterpret_cast<char*>(m_offsets.data()), offsets_size * sizeof(uint32_t));
    //             dictionary_file.read(reinterpret_cast<char*>(m_table.data()), table_size * sizeof(uint32_t));
    //         }

    //         bool full() {
    //             return m_size == num_entries;
    //         }

    //         bool append(uint32_t const* entry, uint32_t entry_size, uint64_t freq)
    //         {
    //             if (entry_size > max_entry_size) {
    //                 throw std::runtime_error("entry not allowed");
    //             }

    //             if (full()) {
    //                 return false;
    //             }

    //             m_offsets.push_back(m_table.size());
    //             m_table.push_back(entry_size);
    //             for (uint32_t i = 0; i < entry_size; ++i, ++entry) {
    //                 m_table.push_back(*entry);
    //             }
    //             ++m_size;

    //             // logging
    //             double cost_saving = compute_saving(entry_size, freq, m_total_integers);
    //             m_total_coverage += freq * entry_size * 100.0 / m_total_integers;
    //             m_final_bpi -= cost_saving;
    //             if (m_size % 1000 == 0) {
    //                 logger() << "entries in dictionary " << m_size << "/" << num_entries << std::endl;
    //                 logger() << "current bits x integer: " << m_final_bpi << std::endl;
    //                 logger() << "covering " << m_total_coverage << "% of integers" << std::endl;
    //             }

    //             return true;
    //         }

    //         void prepare_for_encoding() {
    //             std::vector<uint32_t> run(256, 0);
    //             uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
    //             uint32_t i = 1;
    //             for (uint32_t n = 256; n >= max_entry_size; n /= 2, ++i) {
    //                 uint64_t hash = hash_bytes64(byte_range(ptr, ptr + n * sizeof(uint32_t)));
    //                 m_map[hash] = i;
    //             }
    //             for (; i < num_entries; ++i) {
    //                 uint8_t const* ptr = reinterpret_cast<uint8_t const*>(get(i));
    //                 uint32_t entry_size = size(i);
    //                 uint64_t hash = hash_bytes64(byte_range(ptr, ptr + entry_size * sizeof(uint32_t)));
    //                 m_map[hash] = i;
    //             }
    //         }

    //         // Giulio: return the index of the pattern if found in the table,
    //         // otherwise return the invalid_index value
    //         uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
    //         {
    //             uint8_t const* b = reinterpret_cast<uint8_t const*>(begin);
    //             uint8_t const* e = b + entry_size * sizeof(uint32_t);
    //             uint64_t hash = hash_bytes64(byte_range(b, e));
    //             auto it = m_map.find(hash);
    //             if (it != m_map.end()) {
    //                 assert((*it).second <= num_entries);
    //                 return (*it).second;
    //             }
    //             return invalid_index;
    //         }

    //         void build(dictionary& dict) {
    //             std::swap(m_size, dict.m_size);
    //             dict.m_offsets.steal(m_offsets);
    //             dict.m_table.steal(m_table);
    //             builder().swap(*this);
    //         }

    //         void swap(builder& other) {
    //             std::swap(m_size, other.m_size);
    //             std::swap(m_final_bpi, other.m_final_bpi);
    //             std::swap(m_total_coverage, other.m_total_coverage);
    //             std::swap(m_total_integers, other.m_total_integers);
    //             m_offsets.swap(other.m_offsets);
    //             m_table.swap(other.m_table);
    //             m_map.swap(other.m_map);
    //         }

    //         uint32_t size() const {
    //             return m_size;
    //         }

    //         double bpi() const {
    //             return m_final_bpi;
    //         }

    //         double coverage() const {
    //             return m_total_coverage;
    //         }

    //         static std::string type() {
    //             return "packed";
    //         }

    //     private:
    //         std::string m_type;
    //         uint32_t m_size;
    //         double m_final_bpi;
    //         double m_total_coverage;
    //         uint64_t m_total_integers;
    //         std::vector<uint32_t> m_offsets;
    //         std::vector<uint32_t> m_table;

    //         // map from hash codes to table indexes, used during encoding
    //         std::unordered_map<uint64_t, uint32_t> m_map;

    //         uint32_t const* get(uint32_t i) const {
    //             assert(i < size());
    //             uint32_t offset = m_offsets[i];
    //             return &m_table[offset + 1]; // skip size
    //         }

    //         uint32_t size(uint32_t i) const {
    //             assert(i < size());
    //             return m_table[m_offsets[i]];
    //         }
    //     };

    //     dictionary()
    //         : m_size(0)
    //     {}

    //     uint32_t copy(uint32_t i, uint32_t* out) const
    //     {
    //         assert(i < size());
    //         uint32_t offset = m_offsets[i];
    //         uint32_t size = m_table[offset];
    //         uint32_t const* ptr = &m_table[offset + 1];
    //         memcpy(out, ptr, max_entry_size * sizeof(uint32_t));
    //         return size;
    //     }

    //     size_t size() const {
    //         return m_size;
    //     }

    //     void swap(dictionary& other) {
    //         std::swap(m_size, other.m_size);
    //         m_offsets.swap(other.m_offsets);
    //         m_table.swap(other.m_table);
    //     }

    //     template<typename Visitor>
    //     void map(Visitor& visit)
    //     {
    //         visit
    //             (m_size, "m_size")
    //             (m_offsets, "m_offsets")
    //             (m_table, "m_table")
    //             ;
    //     }

    // private:
    //     uint32_t m_size;
    //     succinct::mapper::mappable_vector<uint32_t> m_offsets;
    //     succinct::mapper::mappable_vector<uint32_t> m_table;
    // };


    // NOTE: as used by LSO
    // it's a compact dictionary where
    // all targets' offsets are multiple of max_entry_size.
    template<uint32_t t_num_entries,
             uint32_t t_max_entry_size>
    struct dictionary
    {
        // static_assert(is_power_of_two(t_max_entry_size));
        static const uint32_t num_entries = t_num_entries;
        static const uint32_t max_entry_size = t_max_entry_size;
        static const uint32_t invalid_index = uint32_t(-1);
        static const uint32_t reserved = 6; // 1 for exceptions
                                            // 5 for runs: 256, 128, 64, 32 and 16

        struct builder
        {
            static const uint32_t num_entries = dictionary::num_entries;
            static const uint32_t max_entry_size = dictionary::max_entry_size;
            static const uint32_t invalid_index = dictionary::invalid_index;
            static const uint32_t reserved = dictionary::reserved;

            uint64_t codewords = 0;
            uint64_t exceptions = 0;

            builder()
                : m_size(reserved)
                , m_final_bpi(0.0)
                , m_total_coverage(0.0)
                , m_total_integers(0)
            {}

            void init(uint64_t total_integers = 0) {
                m_size = reserved;
                // at the beginning, everything is exception(al)
                m_final_bpi = constants::initial_bpi;
                m_total_coverage = 0.0;
                m_total_integers = total_integers;
                m_table.reserve(num_entries);
                m_offsets.reserve(2 * num_entries);
                for (uint32_t i = 0; i < reserved; ++i) { // unused offsets
                    m_offsets.push_back(0);
                    m_offsets.push_back(0);
                }
            }

            void load_from_file(std::string dict_file) {
                std::ifstream ifs(dict_file);
                load(ifs);
            }

            bool try_store_to_file(std::string dict_file) const {
                std::ofstream ofs(dict_file);
                if (ofs) {
                    write(ofs);
                    return true;
                }
                return false;
            }

            void write(std::ofstream& dictionary_file) const {
                // logger() << "saving " << m_size << " entries" << std::endl;
                uint32_t offsets_size = m_offsets.size();
                uint32_t table_size = m_table.size();
                dictionary_file.write(reinterpret_cast<char const*>(&m_size), sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(&offsets_size), sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(&table_size), sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(m_offsets.data()), offsets_size * sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), table_size * sizeof(uint32_t));
            }

            void load(std::ifstream& dictionary_file) {
                uint32_t offsets_size = 0;
                uint32_t table_size = 0;
                dictionary_file.read(reinterpret_cast<char*>(&m_size), sizeof(uint32_t));
                // logger() << "loading " << m_size << " entries" << std::endl;
                dictionary_file.read(reinterpret_cast<char*>(&offsets_size), sizeof(uint32_t));
                dictionary_file.read(reinterpret_cast<char*>(&table_size), sizeof(uint32_t));
                // logger() << "offsets_size " << offsets_size << std::endl;
                // logger() << "table_size " << table_size << std::endl;
                m_table.resize(table_size);
                m_offsets.resize(offsets_size);
                dictionary_file.read(reinterpret_cast<char*>(m_offsets.data()), offsets_size * sizeof(uint32_t));
                dictionary_file.read(reinterpret_cast<char*>(m_table.data()), table_size * sizeof(uint32_t));

                // for (auto x: m_offsets) {
                //     std::cerr << x << " ";
                // }
                // std::cerr << std::endl;

                // for (uint32_t i = 2 * reserved, j = reserved; i < 2 * m_size; i += 2, ++j) {
                //     uint32_t s = m_offsets[i];
                //     std::cerr << "size: " << s << std::endl;
                //     if (s != size(j)) {
                //         std::cerr << "ERROR: got " << size(j) << ", but expected " << s << std::endl;
                //     }
                //     uint32_t offset = m_offsets[i + 1];
                //     std::cerr << "offset: " << offset << std::endl;
                //     auto ptr = get(j);
                //     for (uint32_t k = 0; k < s; ++k) {
                //         std::cerr << *ptr << " ";
                //         ++ptr;
                //     }
                //     std::cerr << std::endl;
                // }
            }

            bool full() {
                return m_size == num_entries;
            }

            bool append(uint32_t const* entry, uint32_t entry_size, uint64_t freq)
            {
                (void) freq;

                if (entry_size > max_entry_size) {
                    throw std::runtime_error("entry not allowed");
                }

                if (full()) return false;

                // 1. add the current entry
                uint32_t entry_offset = m_table.size();
                m_offsets.push_back(entry_size);
                m_offsets.push_back(entry_offset);
                m_set.insert(hash_bytes64(entry, entry_size));
                for (auto ptr = entry; ptr != entry + entry_size; ++ptr) {
                    m_table.push_back(*ptr);
                }
                assert(m_table.size() % max_entry_size == 0);

                ++m_size;
                if (full()) return false;

                // 2. loop through the prefixes and check if they appear already in the dictionary
                // by using the unordered_set; if not, just add a new entry with the same address as the current entry
                // and proper target size
                assert(constants::target_sizes[0] == max_entry_size);
                for (uint32_t s = 1; s < constants::max_fractal_steps; ++s) {
                    uint32_t target_size = constants::target_sizes[s];
                    uint64_t hash = hash_bytes64(entry, target_size);
                    if (m_set.find(hash) == m_set.cend()) {
                        m_offsets.push_back(target_size);
                        m_offsets.push_back(entry_offset);
                        m_set.insert(hash);
                        ++m_size;
                        if (full()) return false;
                    }
                }

                return true;
            }

            void prepare_for_encoding() {
                std::vector<uint32_t> run(256, 0);
                uint32_t i = 1;
                for (uint32_t n = 256; n >= max_entry_size; n /= 2, ++i) {
                    // logger() << "run of size " << n << " with codeword " << i << std::endl;
                    uint64_t hash = hash_bytes64(run.data(), n);
                    m_map[hash] = i;
                }
                for (; i < num_entries; ++i) {
                    uint32_t entry_size = size(i);
                    // logger() << "entry of size " << entry_size << ": ";
                    // auto ptr = get(i);
                    // for (uint32_t k = 0; k < entry_size; ++k) {
                    //     std::cerr << *ptr << " ";
                    //     ++ptr;
                    // }
                    // std::cerr << std::endl;
                    uint64_t hash = hash_bytes64(get(i), entry_size);
                    m_map[hash] = i;
                }
            }

            // Giulio: return the index of the pattern if found in the table,
            // otherwise return the invalid_index value
            uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
            {
                uint64_t hash = hash_bytes64(begin, entry_size);
                auto it = m_map.find(hash);
                if (it != m_map.end()) {
                    assert((*it).second <= num_entries);
                    return (*it).second;
                }
                return invalid_index;
            }

            void build(dictionary& dict) {
                std::swap(m_size, dict.m_size);
                dict.m_offsets.steal(m_offsets);
                dict.m_table.steal(m_table);
                builder().swap(*this);
            }

            void swap(builder& other) {
                std::swap(m_size, other.m_size);
                std::swap(m_final_bpi, other.m_final_bpi);
                std::swap(m_total_coverage, other.m_total_coverage);
                std::swap(m_total_integers, other.m_total_integers);
                m_offsets.swap(other.m_offsets);
                m_table.swap(other.m_table);
                m_map.swap(other.m_map);
            }

            uint32_t size() const {
                return m_size;
            }

            double bpi() const {
                return m_final_bpi;
            }

            double coverage() const {
                return m_total_coverage;
            }

            static std::string type() {
                return "packed-LSO";
            }

        private:
            std::string m_type;
            uint32_t m_size;
            double m_final_bpi;
            double m_total_coverage;
            uint64_t m_total_integers;

            std::vector<uint32_t> m_offsets;
            std::vector<uint32_t> m_table;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;

            // set to check if a string is already in the dictionary or not
            std::unordered_set<uint64_t> m_set;

            uint32_t const* get(uint32_t i) const {
                assert(i < size());
                return &m_table[m_offsets[i * 2 + 1]];
            }

            uint32_t size(uint32_t i) const {
                assert(i < size());
                return m_offsets[i * 2];
            }
        };

        dictionary()
            : m_size(0)
        {}

        uint32_t copy(uint32_t i, uint32_t* out) const
        {
            assert(i < size());
            uint32_t size = m_offsets[i * 2];
            uint32_t offset = m_offsets[i * 2 + 1];
            uint32_t const* ptr = &m_table[offset];

            // NOTE: if we copy max_entry_size, then we cannot
            // skip anymore because we will not have 0s: must be overwritten
            // this is also true for the packed variant as well

            memcpy(out, ptr, max_entry_size * sizeof(uint32_t)); // NOTE: this will decode uncorrectly, since we will copy some garbage
            // memcpy(out, ptr, size * sizeof(uint32_t));

            return size;
        }

        size_t size() const {
            return m_size;
        }

        void swap(dictionary& other) {
            std::swap(m_size, other.m_size);
            m_offsets.swap(other.m_offsets);
            m_table.swap(other.m_table);
        }

        template<typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_size, "m_size")
                (m_offsets, "m_offsets")
                (m_table, "m_table")
                ;
        }

    private:
        uint32_t m_size;
        succinct::mapper::mappable_vector<uint32_t> m_offsets;
        succinct::mapper::mappable_vector<uint32_t> m_table;
    };

}
