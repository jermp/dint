#pragma once

#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <limits>

#include <succinct/mappable_vector.hpp>

#include "dint_configuration.hpp"
#include "model_build_utils.hpp"
#include "hash_utils.hpp"
#include "util.hpp"

// #include <x86intrin.h>
// #include <immintrin.h> // AVX2

#define EXCEPTIONS 2
#define INF std::numeric_limits<uint32_t>::max()

namespace ds2i {

    // RECTANGULAR
    template<uint32_t t_num_entries,
             uint32_t t_max_entry_size>
    struct dictionary
    {
        static_assert(is_power_of_two(t_max_entry_size));
        static const uint32_t num_entries = t_num_entries;
        static const uint32_t max_entry_size = t_max_entry_size;
        static const uint32_t invalid_index = uint32_t(-1);
        static const uint32_t reserved = EXCEPTIONS + 5;

        struct builder
        {
            static const uint32_t num_entries = dictionary::num_entries;
            static const uint32_t max_entry_size = dictionary::max_entry_size;
            static const uint32_t invalid_index = dictionary::invalid_index;
            static const uint32_t reserved = dictionary::reserved;

            uint64_t codewords = 0;
            uint64_t small_exceptions = 0;
            uint64_t large_exceptions = 0;

            uint64_t block_encoded_with_large_dict = 0;
            uint64_t block_encoded_with_small_dict = 0;

            builder()
                : m_pos(0)
                , m_size(reserved)
                , m_table(0, 0)
            {}

            size_t num_bytes_for(uint64_t n) const {
                return n * (max_entry_size + 1) * sizeof(m_table.front());
            }

            void init() {
                m_pos = reserved * (max_entry_size + 1);
                m_size = reserved;
                m_table.resize(num_entries * (max_entry_size + 1), 0);

                uint32_t pos = max_entry_size + 1;
                for (int i = 0; i < EXCEPTIONS; ++i, pos += max_entry_size + 1) {
                    m_table[pos - 1] = 1;
                }
                for (int i = 0, size = 256; i < 5; ++i, pos += max_entry_size + 1, size /= 2) {
                    m_table[pos - 1] = size;
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
                dictionary_file.write(reinterpret_cast<char const*>(&m_size), sizeof(uint32_t));
                dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), num_bytes_for(m_size));
            }

            void load(std::ifstream& dictionary_file) {
                uint32_t size = 0;
                dictionary_file.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
                init();
                m_size = size;
                dictionary_file.read(reinterpret_cast<char*>(m_table.data()), num_bytes_for(m_size));
            }

            bool full() {
                return m_size == num_entries;
            }

            bool append(uint32_t const* entry, uint32_t entry_size)
            {
                if (full()) return false;
                assert(m_pos % (max_entry_size + 1) == 0);
                std::copy(entry, entry + entry_size, &m_table[m_pos]);
                m_pos += max_entry_size + 1;
                m_table[m_pos - 1] = entry_size;
                ++m_size;
                return true;
            }

            void prepare_for_encoding()
            {
                std::vector<uint32_t> run(256, 0);
                uint32_t i = EXCEPTIONS;
                for (uint32_t n = 256; n >= 16; n /= 2, ++i) {
                    uint64_t hash = hash_bytes64(run.data(), n);
                    m_map[hash] = i;
                }
                for (; i < size(); ++i) {
                    // auto ptr = get(i);
                    // std::cout << i << ": " << std::endl;
                    // for (uint32_t k = 0; k < size(i); ++k) {
                    //     std::cout << *ptr << " ";
                    //     ++ptr;
                    // }
                    // std::cout << std::endl;
                    uint64_t hash = hash_bytes64(get(i), size(i));
                    m_map[hash] = i;
                }
            }

            // NOTE: return the index of the pattern if found in the table,
            // otherwise return the invalid_index value
            uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
            {
                uint64_t hash = hash_bytes64(begin, entry_size);
                auto it = m_map.find(hash);
                if (it != m_map.end()) {
                    assert((*it).second < num_entries);
                    return (*it).second;
                }
                return invalid_index;
            }

            void build(dictionary& dict) {
                dict.m_table.steal(m_table);
                builder().swap(*this);
            }

            void swap(builder& other) {
                std::swap(m_pos, other.m_pos);
                std::swap(m_size, other.m_size);
                m_table.swap(other.m_table);
                m_map.swap(other.m_map);
            }

            uint32_t size() const {
                return m_size;
            }

            static std::string type() {
                return "rectangular";
            }

            // print vocabulary entries usage
            void print_usage() {
                std::vector<uint32_t> sizes;
                sizes.push_back(EXCEPTIONS);
                for (uint32_t i = 0; i < constants::num_target_sizes; ++i) {
                    sizes.push_back(0);
                }
                sizes.push_back(5); // for the runs

                for (uint64_t i = 0, k = 0; i < m_table.size(); i += max_entry_size + 1, ++k) {
                    uint32_t size = m_table[i + max_entry_size];
                    // std::cout << k << ": " << size << " - [";
                    // for (uint32_t k = 0; k < max_entry_size; ++k) {
                    //     std::cout << m_table[i + k];
                    //     if (k != max_entry_size - 1) std::cout << "|";
                    // }
                    // std::cout << "]" << std::endl;
                    sizes[ceil_log2(size) + 1] += 1;
                }

                std::cout << "rare: " << EXCEPTIONS << " (" << EXCEPTIONS * 100.0 / num_entries << "%)" << std::endl;
                for (uint32_t i = 0; i < constants::num_target_sizes; ++i) {
                    std::cout << "entries of size "
                              << (uint32_t(1) << i) << ": "
                              << sizes[i + 1]
                              << "(" << sizes[i + 1] * 100.0 / num_entries << "%)"
                              << std::endl;
                }
                std::cout << "freq.: 5 (" << 5 * 100.0 / num_entries << "%)" << std::endl;
            }

            uint32_t size(uint32_t i) const {
                uint32_t begin = i * (t_max_entry_size + 1);
                uint32_t end = begin + t_max_entry_size;
                return m_table[end];
            }

            uint32_t const* get(uint32_t i) const {
                uint32_t begin = i * (t_max_entry_size + 1);
                return &m_table[begin];
            }

        private:
            uint32_t m_pos;
            uint32_t m_size;
            std::vector<uint32_t> m_table;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;
        };

        dictionary()
        {}

        uint32_t copy(uint32_t i, uint32_t* out) const
        {
            assert(i < num_entries);
            uint32_t begin = i * (max_entry_size + 1);
            uint32_t const* ptr = &m_table[begin];

            // APPROACH 1: always copy max_entry_size * sizeof(uint32_t) bytes,
            // regardless the fact that most entries need less bytes
            // std::cout << "i " << i << "/" << num_entries << std::endl;
            memcpy(out, ptr, max_entry_size * sizeof(uint32_t));
            // std::cout << "done" << std::endl;

            // NOTE: slower than memcpy
            // *(out + 0) = *(ptr + 0);
            // *(out + 1) = *(ptr + 1);
            // *(out + 2) = *(ptr + 2);
            // *(out + 3) = *(ptr + 3);
            // *(out + 4) = *(ptr + 4);
            // *(out + 5) = *(ptr + 5);
            // *(out + 6) = *(ptr + 6);
            // *(out + 7) = *(ptr + 7);

            // NOTE: faster than previous one, but still a bit slower than memcpy
            // uint64_t* _out = reinterpret_cast<uint64_t*>(out);
            // uint64_t const* _in = reinterpret_cast<uint64_t const*>(ptr);
            // *(_out + 0) = *(_in + 0);
            // *(_out + 1) = *(_in + 1);

            // APPROACH 2: split the copy into 2 parts and copy
            // the rest only when necessary
            // NOTE: slower than memcpy
            // memcpy(out, ptr, 16);
            // if (size == 8) {
            //     memcpy(out + 4, ptr + 4, 16);
            // }
            // *(reinterpret_cast<uint64_t*>(out)) = *(reinterpret_cast<uint64_t const*>(ptr));
            // if (size == 4) {
            //     *(reinterpret_cast<uint64_t*>(out) + 1) = *(reinterpret_cast<uint64_t const*>(ptr) + 1);
            // }

            // APPROACH 3: SIMD
            // NOTE: equal to memcpy
            // _mm_storeu_si128((__m128i*) out, _mm_loadu_si128((__m128i const*) ptr)); // l = 4
            // _mm256_storeu_si256((__m256i*) out, _mm256_lddqu_si256((__m256i const*) ptr)); // l = 8
            // _mm512_storeu_si512(out, _mm512_loadu_si512(ptr)); // l = 16

            // __m128i tmp_0 = _mm_loadu_si128((__m128i*) ptr);
            // _mm_storeu_si128((__m128i*) out, tmp_0);
            // if (size == 8) {
            //     __m128i tmp_4 = _mm_loadu_si128((__m128i*) (ptr + 4));
            //     _mm_storeu_si128((__m128i*)(out + 4), tmp_4);
            // }

            // uint32_t end = begin + max_entry_size;
            uint32_t size = *(ptr + max_entry_size); // m_table[end];

            // // APPROACH 4: copy exactly the needed bytes
            // memcpy(out, ptr, size * sizeof(uint32_t));

            return size;
        }

        // uint32_t size(uint32_t i) const {
        //     uint32_t begin = i * (t_max_entry_size + 1);
        //     uint32_t end = begin + t_max_entry_size;
        //     return m_table[end];
        // }

        void swap(dictionary& other) {
            m_table.swap(other.m_table);
        }

        template<typename Visitor>
        void map(Visitor& visit) {
            visit(m_table, "m_table");
        }

        // void print(uint32_t i) {
        //     uint32_t begin = i * (max_entry_size + 1);
        //     uint32_t end = begin + max_entry_size;
        //     uint32_t size = m_table[end];
        //     uint32_t const* ptr = &m_table[begin];
        //     std::cout << "[";
        //     for (uint32_t k = 0; k < size; ++k) {
        //         std::cout << *ptr;
        //         ++ptr;
        //         if (k != size - 1) std::cout << "|";
        //     }
        //     std::cout << "]" << std::endl;
        // }

    private:
        succinct::mapper::mappable_vector<uint32_t> m_table;
    };

    using large_dictionary_type = dictionary
                                <constants::num_entries,
                                 constants::max_entry_size>;

    using small_dictionary_type = dictionary
                                <256,
                                 large_dictionary_type::max_entry_size>;



    // PACKED
    template<uint32_t t_num_entries,
             uint32_t t_max_entry_size>
    struct dictionary_packed
    {
        static_assert(is_power_of_two(t_max_entry_size));
        static const uint32_t num_entries = t_num_entries;
        static const uint32_t max_entry_size = t_max_entry_size;
        static const uint32_t invalid_index = uint32_t(-1);
        static const uint32_t reserved = EXCEPTIONS + 5;

        struct builder
        {
            static const uint32_t num_entries = dictionary_packed::num_entries;
            static const uint32_t max_entry_size = dictionary_packed::max_entry_size;
            static const uint32_t invalid_index = dictionary_packed::invalid_index;
            static const uint32_t reserved = dictionary_packed::reserved;

            uint64_t codewords = 0;
            uint64_t small_exceptions = 0;
            uint64_t large_exceptions = 0;

            uint64_t block_encoded_with_large_dict = 0;
            uint64_t block_encoded_with_small_dict = 0;

            builder()
                : m_size(reserved)
            {}

            void init() {
                m_size = reserved;
                m_offsets.reserve(num_entries);

                // NOTE: push [max_entry_size] 0s at the beginning of the table
                // to be copied in case of a run, thus the offset corresponding to
                // the indexes 2, 3, 4, 5, 6 (i.e., runs) is always 0
                for (uint32_t i = 0; i != max_entry_size; ++i) {
                    m_table.push_back(0);
                }
                for (uint32_t i = 0; i != EXCEPTIONS; ++i) {
                    m_offsets.push_back(0);
                }
                for (uint32_t i = 0, size = 256; i != 5; ++i, size /= 2) {
                    uint32_t size_and_offset = (size - 1) << 24; // offset is 0
                    m_offsets.push_back(size_and_offset);
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
                dictionary_file.read(reinterpret_cast<char*>(&offsets_size), sizeof(uint32_t));
                dictionary_file.read(reinterpret_cast<char*>(&table_size), sizeof(uint32_t));
                m_table.resize(table_size);
                m_offsets.resize(offsets_size);
                dictionary_file.read(reinterpret_cast<char*>(m_offsets.data()), offsets_size * sizeof(uint32_t));
                dictionary_file.read(reinterpret_cast<char*>(m_table.data()), table_size * sizeof(uint32_t));
            }

            bool full() {
                return m_size == num_entries;
            }

            bool append(uint32_t const* entry, uint32_t entry_size)
            {
                assert(entry_size > 0);
                if (entry_size > max_entry_size) {
                    throw std::runtime_error("entry not allowed");
                }
                if (full()) {
                    return false;
                }

                // Option 1: always push a new entry
                // uint32_t size_and_offset = ((entry_size - 1) << 24) | uint32_t(m_table.size());
                // m_offsets.push_back(size_and_offset);
                // for (uint32_t i = 0; i != entry_size; ++i, ++entry) {
                //     m_table.push_back(*entry);
                // }

                // Option 2: search for a prefix
                uint32_t S = (entry_size - 1) << 24;
                uint32_t O = 0;

                bool run = true;
                for (uint32_t k = 0; k != entry_size; ++k) {
                    if (entry[k] != 0) {
                        run = false;
                        break;
                    }
                }

                if (run) {

                    // std::cout << "searching for entry:" << std::endl;
                    // for (uint32_t k = 0; k != entry_size; ++k) {
                    //     std::cout << entry[k] << " ";
                    // }
                    // std::cout << std::endl;
                    // std::cout << "RUN!" << std::endl;

                    uint32_t size_and_offset = S;
                    m_offsets.push_back(size_and_offset);

                } else {
                    bool found = false;
                    O = max_entry_size;
                    for (uint32_t i = reserved; i != m_size; ++i) {
                        auto* ptr = get(i);
                        uint32_t s = size(i);
                        O = offset(i);

                        if (s > entry_size) {
                            found = true;
                            for (uint32_t k = 0; k != entry_size; ++k) {
                                if (ptr[k] != entry[k]) {
                                    found = false;
                                    break;
                                }
                            }

                            if (found) {

                                // std::cout << "searching for entry:" << std::endl;
                                // for (uint32_t k = 0; k != entry_size; ++k) {
                                //     std::cout << entry[k] << " ";
                                // }
                                // std::cout << std::endl;

                                // std::cout << "s " << s << std::endl;
                                // std::cout << "\tFOUND as prefix of:\n\t";
                                // for (uint32_t k = 0; k != s; ++k) {
                                //     std::cout << m_table[O + k] << " ";
                                // }
                                // std::cout << std::endl;
                                break;
                            }
                        }
                    }

                    if (found) {
                        uint32_t size_and_offset = S | O;
                        m_offsets.push_back(size_and_offset);
                    } else {
                        O = m_table.size();
                        uint32_t size_and_offset = S | O;
                        m_offsets.push_back(size_and_offset);
                        for (uint32_t i = 0; i != entry_size; ++i) {
                            m_table.push_back(entry[i]);
                        }
                    }

                }

                ++m_size;

                return true;
            }

            void prepare_for_encoding()
            {
                std::vector<uint32_t> run(256, 0);
                uint32_t i = EXCEPTIONS;
                for (uint32_t n = 256; n >= 16; n /= 2, ++i) {
                    uint64_t hash = hash_bytes64(run.data(), n);
                    m_map[hash] = i;
                }
                for (; i < size(); ++i) {
                    // auto ptr = get(i);
                    // std::cout << i << ": " << std::endl;
                    // for (uint32_t k = 0; k < size(i); ++k) {
                    //     std::cout << *ptr << " ";
                    //     ++ptr;
                    // }
                    // std::cout << std::endl;
                    uint64_t hash = hash_bytes64(get(i), size(i));
                    m_map[hash] = i;
                }
            }

            // NOTE: return the index of the pattern if found in the table,
            // otherwise return the invalid_index value
            uint32_t lookup(uint32_t const* begin, uint32_t entry_size) const
            {
                uint64_t hash = hash_bytes64(begin, entry_size);
                auto it = m_map.find(hash);
                if (it != m_map.end()) {
                    assert((*it).second < num_entries);
                    return (*it).second;
                }
                return invalid_index;
            }

            void build(dictionary_packed& dict) {
                dict.m_offsets.steal(m_offsets);
                dict.m_table.steal(m_table);
                builder().swap(*this);
            }

            void swap(builder& other) {
                std::swap(m_size, other.m_size);
                m_offsets.swap(other.m_offsets);
                m_table.swap(other.m_table);
                m_map.swap(other.m_map);
            }

            uint32_t size() const {
                return m_size;
            }

            static std::string type() {
                return "packed";
            }

            // print vocabulary entries usage
            void print_usage() {
                // TODO
                // std::vector<uint32_t> sizes;
                // sizes.push_back(EXCEPTIONS);
                // for (uint32_t i = 0; i < constants::num_target_sizes; ++i) {
                //     sizes.push_back(0);
                // }
                // sizes.push_back(5); // for the runs

                // for (uint64_t i = 0, k = 0; i < m_table.size(); i += max_entry_size + 1, ++k) {
                //     uint32_t size = m_table[i + max_entry_size];
                //     // std::cout << k << ": " << size << " - [";
                //     // for (uint32_t k = 0; k < max_entry_size; ++k) {
                //     //     std::cout << m_table[i + k];
                //     //     if (k != max_entry_size - 1) std::cout << "|";
                //     // }
                //     // std::cout << "]" << std::endl;
                //     sizes[ceil_log2(size) + 1] += 1;
                // }

                // std::cout << "rare: " << EXCEPTIONS << " (" << EXCEPTIONS * 100.0 / num_entries << "%)" << std::endl;
                // for (uint32_t i = 0; i < constants::num_target_sizes; ++i) {
                //     std::cout << "entries of size "
                //               << (uint32_t(1) << i) << ": "
                //               << sizes[i + 1]
                //               << "(" << sizes[i + 1] * 100.0 / num_entries << "%)"
                //               << std::endl;
                // }
                // std::cout << "freq.: 5 (" << 5 * 100.0 / num_entries << "%)" << std::endl;
            }

            uint32_t size(uint32_t i) const {
                assert(i < size());
                return (m_offsets[i] >> 24) + 1;
            }

            uint32_t offset(uint32_t i) const {
                assert(i < size());
                return m_offsets[i] & 0xFFFFFF;
            }

            uint32_t const* get(uint32_t i) const {
                return &m_table[offset(i)];
            }

        private:
            uint32_t m_size;
            std::vector<uint32_t> m_offsets;
            std::vector<uint32_t> m_table;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;
        };

        dictionary_packed()
        {}

        uint32_t copy(uint32_t i, uint32_t* out) const
        {
            assert(i < num_entries);
            uint32_t size_and_offset = m_offsets[i];
            uint32_t offset = size_and_offset & 0xFFFFFF;
            uint32_t size = (size_and_offset >> 24) + 1;
            uint32_t const* ptr = &m_table[offset];
            memcpy(out, ptr, max_entry_size * sizeof(uint32_t));
            return size;
        }

        void swap(dictionary_packed& other) {
            m_offsets.swap(other.m_offsets);
            m_table.swap(other.m_table);
        }

        template<typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_offsets, "m_offsets")
                (m_table, "m_table")
                ;
        }

    private:
        succinct::mapper::mappable_vector<uint32_t> m_offsets;
        succinct::mapper::mappable_vector<uint32_t> m_table;
    };


    // using large_dictionary_type = dictionary_packed
    //                             <constants::num_entries,
    //                              constants::max_entry_size>;

    // using small_dictionary_type = dictionary_packed
    //                             <256,
    //                              large_dictionary_type::max_entry_size>;



}
