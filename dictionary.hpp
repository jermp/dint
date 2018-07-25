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

#include <x86intrin.h>
#include <immintrin.h> // AVX2

#define EXCEPTIONS 2
#define INF std::numeric_limits<uint32_t>::max()

namespace ds2i {

    template<uint32_t t_num_entries,
             uint32_t t_max_entry_size>
    struct dictionary_rect
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

            uint64_t block_encoded_with_large_dict = 0;
            uint64_t block_encoded_with_small_dict = 0;

            builder()
                : m_pos(0)
                , m_size(reserved)
                , m_table(0, 0)
            {
            }

            void init() {
                m_pos = reserved * (max_entry_size + 1);
                m_size = reserved;
                m_table.resize(num_entries * (max_entry_size + 1), 0);
                // m_freqs.reserve(num_entries);

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
                dictionary_file.write(reinterpret_cast<char const*>(m_table.data()), m_table.size() * sizeof(uint32_t));
            }

            void load(std::ifstream& dictionary_file) {
                uint32_t size = 0;
                dictionary_file.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
                init();
                dictionary_file.read(reinterpret_cast<char*>(m_table.data()), m_table.size() * sizeof(uint32_t));
                m_size = size;
            }

            bool full() {
                return m_size == num_entries;
            }

            bool append(uint32_t const* entry, uint32_t entry_size, uint64_t /*freq*/)
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
                    uint64_t hash = hash_bytes64(get(i), size(i));
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
            void print() {
                std::vector<uint32_t> sizes;
                sizes.push_back(EXCEPTIONS);
                for (uint32_t i = 0; i < constants::num_target_sizes; ++i) {
                    sizes.push_back(0);
                }
                sizes.push_back(5); // for the runs

                for (uint64_t i = 0, k = 0; i < m_table.size(); i += max_entry_size + 1, ++k) {
                    uint32_t size = m_table[i + max_entry_size];
                    std::cout << k << ": " << size << " - [";
                    for (uint32_t k = 0; k < max_entry_size; ++k) {
                        std::cout << m_table[i + k];
                        if (k != max_entry_size - 1) std::cout << "|";
                    }
                    std::cout << "]" << std::endl;
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

        dictionary_rect()
        {}

        int const* data() const {
            return reinterpret_cast<int const*>(m_table.data());
        }

        uint32_t copy(uint32_t i, uint32_t* out) const
        {
            assert(i < num_entries);
            uint32_t begin = i * (max_entry_size + 1);
            uint32_t const* ptr = &m_table[begin];

            memcpy(out, ptr, max_entry_size * sizeof(uint32_t));

            uint32_t size = *(ptr + max_entry_size); // m_table[end];

            return size;
        }

        uint32_t size(uint32_t i) const {
            uint32_t begin = i * (t_max_entry_size + 1);
            uint32_t end = begin + t_max_entry_size;
            return m_table[end];
        }

        void swap(dictionary_rect& other) {
            m_table.swap(other.m_table);
        }

        template<typename Visitor>
        void map(Visitor& visit) {
            visit(m_table, "m_table");
        }

    private:
        succinct::mapper::mappable_vector<uint32_t> m_table;
    };

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
            static const uint32_t num_entries = dictionary::num_entries;
            static const uint32_t max_entry_size = dictionary::max_entry_size;
            static const uint32_t invalid_index = dictionary::invalid_index;
            static const uint32_t reserved = dictionary::reserved;

            builder()
                : m_size(reserved)
            {}

            void init(uint64_t total_integers = 0) {
                m_size = reserved;
                m_offsets.reserve(num_entries);
                for (uint32_t i = 0; i < reserved; ++i) { // unused offsets
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

            bool append(uint32_t const* entry, uint32_t entry_size, uint64_t freq)
            {
                if (entry_size > max_entry_size) {
                    throw std::runtime_error("entry not allowed");
                }
                if (full()) {
                    return false;
                }

                uint32_t size_and_offset = (entry_size << 16) + uint32_t(m_table.size());
                m_offsets.push_back(size_and_offset);
                for (uint32_t i = 0; i < entry_size; ++i, ++entry) {
                    m_table.push_back(*entry);
                }
                ++m_size;
                return true;
            }

            void prepare_for_encoding() {
                std::vector<uint32_t> run(256, 0);
                uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
                uint32_t i = 1;
                for (uint32_t n = 256; n >= max_entry_size; n /= 2, ++i) {
                    uint64_t hash = hash_bytes64(byte_range(ptr, ptr + n * sizeof(uint32_t)));
                    m_map[hash] = i;
                }
                for (; i < num_entries; ++i) {
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

        private:
            uint32_t m_size;
            std::vector<uint32_t> m_offsets;
            std::vector<uint32_t> m_table;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;

            uint32_t const* get(uint32_t i) const {
                assert(i < size());
                uint32_t offset = m_offsets[i] & 0xFFFF;
                return &m_table[offset]; // skip size
            }

            uint32_t size(uint32_t i) const {
                assert(i < size());
                return m_offsets[i] >> 16;
            }
        };

        dictionary_packed()
            : m_size(0)
        {}

        uint32_t copy(uint32_t i, uint32_t* out) const
        {
            assert(i < size());
            uint32_t offset = m_offsets[i] & 0xFFFF;
            uint32_t size = m_offsets[i] >> 16;
            uint32_t const* ptr = &m_table[offset];
            memcpy(out, ptr, max_entry_size * sizeof(uint32_t));
            return size;
        }

        size_t size() const {
            return m_size;
        }

        void swap(dictionary_packed& other) {
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


    template<uint32_t t_num_entries,
             uint32_t t_max_entry_size>
    struct dictionary_overlap
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

            builder()
                : m_size(reserved)
            {}

            void init(uint64_t total_integers = 0) {
                m_size = reserved;
                m_offsets.reserve(num_entries);
                for (uint32_t i = 0; i < reserved; ++i) { // unused offsets
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

            bool append(uint32_t const* entry, uint32_t entry_size, uint64_t freq)
            {
                if (entry_size > max_entry_size) {
                    throw std::runtime_error("entry not allowed");
                }
                if (full()) {
                    return false;
                }

                uint32_t size_and_offset = (entry_size << 16) + uint32_t(m_table.size());
                m_offsets.push_back(size_and_offset);
                for (uint32_t i = 0; i < entry_size; ++i, ++entry) {
                    m_table.push_back(*entry);
                }
                ++m_size;
                return true;
            }

            struct target_node_t {
                std::vector<uint32_t> sizes;
                std::vector<uint32_t> entry;
                std::vector<target_node_t> contained_targets;
                bool operator<(const target_node_t& other) const {
                    return sizes[0] < other.sizes[0];
                }
            };

            template <typename Container>
            bool is_contained(const Container& cont, const std::Container& elem)
            {
                if(elem.size() >= cont.size()) return false;
                return std::search(cont.begin(), cont.end(), elem.begin(), elem.end()) != cont.end();
            }

            void find_substring_overlaps(std::vector<target_node_t>& entries) {
                std::sort(entries.begin(),entries.end());
                for(size_t i=0;i<entries.size();i++) {
                    auto& cur = entries[i];
                    int max_contained = -1;
                    for(size_t j=i;j<entries.size();j++) {
                        if(is_contained(entries[j],cur)) max_contained = j;
                    }
                    if(max_contained != -1) {
                        auto& dest = entries[max_contained];
                        for(auto& ct : cur.contained_targets) {
                            dest.contained_targets.push_back(ct);
                        }
                        cur.contained_targets.clear();
                        dest.contained_targets.push_back(cur);
                        cur.sizes[0] = -1; // mark as invalid
                    }
                }
                auto itr = entries.begin();
                auto end = entries.end();
                while(itr != end) {
                    if(itr->sizes[0] == -1) {
                        itr = entries.erase(itr);
                    } else {
                        ++itr;
                    }
                }
            }

            size_t compute_overlap(const std::vector<uint32_t>& A,const std::vector<uint32_t>& B) const {
                auto res = std::mismatch(A.rbegin(),A.rend(),B.begin(),B.end());
                return std::distance(B.begin(),res.second);
            }

            std::tuple<size_t,size_t,size_t> find_largest_overlap(std::vector<target_node_t>& entries) const {
                size_t max_overlap = 0;
                size_t max_index_left = 0;
                size_t max_index_right = 0;
                for(size_t i=0;i<entries.size();i++) {
                    for(size_t j=0;j<entries.size();j++) {
                        if(i != j) {
                            auto overlap = compute_overlap(entries[i].entry,entries[j].entry);
                            if(overlap > max_overlap) {
                                max_overlap = overlap;
                                max_index_left = i;
                                max_index_right = j;
                            }
                        }
                    }
                }
                return make_tuple(max_overlap,max_index_left,max_index_right);
            }

            void perform_greedy_prefixsuffix_overlap(std::vector<target_node_t>& entries) {
                size_t overlap,id_left,id_right;
                std::tie(overlap,id_left,id_right) = find_largest_overlap(entries);
                while(overlap != 0) {
                    // (a) merge the nodes
                    auto& node_left = entries[id_left];
                    auto copy_itr = entries[id_right].entry.begin() + overlap;
                    std::copy(copy_itr,entries[id_right].entry.end(),std::back_inserter(node_left.entry));
                    std::copy(entries[id_right].sizes.begin(),entries[id_right].sizes.end(),std::back_inserter(node_left.sizes));
                    // (b) remove the right node
                    entries.erase(entries.begin()+id_right);
                    // (c) search again
                    std::tie(overlap,id_left,id_right) = find_largest_overlap(entries);
                }

            }

            void find_overlaps_and_compact() {
                // (1) create initial 'nodes'
                std::vector<target_node_t> entries;
                for(size_t i=0;i<m_size;i++) {
                    target_node_t node;
                    auto start = m_offsets[i] & 0xFFFF;
                    auto size = m_offsets[i] >> 16;
                    auto itr = m_table.begin() + start + 1;
                    std::copy(itr,itr + size,std::back_inserter(node.entry));
                    node.sizes.push_back(size);
                    entries.push_back(node);
                }

                // (2) try to find overlaps inside targets
                find_substring_overlaps(entries);

                // (3) perform the greedy prefix-suffix overlap 
                perform_greedy_prefixsuffix_overlap(entries);

                // (4) create compacted dictionary string and offset array
                m_table.clear();
                m_offsets.clear();
                for(auto& cur : entries) {
                    auto itr = cur.entry.begin();
                    auto size_before = m_table.size();
                    for(auto size : cur.sizes) {
                        uint32_t size_and_offset = (size << 16) + uint32_t(m_table.size());
                        m_offsets.push_back(size_and_offset);
                        std::copy(itr,itr + size,std::back_inserter(m_table));
                        itr += size;
                    }
                    auto entry_start = m_table.begin() + size_before;
                    for(auto& ct : cur.contained_targets) {
                        auto itr = std::search(entry_start,m_table.end(),ct.entry.begin(),ct.entry.end());
                        uint32_t offset = std::distance(m_table.begin(),itr);
                        uint32_t size_and_offset = (ct.sizes[0] << 16) + uint32_t(offset);
                        m_offsets.push_back(size_and_offset);
                    }
                }
            }

            void prepare_for_encoding() {
                find_overlaps_and_compact();

                std::vector<uint32_t> run(256, 0);
                uint8_t const* ptr = reinterpret_cast<uint8_t const*>(run.data());
                uint32_t i = 1;
                for (uint32_t n = 256; n >= max_entry_size; n /= 2, ++i) {
                    uint64_t hash = hash_bytes64(byte_range(ptr, ptr + n * sizeof(uint32_t)));
                    m_map[hash] = i;
                }
                for (; i < num_entries; ++i) {
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
                    assert((*it).second <= num_entries);
                    return (*it).second;
                }
                return invalid_index;
            }

            void build(dictionary_overlap& dict) {
                std::swap(m_size, dict.m_size);
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
                return "overlap";
            }

        private:
            uint32_t m_size;
            std::vector<uint32_t> m_offsets;
            std::vector<uint32_t> m_table;

            // map from hash codes to table indexes, used during encoding
            std::unordered_map<uint64_t, uint32_t> m_map;

            uint32_t const* get(uint32_t i) const {
                assert(i < size());
                uint32_t offset = m_offsets[i] & 0xFFFF;
                return &m_table[offset]; // skip size
            }

            uint32_t size(uint32_t i) const {
                assert(i < size());
                return m_offsets[i] >> 16;
            }
        };

        dictionary_overlap()
            : m_size(0)
        {}

        uint32_t copy(uint32_t i, uint32_t* out) const
        {
            assert(i < size());
            uint32_t offset = m_offsets[i] & 0xFFFF;
            uint32_t size = m_offsets[i] >> 16;
            uint32_t const* ptr = &m_table[offset];
            memcpy(out, ptr, max_entry_size * sizeof(uint32_t));
            return size;
        }

        size_t size() const {
            return m_size;
        }

        void swap(dictionary_packed& other) {
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
