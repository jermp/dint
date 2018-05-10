#pragma once

#include <succinct/mappable_vector.hpp>
#include <succinct/bit_vector.hpp>

#include "compact_elias_fano.hpp"
#include "dict_posting_list.hpp"

namespace ds2i {

    template<uint32_t t_dict_entries = 65536,uint32_t t_dict_entry_width = 16>
    struct dict_builder_alistair {
        static const uint32_t dict_entries = t_dict_entries;
        static const uint32_t dict_entry_width = t_dict_entry_width;

        static void build(std::vector<uint32_t>& dict,std::vector<uint32_t>& block_stats)
        {
            const size_t dict_size = dict_entries * dict_entry_width;
            dict.resize(dict_size);
        }
    };

    template<uint32_t t_dict_entries = 65536,uint32_t t_dict_entry_width = 16>
    struct dict_builder_giulio {
        static const uint32_t dict_entries = t_dict_entries;
        static const uint32_t dict_entry_width = t_dict_entry_width;

        static void build(std::vector<uint32_t>& dict,std::vector<uint32_t>& block_stats)
        {
            const size_t dict_size = dict_entries * dict_entry_width;
            dict.resize(dict_size);
        }
    };

    template<uint32_t t_dict_entry_width = 16>
    struct dict_block_coder_greedy {
        static const uint32_t dict_entry_width = t_dict_entry_width;
        static const uint64_t block_size = 256;
        static const uint64_t overflow = 512; // dict coder can potentially overshoot...?
        static void encode(uint32_t const *dict,uint32_t const *in, uint32_t sum_of_values, size_t n,
                            std::vector<uint8_t> &out)
        {
            // TODO USE DICT TO ENCODE BLOCK
        }

        static uint8_t const *decode(uint32_t const *dict,uint8_t const *in, uint32_t *out,
                               uint32_t sum_of_values, size_t n)
        {
            // TODO USE DICT TO DECODE BLOCK
        }
    };

    template <typename DictBuilder,typename DictBlockCoder>
    class dict_freq_index {
    public:
        dict_freq_index()
            : m_size(0)
        {}

        class builder {
        public:
            builder(uint64_t num_docs, global_parameters const& params)
                : m_params(params)
            {
                m_num_docs = num_docs;
                m_endpoints.push_back(0);
            }

            template <typename DocsIterator, typename FreqsIterator>
            void add_posting_list(uint64_t n, DocsIterator docs_begin,
                                  FreqsIterator freqs_begin, uint64_t /* occurrences */)
            {
                if (!n) throw std::invalid_argument("List must be nonempty");
                dict_posting_list<DictBlockCoder>::write(m_doc_dict,m_freq_dict,m_lists, n,docs_begin, freqs_begin);
                m_endpoints.push_back(m_lists.size());
            }

            template <typename BlockDataRange>
            void add_posting_list(uint64_t n, BlockDataRange const& blocks)
            {
                if (!n) throw std::invalid_argument("List must be nonempty");
                dict_posting_list<DictBlockCoder>::write_blocks(m_lists, n, blocks);
                m_endpoints.push_back(m_lists.size());
            }

            template <typename BytesRange>
            void add_posting_list(BytesRange const& data)
            {
                m_lists.insert(m_lists.end(), std::begin(data), std::end(data));
                m_endpoints.push_back(m_lists.size());
            }

            void build(dict_freq_index& sq)
            {
                sq.m_params = m_params;
                sq.m_size = m_endpoints.size() - 1;
                sq.m_num_docs = m_num_docs;
                sq.m_lists.steal(m_lists);
                sq.m_doc_dict.steal(m_doc_dict);
                sq.m_freq_dict.steal(m_freq_dict);

                succinct::bit_vector_builder bvb;
                compact_elias_fano::write(bvb, m_endpoints.begin(),
                                          sq.m_lists.size(), sq.m_size,
                                          m_params); // XXX
                succinct::bit_vector(&bvb).swap(sq.m_endpoints);
            }

            void build_dict()
            {
                DictBuilder::build(m_doc_dict,m_block_docs_stats);
                DictBuilder::build(m_freq_dict,m_block_freqs_stats);
            }

            template <typename DocsIterator, typename FreqsIterator>
            void model_posting_list(uint64_t, DocsIterator, FreqsIterator, uint64_t)
            {
                // TODO gather m_block_docs_stats && m_block_freqs_stats using hashtable
            }

        private:
            global_parameters m_params;
            size_t m_num_docs;
            std::vector<uint64_t> m_endpoints;
            std::vector<uint8_t> m_lists;
            std::vector<uint32_t> m_doc_dict;
            std::vector<uint32_t> m_freq_dict;

            // TODO
            std::vector<uint32_t> m_block_docs_stats;
            std::vector<uint32_t> m_block_freqs_stats;
        };

        size_t size() const
        {
            return m_size;
        }

        uint64_t num_docs() const
        {
            return m_num_docs;
        }

        uint32_t dict_entries() const {
            return DictBuilder::dict_entries;
        }

        uint32_t dict_entry_width() const {
            return DictBuilder::dict_entry_width;
        }

        uint32_t block_size() const {
            return DictBlockCoder::block_size;
        }

        typedef typename dict_posting_list<DictBlockCoder>::document_enumerator document_enumerator;

        document_enumerator operator[](size_t i) const
        {
            assert(i < size());
            compact_elias_fano::enumerator endpoints(m_endpoints, 0,
                                                     m_lists.size(), m_size,
                                                     m_params);

            auto endpoint = endpoints.move(i).second;
            return document_enumerator(m_doc_dict.data(),m_freq_dict.data(),m_lists.data() + endpoint, num_docs(), i);
        }

        void warmup(size_t i) const
        {
            assert(i < size());
            compact_elias_fano::enumerator endpoints(m_endpoints, 0,
                                                     m_lists.size(), m_size,
                                                     m_params);

            auto begin = endpoints.move(i).second;
            auto end = m_lists.size();
            if (i + 1 != size()) {
                end = endpoints.move(i + 1).second;
            }

            volatile uint32_t tmp;
            for (size_t i = begin; i != end; ++i) {
                tmp = m_lists[i];
            }
            (void)tmp;
        }

        void swap(dict_freq_index& other)
        {
            std::swap(m_params, other.m_params);
            std::swap(m_size, other.m_size);
            m_endpoints.swap(other.m_endpoints);
            m_lists.swap(other.m_lists);
            m_doc_dict.swap(other.m_doc_dict);
            m_freq_dict.swap(other.m_freq_dict);
        }

        template <typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_params, "m_params")
                (m_size, "m_size")
                (m_num_docs, "m_num_docs")
                (m_endpoints, "m_endpoints")
                (m_lists, "m_lists")
                (m_doc_dict, "m_doc_dict")
                (m_freq_dict, "m_freq_dict")
                ;
        }

    private:
        global_parameters m_params;
        size_t m_size;
        size_t m_num_docs;
        uint32_t m_dict_entries;
        uint32_t m_dict_entry_width;
        succinct::bit_vector m_endpoints;
        succinct::mapper::mappable_vector<uint8_t> m_lists;
        succinct::mapper::mappable_vector<uint32_t> m_doc_dict;
        succinct::mapper::mappable_vector<uint32_t> m_freq_dict;
    };
}
