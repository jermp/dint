#pragma once

#include <succinct/mappable_vector.hpp>
#include <succinct/bit_vector.hpp>

#include "util.hpp"
#include "dictionary.hpp"
#include "dictionary_types.hpp"
#include "compact_elias_fano.hpp"
#include "dict_posting_list.hpp"
#include "block_statistics.hpp"

namespace ds2i {

    template<typename DictionaryStrategy,
             typename Encoder>
    struct dict_freq_index {
        using dict_strategy = DictionaryStrategy;
        using dict_type = typename dict_strategy::dict_type;
        using dict_builder_type = typename dict_strategy::builder_type;
        using encoder_type = Encoder;

        struct builder {
            builder(uint64_t num_docs, global_parameters const& params)
                : m_params(params)
            {
                m_num_docs = num_docs;
                m_endpoints.push_back(0);
            }

            template<typename DocsIterator, typename FreqsIterator>
            void add_posting_list(uint64_t n, DocsIterator docs_begin,
                                              FreqsIterator freqs_begin,
                                  uint64_t /* occurrences */)
            {
                if (!n) throw std::invalid_argument("List must be nonempty");
                dict_posting_list<dict_type,encoder_type>::write(m_docs_dict_builder,
                                                   m_freqs_dict_builder,
                                                   m_lists, n,
                                                   docs_begin, freqs_begin);
                m_endpoints.push_back(m_lists.size());
            }

            void build_model(std::string const& prefix_name)
            {
                {
                    DS2I_LOG << "building or loading dictionary for docs...";
                    build_or_load_dict<dict_strategy>(m_docs_dict_builder,prefix_name,dint_data_type::docs);
                }
                {
                    DS2I_LOG << "building or loading dictionary for freqs...";
                    build_or_load_dict<dict_strategy>(m_freqs_dict_builder,prefix_name,dint_data_type::freqs);
                }

                m_docs_dict_builder.prepare_for_encoding();
                m_freqs_dict_builder.prepare_for_encoding();
            }

            void build(dict_freq_index& dfi)
            {
                dfi.m_params = m_params;
                dfi.m_size = m_endpoints.size() - 1;
                dfi.m_num_docs = m_num_docs;
                dfi.m_lists.steal(m_lists);

                m_docs_dict_builder.build(dfi.m_docs_dict);
                m_freqs_dict_builder.build(dfi.m_freqs_dict);

                succinct::bit_vector_builder bvb;
                compact_elias_fano::write(bvb, m_endpoints.begin(),
                                          dfi.m_lists.size(),
                                          dfi.m_size, m_params);
                succinct::bit_vector(&bvb).swap(dfi.m_endpoints);
            }

        private:
            global_parameters m_params;
            size_t m_num_docs;
            std::vector<uint64_t> m_endpoints;
            std::vector<uint8_t> m_lists;

            dict_builder_type m_docs_dict_builder;
            dict_builder_type m_freqs_dict_builder;
        };

        dict_freq_index()
            : m_size(0)
        {}

        // # of terms in the collection
        size_t size() const {
            return m_size;
        }

        // # of docs in the collection
        uint64_t num_docs() const {
            return m_num_docs;
        }

        dict_type const& docs_dict() const {
            return m_docs_dict;
        }

        dict_type const& freqs_dict() const {
            return m_freqs_dict;
        }

        typedef typename dict_posting_list<dict_type,encoder_type>::document_enumerator document_enumerator;

        document_enumerator operator[](size_t i) const
        {
            assert(i < size());
            compact_elias_fano::enumerator endpoints(m_endpoints, 0,
                                                     m_lists.size(), m_size,
                                                     m_params);
            auto endpoint = endpoints.move(i).second;
            return document_enumerator(&m_docs_dict, &m_freqs_dict,
                                       m_lists.data() + endpoint,
                                       num_docs(), i);
        }

        void warmup(size_t i)
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
            m_docs_dict.swap(other.m_docs_dict);
            m_freqs_dict.swap(other.m_freqs_dict);
        }

        template<typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_params, "m_params")
                (m_size, "m_size")
                (m_num_docs, "m_num_docs")
                (m_endpoints, "m_endpoints")
                (m_lists, "m_lists")
                (m_docs_dict, "m_docs_dict")
                (m_freqs_dict, "m_freqs_dict")
                ;
        }

    private:
        global_parameters m_params;
        size_t m_size;
        size_t m_num_docs;
        succinct::bit_vector m_endpoints;
        succinct::mapper::mappable_vector<uint8_t> m_lists;
        dict_type m_docs_dict;
        dict_type m_freqs_dict;
    };
}

