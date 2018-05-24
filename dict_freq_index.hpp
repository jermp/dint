#pragma once

#include <succinct/mappable_vector.hpp>
#include <succinct/bit_vector.hpp>

#include "util.hpp"
#include "dictionary.hpp"
#include "dictionary_builders.hpp"
#include "compact_elias_fano.hpp"
#include "dict_posting_list.hpp"
#include "blocks_statistics.hpp"

namespace ds2i {

    template<typename DictionaryBuilder,
             typename Encoder>
    struct dict_freq_index {

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
                dict_posting_list<Encoder>::write(&m_docs_dict_builder,
                                                  &m_freqs_dict_builder,
                                                   m_lists, n,
                                                   docs_begin, freqs_begin);
                m_endpoints.push_back(m_lists.size());
            }

            template<typename InputCollection>
            void build_model(InputCollection const& input, std::string const& prefix_name)
            {
                DS2I_LOG << "Building dictionary for docs...";
                DictBuilder::build(m_docs_dict_builder,input,dict_type::docs);

                DS2I_LOG << "Building dictionary for freqs...";
                DictBuilder::build(m_freqs_dict_builder,input,dict_type::freqs);

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

            dictionary::builder m_docs_dict_builder;
            dictionary::builder m_freqs_dict_builder;
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

        dictionary const& docs_dict() const {
            return m_docs_dict;
        }

        dictionary const& freqs_dict() const {
            return m_freqs_dict;
        }

        typedef typename dict_posting_list<Encoder>::document_enumerator document_enumerator;

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
        dictionary m_docs_dict;
        dictionary m_freqs_dict;
    };
}


                // logger() << "Collecting statistics..." << std::endl;

                // // step 1. collect statistics
                // std::vector<uint32_t> gaps;
                // uint64_t processed_lists = 0;
                // uint64_t total_integers = 0;

                // // for (uint32_t block_size = 1 /*16*/; block_size != 0; block_size /= 2)
                // // {
                // //     blocks_statistics docs_blocks_stats(block_size);
                // //     blocks_statistics freqs_blocks_stats(block_size);
                // //     total_integers = 0;
                // //     processed_lists = 0;

                // //     for (auto const& plist: input)
                // //     {
                // //         size_t n = plist.docs.size();
                // //         if (n > MIN_SIZE)
                // //         {
                // //             total_integers += n;
                // //             if (!n) throw std::invalid_argument("List must be nonempty");

                // //             gaps.reserve(n);
                // //             auto docs_begin = plist.docs.begin();
                // //             auto docs_end = docs_begin + n;
                // //             uint32_t prev = 0;
                // //             while (docs_begin != docs_end) {
                // //                 gaps.push_back(*docs_begin - prev);
                // //                 prev = *docs_begin;
                // //                 ++docs_begin;
                // //             }
                // //             assert(gaps.size() == n);

                // //             docs_blocks_stats.process(gaps.data(), n);
                // //             freqs_blocks_stats.process(plist.freqs.begin(), n); // do not take gaps

                // //             gaps.clear();
                // //             ++processed_lists;

                // //             if (processed_lists and processed_lists % 10000 == 0) {
                // //                 logger() << "processed " << processed_lists << " lists" << std::endl;
                // //                 logger() << "processed " << total_integers << " integers" << std::endl;
                // //             }
                // //         }
                // //     }

                // //     logger() << "processed " << processed_lists << " lists" << std::endl;
                // //     logger() << "processed " << total_integers << " integers" << std::endl;

                // //     // write blocks statistics to the disk
                // //     logger() << "Writing blocks statistics to the disk..." << std::endl;
                // //     std::string docs_output_filename("./docs.blocks_stats." + std::to_string(block_size) + ".bin");
                // //     docs_blocks_stats.sort_and_write(docs_output_filename);
                // //     std::string freqs_output_filename("./freqs.blocks_stats." + std::to_string(block_size) + ".bin");
                // //     freqs_blocks_stats.sort_and_write(freqs_output_filename);
                // // }

                // total_integers = 5406586692;
