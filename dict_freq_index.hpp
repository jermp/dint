#pragma once

#include <succinct/mappable_vector.hpp>
#include <succinct/bit_vector.hpp>

#include "util.hpp"
#include "dictionary.hpp"
#include "dictionary_builders.hpp"
#include "compact_elias_fano.hpp"
#include "dict_posting_list.hpp"
#include "block_statistics.hpp"

namespace ds2i {

    template<typename DictionaryBuilder, typename Coder>
    struct dict_freq_index
    {
        using dictionary_builder = DictionaryBuilder;
        using large_dictionary_type = typename dictionary_builder::large_dictionary_type;
        // using small_dictionary_type = typename dictionary_builder::small_dictionary_type;
        using coder_type = Coder;

        typedef dict_posting_list<large_dictionary_type/*, small_dictionary_type*/, coder_type> sequence_type;

        struct builder {
            builder(uint64_t num_docs, global_parameters const& params)
                : m_num_docs(num_docs)
                , m_params(params)
                , m_queue(1 << 24)
                , m_docs_large_dict_builders(constants::num_selectors)
                , m_freqs_large_dict_builders(constants::num_selectors)
            {
                m_endpoints.push_back(0);
            }

            template<typename DocsIterator, typename FreqsIterator>
            void add_posting_list(uint64_t n, DocsIterator docs_begin,
                                              FreqsIterator freqs_begin,
                                  uint64_t /* occurrences */)
            {
                if (!n) throw std::invalid_argument("List must be nonempty");

                // sequence_type::write(m_docs_large_dict_builder, m_freqs_large_dict_builder,
                //                      m_docs_small_dict_builder, m_freqs_small_dict_builder,
                //                      m_lists, n, docs_begin, freqs_begin);
                // m_endpoints.push_back(m_lists.size());

                std::shared_ptr<list_adder<DocsIterator, FreqsIterator>>
                    ptr(new list_adder<DocsIterator, FreqsIterator>(*this, docs_begin, freqs_begin, n));
                m_queue.add_job(ptr, 2 * n);
            }

            void build_model(std::string const& prefix_name)
            {
                logger() << "building or loading dictionaries for docs..." << std::endl;
                build_or_load_dict(m_docs_large_dict_builders, prefix_name, data_type::docs);
                logger() << "DONE" << std::endl;

                logger() << "building or loading dictionaries for freqs..." << std::endl;
                build_or_load_dict(m_freqs_large_dict_builders, prefix_name, data_type::freqs);
                logger() << "DONE" << std::endl;

                for (int s = 0; s != constants::num_selectors; ++s) {
                    m_docs_large_dict_builders[s].prepare_for_encoding();
                    m_freqs_large_dict_builders[s].prepare_for_encoding();
                }
            }

            void build(dict_freq_index& dfi)
            {
                m_queue.complete();

                dfi.m_params = m_params;
                dfi.m_size = m_endpoints.size() - 1;
                dfi.m_num_docs = m_num_docs;
                dfi.m_lists.steal(m_lists);

                // std::cout << "BLOCK STATISTICS" << std::endl;
                // uint64_t total_blocks = m_docs_large_dict_builder.block_encoded_with_large_dict
                //                       + m_docs_large_dict_builder.block_encoded_with_small_dict;

                // std::cout << "blocks encoded with large dict on docs: "
                //           << m_docs_large_dict_builder.block_encoded_with_large_dict << "/" << total_blocks << " "
                //           << "(" << m_docs_large_dict_builder.block_encoded_with_large_dict * 100.0 / total_blocks << "%)" << std::endl;

                // std::cout << "blocks encoded with small dict on docs: "
                //           << m_docs_large_dict_builder.block_encoded_with_small_dict << "/" << total_blocks << " "
                //           << "(" << m_docs_large_dict_builder.block_encoded_with_small_dict * 100.0 / total_blocks << "%)" << std::endl;

                // total_blocks = m_freqs_large_dict_builder.block_encoded_with_large_dict
                //              + m_freqs_large_dict_builder.block_encoded_with_small_dict;

                // std::cout << "blocks encoded with large dict on freqs: "
                //           << m_freqs_large_dict_builder.block_encoded_with_large_dict << "/" << total_blocks << " "
                //           << "(" << m_freqs_large_dict_builder.block_encoded_with_large_dict * 100.0 / total_blocks << "%)" << std::endl;

                // std::cout << "blocks encoded with small dict on freqs: "
                //           << m_freqs_large_dict_builder.block_encoded_with_small_dict << "/" << total_blocks << " "
                //           << "(" << m_freqs_large_dict_builder.block_encoded_with_small_dict * 100.0 / total_blocks << "%)" << std::endl;


                // std::cout << "LARGE DICTIONARY" << std::endl;
                // std::cout << "docs codewords: " << m_docs_large_dict_builder.codewords << std::endl;
                // std::cout << "docs small exceptions: " << m_docs_large_dict_builder.small_exceptions << std::endl;
                // std::cout << "docs large exceptions: " << m_docs_large_dict_builder.large_exceptions << std::endl;
                // std::cout << "freqs codewords: " << m_freqs_large_dict_builder.codewords << std::endl;
                // std::cout << "freqs small exceptions: " << m_freqs_large_dict_builder.small_exceptions << std::endl;
                // std::cout << "freqs large exceptions: " << m_freqs_large_dict_builder.large_exceptions << std::endl;

                // std::cout << "SMALL DICTIONARY" << std::endl;
                // std::cout << "docs codewords: " << m_docs_small_dict_builder.codewords << std::endl;
                // std::cout << "docs small exceptions: " << m_docs_small_dict_builder.small_exceptions << std::endl;
                // std::cout << "docs large exceptions: " << m_docs_small_dict_builder.large_exceptions << std::endl;
                // std::cout << "freqs codewords: " << m_freqs_small_dict_builder.codewords << std::endl;
                // std::cout << "freqs small exceptions: " << m_freqs_small_dict_builder.small_exceptions << std::endl;
                // std::cout << "freqs large exceptions: " << m_freqs_small_dict_builder.large_exceptions << std::endl;

                for (int s = 0; s != constants::num_selectors; ++s) {
                    m_docs_large_dict_builders[s].build(dfi.m_docs_large_dicts[s]);
                    m_freqs_large_dict_builders[s].build(dfi.m_freqs_large_dicts[s]);
                }

                // m_docs_large_dict_builder.build(dfi.m_docs_large_dict);
                // m_freqs_large_dict_builder.build(dfi.m_freqs_large_dict);

                // m_docs_small_dict_builder.build(dfi.m_docs_small_dict);
                // m_freqs_small_dict_builder.build(dfi.m_freqs_small_dict);


                succinct::bit_vector_builder bvb;
                compact_elias_fano::write(bvb, m_endpoints.begin(),
                                          dfi.m_lists.size(),
                                          dfi.m_size, m_params);
                succinct::bit_vector(&bvb).swap(dfi.m_endpoints);
            }

        private:

            template<typename DocsIterator, typename FreqsIterator>
            struct list_adder : semiasync_queue::job {
                list_adder(builder& b,
                           DocsIterator docs_begin,
                           FreqsIterator freqs_begin,
                           uint64_t n)
                    : b(b)
                    , docs_begin(docs_begin)
                    , freqs_begin(freqs_begin)
                    , n(n)
                {}

                virtual void prepare() {
                    sequence_type::write(b.m_docs_large_dict_builders, b.m_freqs_large_dict_builders,
                                         lists, n, docs_begin, freqs_begin);
                }

                virtual void commit() {
                    b.m_lists.insert(b.m_lists.end(), lists.begin(), lists.end());
                    b.m_endpoints.push_back(b.m_lists.size());
                }

                builder& b;
                std::vector<uint8_t> lists;
                DocsIterator docs_begin;
                FreqsIterator freqs_begin;
                uint64_t n;
            };

            uint64_t m_num_docs;
            global_parameters m_params;
            semiasync_queue m_queue;
            std::vector<uint64_t> m_endpoints;
            std::vector<uint8_t> m_lists;
            std::vector<typename large_dictionary_type::builder> m_docs_large_dict_builders;
            std::vector<typename large_dictionary_type::builder> m_freqs_large_dict_builders;

            // typename large_dictionary_type::builder m_docs_large_dict_builder;
            // typename large_dictionary_type::builder m_freqs_large_dict_builder;
            // typename small_dictionary_type::builder m_docs_small_dict_builder;
            // typename small_dictionary_type::builder m_freqs_small_dict_builder;

            void build_or_load_dict(
                                    std::vector<typename large_dictionary_type::builder>& large_dict_builders,
                                    // typename large_dictionary_type::builder& large_dict_builder,
                                    // typename small_dictionary_type::builder& small_dict_builder,
                                    std::string prefix_name, data_type dt)
            {
                std::string file_name = prefix_name + extension(dt);
                using namespace boost::filesystem;
                path p(file_name);
                using ld_type = typename large_dictionary_type::builder;

                bool already_built = true;
                std::vector<std::string> filenames;
                for (int s = 0; s != constants::num_selectors; ++s) {
                    uint32_t selector_code = constants::selector_codes[s];
                    std::string dictionary_filename = "./dict."
                                                + p.filename().string() + "."
                                                + ld_type::type() + "."
                                                + dictionary_builder::type() + "."
                                                + std::to_string(selector_code)
                                                ;
                    already_built = already_built and boost::filesystem::exists(dictionary_filename);
                    filenames.push_back(dictionary_filename);
                }

                // std::cout << "***" << std::endl;
                // for (auto const& s: filenames) {
                //     std::cout << s << std::endl;
                // }

                if (already_built) {
                    for (int s = 0; s != constants::num_selectors; ++s) {
                        large_dict_builders[s].load_from_file(filenames[s]);
                    }
                } else {
                    using statistics_type = typename dictionary_builder::statistics_type;
                    auto statistics =
                        statistics_type::create_or_load(prefix_name, dt,
                                                        dictionary_builder::filter());
                    dictionary_builder::build(large_dict_builders, statistics, dt);
                    for (int s = 0; s != constants::num_selectors; ++s) {
                        if (not large_dict_builders[s].try_store_to_file(filenames[s])) {
                            logger() << "cannot write dictionary to file '" << filenames[s] << "'" << std::endl;
                        }
                    }
                }

                // builder.print_entries(dictionary_file + ".log");
                // builder.print_indexes(dictionary_file + ".rbo_input.estimated");
            }
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

        // large_dictionary_type const& docs_large_dict() const {
        //     return m_docs_large_dict;
        // }

        // large_dictionary_type const& freqs_large_dict() const {
        //     return m_freqs_large_dict;
        // }

        // small_dictionary_type const& docs_small_dict() const {
        //     return m_docs_small_dict;
        // }

        // small_dictionary_type const& freqs_small_dict() const {
        //     return m_freqs_small_dict;
        // }

        typedef typename dict_posting_list<large_dictionary_type/*, small_dictionary_type*/, coder_type>::document_enumerator document_enumerator;

        document_enumerator operator[](size_t i) const
        {
            assert(i < size());
            compact_elias_fano::enumerator endpoints(m_endpoints, 0,
                                                     m_lists.size(), m_size,
                                                     m_params);
            auto endpoint = endpoints.move(i).second;
            return document_enumerator(
                                       // &m_docs_large_dict, &m_freqs_large_dict,
                                       // &m_docs_small_dict, &m_freqs_small_dict,
                                       m_docs_large_dicts, m_freqs_large_dicts,

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

        void swap(dict_freq_index& other) {
            std::swap(m_params, other.m_params);
            std::swap(m_size, other.m_size);
            m_endpoints.swap(other.m_endpoints);
            m_lists.swap(other.m_lists);

            for (int s = 0; s != constants::num_selectors; ++s) {
                m_docs_large_dicts[s].swap(other.m_docs_large_dicts[s]);
                m_freqs_large_dicts[s].swap(other.m_freqs_large_dicts[s]);
            }

            // m_docs_large_dict.swap(other.m_docs_large_dict);
            // m_freqs_large_dict.swap(other.m_freqs_large_dict);
            // m_docs_small_dict.swap(other.m_docs_small_dict);
            // m_freqs_small_dict.swap(other.m_freqs_small_dict);
        }

        template<typename Visitor>
        void map(Visitor& visit) {
            visit
                (m_params, "m_params")
                (m_size, "m_size")
                (m_num_docs, "m_num_docs")
                (m_endpoints, "m_endpoints")
                (m_lists, "m_lists")
                // (m_docs_large_dict, "m_docs_large_dict")
                // (m_freqs_large_dict, "m_freqs_large_dict")
                // (m_docs_small_dict, "m_docs_small_dict")
                // (m_freqs_small_dict, "m_freqs_small_dict")
                ;

            for (int s = 0; s != constants::num_selectors; ++s) {
                std::string docs_dict_name = "m_docs_large_dicts[" + std::to_string(s) + "]";
                std::string freqs_dict_name = "m_freqs_large_dicts[" + std::to_string(s) + "]";
                visit
                    (m_docs_large_dicts[s], docs_dict_name.c_str())
                    (m_freqs_large_dicts[s], freqs_dict_name.c_str())
                    ;
            }

        }

    private:
        global_parameters m_params;
        size_t m_size;
        size_t m_num_docs;
        succinct::bit_vector m_endpoints;
        succinct::mapper::mappable_vector<uint8_t> m_lists;

        std::vector<large_dictionary_type> m_docs_large_dicts;
        std::vector<large_dictionary_type> m_freqs_large_dicts;

        // large_dictionary_type m_docs_large_dict;
        // large_dictionary_type m_freqs_large_dict;
        // small_dictionary_type m_docs_small_dict;
        // small_dictionary_type m_freqs_small_dict;
    };
}

