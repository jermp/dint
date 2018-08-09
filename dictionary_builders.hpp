#pragma once

#include <boost/progress.hpp>
#include <boost/filesystem.hpp>

#include <unordered_map>
#include <queue>

#include "dint_configuration.hpp"
#include "statistics_collectors.hpp"
#include "binary_blocks_collection.hpp"
#include "dictionary.hpp"
#include "heap.hpp"
#include "hash_utils.hpp"
#include "util.hpp"
#include "model_build_utils.hpp"

namespace ds2i {

    struct cost_filter {
        cost_filter(double threshold = constants::eps)
            : m_threshold(threshold)
        {}

        bool operator()(block_type const& block, uint64_t total_integers) const {
            return compute_saving(block.data.size(),
                                  block.freq,
                                  total_integers) > m_threshold;
        }

    private:
        double m_threshold;
    };

    template<typename LargeDictionary,
             typename SmallDictionary,
             typename Statistics>
    struct decreasing_static_frequencies
    {
        typedef LargeDictionary dictionary_type;
        typedef LargeDictionary large_dictionary_type;
        typedef SmallDictionary small_dictionary_type;
        typedef Statistics statistics_type;

        static std::string type() {
            return "DSF-" + std::to_string(large_dictionary_type::num_entries) +
                      "-" + std::to_string(large_dictionary_type::max_entry_size);
        }

        static auto filter() {
            cost_filter filter(constants::eps / 1000);
            return filter;
        }

        static void build(typename large_dictionary_type::builder& dict_builder,
                          statistics_type& stats)
        {
            logger() << "building " << type() << " dictionary for " << stats.total_integers << std::endl;

            dict_builder.init();
            uint64_t n = large_dictionary_type::num_entries;
            std::cout << "stats.blocks.size() " << stats.blocks.size() << std::endl;
            if (stats.blocks.size() < n) {
                n = stats.blocks.size();
            }

            // NOTE: sort by decreasing length
            std::sort(stats.blocks.begin(), stats.blocks.begin() + n,
                [](auto const& l, auto const& r) {
                    return l.data.size() > r.data.size();
                });

            auto it = stats.blocks.begin();
            for (uint64_t i = 0; i < n; ++i, ++it) {
                auto const& block = *it;
                dict_builder.append(block.data.data(),
                                    block.data.size());
            }

            logger() << "DONE" << std::endl;
        }

        static void build(std::vector<typename large_dictionary_type::builder>& large_dict_builders,
                          std::vector<typename small_dictionary_type::builder>& small_dict_builders,
                          Statistics& stats)
        {
            logger() << "building " << type() << " dictionaries for " << stats.total_integers << std::endl;

            for (int s = 0; s != constants::num_selectors; ++s)
            {
                // small dict
                {
                    small_dict_builders[s].init();
                    uint64_t n = small_dictionary_type::num_entries;
                    if (stats.blocks[s].size() < n) {
                        n = stats.blocks[s].size();
                    }

                    // NOTE: sort by decreasing length
                    std::sort(stats.blocks[s].begin(),
                              stats.blocks[s].begin() + n,
                        [](auto const& l, auto const& r) {
                            return l.data.size() > r.data.size();
                        });

                    auto it = stats.blocks[s].begin();
                    for (uint64_t i = 0; i < n; ++i, ++it) {
                        auto const& block = *it;
                        small_dict_builders[s].append(block.data.data(),
                                                      block.data.size());
                    }
                }

                // large dict
                {
                    large_dict_builders[s].init();
                    uint64_t n = large_dictionary_type::num_entries;
                    if (stats.blocks[s].size() < n) {
                        n = stats.blocks[s].size();
                    }

                    auto it = stats.blocks[s].begin();
                    // std::string file_name = "dict." + type() + "." + std::to_string(constants::selector_codes[s]);
                    // std::ofstream out(file_name.c_str());
                    for (uint64_t i = 0; i < n; ++i, ++it) {
                        auto const& block = *it;

                        // out << i << ": entry [";
                        // for (uint64_t k = 0; k != block.data.size(); ++k) {
                        //     out << block.data[k];
                        //     if (k == block.data.size() - 1) {
                        //         out << "] ";
                        //     } else {
                        //         out << ", ";
                        //     }
                        // }
                        // out << "; freq " << block.freq << "\n";

                        large_dict_builders[s].append(block.data.data(),
                                                      block.data.size());
                    }
                    // out.close();
                }
            }

            logger() << "DONE" << std::endl;
            // std::abort();
        }
    };

    template<typename Dictionary, typename Statistics>
    struct decreasing_static_volume
    {
        typedef Dictionary dictionary_type;
        typedef Statistics statistics_type;

        static std::string type() {
            return "DSV-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static auto filter() {
            cost_filter filter(constants::eps / 1000);
            return filter;
        }

        static void build(typename dictionary_type::builder& builder,
                          statistics_type& stats)
        {
            logger() << "building " << type() << " dictionary for " << stats.total_integers << std::endl;
            builder.init();

            for (auto& block: stats.blocks) {
                block.freq *= block.data.size();
            }

            freq_length_sorter sorter;
            std::sort(stats.blocks.begin(),
                      stats.blocks.end(),
                      sorter);

            uint64_t n = dictionary_type::num_entries;
            if (stats.blocks.size() < n) {
                n = stats.blocks.size();
            }

            auto it = stats.blocks.begin();
            for (uint64_t i = 0; i < n; ++i, ++it) {
                auto const& block = *it;
                builder.append(block.data.data(),
                               block.data.size());
            }

            logger() << "DONE" << std::endl;
        }
    };
}
