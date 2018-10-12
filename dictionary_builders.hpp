#pragma once

#include <boost/progress.hpp>
#include <boost/filesystem.hpp>

#include <unordered_map>
#include <queue>

#include "dint_configuration.hpp"
#include "statistics_collectors.hpp"
#include "binary_blocks_collection.hpp"

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

    template<typename Dictionary,
             typename Statistics>
    struct decreasing_static_frequencies
    {
        typedef Dictionary dictionary_type;
        typedef Statistics statistics_type;

        static std::string type() {
            return "DSF-" + std::to_string(dictionary_type::num_entries) +
                      "-" + std::to_string(dictionary_type::max_entry_size);
        }

        static auto filter() {
            cost_filter filter(constants::eps / 1000);
            return filter;
        }

        static void build(typename dictionary_type::builder& dict_builder,
                          statistics_type& stats)
        {
            logger() << "building " << type() << " dictionary for " << stats.total_integers << std::endl;

            dict_builder.init();
            std::cout << "stats.blocks.size() " << stats.blocks.size() << std::endl;

            for (uint64_t s = 0; s != stats.blocks.size(); ++s)
            {
                uint64_t n = dictionary_type::num_entries;
                if (stats.blocks[s].size() < n) {
                    n = stats.blocks[s].size();
                }

                auto it = stats.blocks[s].begin();
                for (uint64_t i = 0; i != n; ++i, ++it) {
                    auto const& block = *it;
                    dict_builder.append(block.data.data(),
                                        block.data.size(), s);
                }
            }

            dict_builder.build();
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

            for (uint64_t s = 0; s != stats.blocks.size(); ++s)
            {
                for (auto& block: stats.blocks[s]) {
                    block.freq *= block.data.size();
                }

                freq_length_sorter sorter;
                std::sort(stats.blocks[s].begin(),
                          stats.blocks[s].end(),
                          sorter);

                uint64_t n = dictionary_type::num_entries;
                if (stats.blocks[s].size() < n) {
                    n = stats.blocks[s].size();
                }

                auto it = stats.blocks[s].begin();
                for (uint64_t i = 0; i < n; ++i, ++it) {
                    auto const& block = *it;
                    builder.append(block.data.data(),
                                   block.data.size(), s);
                }
            }

            builder.build();
        }
    };
}
