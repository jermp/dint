#pragma once

#include "dint_configuration.hpp"
#include "binary_collection.hpp"
#include "hash_utils.hpp"
#include "util.hpp"
#include "statistics_collectors.hpp"

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <unordered_map>

namespace ds2i {

    template<typename Collector>
    struct block_statistics
    {
        // static_assert(is_power_of_two(Collector::max_block_size));

        static std::string type() {
            return "block_statistics-" + std::to_string(Collector::max_block_size) + "-" + Collector::type();
        }

        template<typename Filter>
        static block_statistics create_or_load(std::string prefix_name, data_type dt,
                                               Filter const& filter)
        {
            std::string file_name = prefix_name + extension(dt);
            using namespace boost::filesystem;
            path p(file_name);
            std::string block_stats_filename = "./" + p.filename().string() + "." + type();

            if (boost::filesystem::exists(block_stats_filename)) {
                return block_statistics(block_stats_filename);
            }

            binary_collection input(file_name.c_str());
            bool compute_gaps = dt == data_type::docs;
            block_statistics stats(input, compute_gaps, filter);
            stats.try_to_store(block_stats_filename);
            return stats;
        }

        // create
        template<typename Filter>
        block_statistics(binary_collection& input, bool compute_gaps,
                         Filter const& filter)
        {
            logger() << "creating block stats (type = " << type() << ")" << std::endl;

            map_type block_map;
            // block_map.max_load_factor(0.01);
            boost::progress_display progress(input.num_postings());
            total_integers = 0;
            std::vector<uint32_t> buf;

            auto it = input.begin();
            if (compute_gaps) { // docs
                ++it; // skip first singleton sequence, containing # of docs
            }

            for (; it != input.end(); ++it) {
                auto const& list = *it;
                size_t n = list.size();
                total_integers += n;
                progress += n + 1;
                buf.reserve(n);
                uint32_t prev = compute_gaps ? -1 : 0;
                auto it = list.begin();
                for (uint32_t i = 0; i < n; ++i, ++it) {
                    buf.push_back(*it - prev - 1);
                    if (compute_gaps) {
                        prev = *it;
                    }
                }
                Collector::collect(buf, block_map);
                buf.clear();
            }

            logger() << "selecting entries..." << std::endl;
            blocks.reserve(block_map.size());
            for (auto& pair: block_map) {
                auto& freq_block = pair.second;
                if (filter(freq_block, total_integers)) {
                    block_type block;
                    block.freq = freq_block.freq;
                    block.data.swap(freq_block.data);
                    blocks.push_back(std::move(block));
                }
            }
            logger() << "DONE" << std::endl;
            logger() << "selected " << blocks.size() << " blocks" << std::endl;

            logger() << "sorting..." << std::endl;
            freq_sorter sorter;
            std::sort(blocks.begin(), blocks.end(), sorter);
            logger() << "DONE" << std::endl;
        }

        block_statistics(std::string file_name)
        {
            std::ifstream in(file_name.c_str());
            std::streamsize bytes = sizeof(uint32_t);
            uint32_t num_blocks, size, freq;
            in.read(reinterpret_cast<char*>(&total_integers), sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&num_blocks), bytes);
            logger() << "reading block stats (total_integers = " << total_integers
                     << "; num_blocks = " << num_blocks << ")" << std::endl;
            blocks.reserve(num_blocks);
            for (uint32_t i = 0; i < num_blocks; ++i) {
                in.read(reinterpret_cast<char*>(&size), bytes);
                in.read(reinterpret_cast<char*>(&freq), bytes);
                block_type block;
                block.freq = freq;
                block.data.resize(size);
                in.read(reinterpret_cast<char*>(block.data.data()), size * bytes);
                blocks.push_back(std::move(block));
            }
            logger() << "DONE" << std::endl;
        }

        void try_to_store(std::string file_name)
        {
            std::ofstream out(file_name.c_str());
            if (out) {
                logger() << "storing stats to disk..." << std::endl;
                out.write(reinterpret_cast<char const*>(&total_integers), sizeof(uint64_t));
                uint32_t num_blocks = blocks.size();
                std::streamsize bytes = sizeof(uint32_t);
                out.write(reinterpret_cast<char const*>(&num_blocks), bytes);
                for (auto const& block: blocks) {
                    uint32_t size = block.data.size();
                    uint32_t freq = block.freq;
                    out.write(reinterpret_cast<char const*>(&size), bytes);
                    out.write(reinterpret_cast<char const*>(&freq), bytes);
                    out.write(reinterpret_cast<char const*>(block.data.data()), size * bytes);
                }
                out.close();
                logger() << "DONE" << std::endl;
            } else {
                logger() << "cannot write block stats. collection directory not writeable" << std::endl;
            }
        }

        uint64_t total_integers;
        std::vector<block_type> blocks;
    };

}

