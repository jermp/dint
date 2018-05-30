#pragma once

#include "binary_collection.hpp"
#include "hash_utils.hpp"
#include "util.hpp"
#include "statistics_collectors.hpp"

#include <boost/progress.hpp>

#include <unordered_map>

namespace ds2i {

    template<typename Collector, uint32_t max_entry_width>
    struct block_statistics
    {
        static_assert(is_power_of_two(max_entry_width));
        using block_type = block_info<max_entry_width>;
        using bm_type = std::unordered_map<uint64_t, block_type>;

        static std::string type() {
            return "block_statistics-" + std::to_string(max_entry_width) + "-" + Collector::type();
        }

        static block_statistics create_or_load(std::string prefix_name, data_type dt) {
            std::string file_name = prefix_name + "." + extension(dt);
            std::string block_stats_file = file_name + "." + type();
            if (boost::filesystem::exists(block_stats_file)) {
                return block_statistics(block_stats_file);
            }
            binary_collection input(file_name.c_str());
            bool compute_gaps = dt == data_type::docs;
            block_statistics block_stats(input, compute_gaps);
            block_stats.try_to_store(block_stats_file);
            return block_stats;
        }

        block_statistics(binary_collection& input, bool compute_gaps) {
            DS2I_LOG << "creating block stats (type = " << type() << ")";
            bm_type block_map;
            block_map.max_load_factor(0.01);
            boost::progress_display progress(input.data_size());
            for (auto const& list: input) {
                size_t n = list.size();
                progress += n+1;
                process(list.begin(), n, compute_gaps, block_map);
            }

            DS2I_LOG << "serializing stats";
            blocks.resize(block_map.size());

            // Giulio: this will double the space at running time
            std::transform(block_map.begin(), block_map.end(),blocks.begin(), [](auto& pair){return pair.second;});

            DS2I_LOG << "sort by freq";
            std::sort(blocks.begin(),blocks.end(),[](const auto& a,const auto& b) {return a.freq > b.freq;});
            DS2I_LOG << "done creating stats";
        }

        block_statistics(std::string file_name) {
            std::ifstream in(file_name.c_str());
            uint64_t num_blocks;
            in.read(reinterpret_cast<char*>(&num_blocks), sizeof(uint64_t));
            DS2I_LOG << "reading block stats (num_blocks = " << num_blocks << ")";
            blocks.resize(num_blocks);
            auto block_data = reinterpret_cast<char*>(blocks.data());
            in.read(block_data, num_blocks * sizeof(block_type));
            DS2I_LOG << "done reading stats from disk";
        }

        template<typename Iterator>
        void process(Iterator it, size_t n, bool compute_gaps, bm_type& block_map)
        {
            thread_local std::vector<uint32_t> buf;
            buf.reserve(n);
            uint32_t prev = compute_gaps ? -1 : 0;
            for (uint32_t i = 0; i < n; ++i, ++it) {
                buf.push_back(*it - prev - 1);
                if (compute_gaps) {
                    prev = *it;
                }
            }
            assert(buf.size() == n);
            Collector::collect(buf, block_map);
        }

        void try_to_store(std::string file_name) {
            std::ofstream out(file_name.c_str());
            if(out) {
                DS2I_LOG << "writing block stats";
                uint64_t num_blocks = blocks.size();
                out.write(reinterpret_cast<char const*>(&num_blocks), sizeof(uint64_t));
                out.write(reinterpret_cast<char const*>(blocks.data()), num_blocks* sizeof(block_type));
                DS2I_LOG << "done writing stats to disk";
            } else {
                DS2I_LOG << "cannot write block stats. collection directory not writeable";
            }
        }

        // debug purposes
        std::string block_string(size_t id) const {
            const auto& b = blocks[id];
            std::string entry_str = "[";
            for(int i=0;i<b.entry_len-1;i++) {
                entry_str += std::to_string(b.entry[i]) + ",";
            }
            entry_str += std::to_string(b.entry[b.entry_len-1]) + "]";

            std::string str = "<id="  + std::to_string(id)
                + ",hash=" + std::to_string(b.hash)
                + ",freq=" + std::to_string(b.freq)
                + ",len=" + std::to_string(b.entry_len)
                + ",entry=" + entry_str + ">";
            return str;
        }

        std::vector<block_type> blocks;
    };

}

