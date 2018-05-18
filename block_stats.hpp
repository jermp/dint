#pragma once

#include "binary_collection.hpp"

#include <boost/progress.hpp>

#include <unordered_map>

const size_t MIN_SIZE_THRES = 4096;

namespace ds2i {

    template<uint32_t max_entry_width>
    struct block_stats_full_stride_geom {

        static std::string type() {
            return "block_stats_full_stride_geom-L" + std::to_string(max_entry_width);
        }

        template<uint32_t width>
        struct block_info {
            size_t freq;
            size_t coverage;
            size_t savings;
            uint16_t entry_len;
            uint32_t entry[width];
        };

        using block_type = block_info<max_entry_width>;

        block_stats_full_stride_geom() {
        }

        block_stats_full_stride_geom(binary_collection& input,bool compute_gaps) {
            logger() << "creating block stats" << std::endl;
            std::unordered_map<uint32_t,uint64_t> block_map;
            block_map.max_load_factor(0.1);
            boost::progress_display progress(input.data_size());
            for (auto const& list: input) {
                size_t n = list.size();
                progress += n+1;
                if (n < MIN_SIZE_THRES) continue;
                process_list(block_map,list,compute_gaps);
            }
        }

        block_stats_full_stride_geom(std::string file_name) {
            std::ifstream in(file_name.c_str());
            uint64_t num_blocks;
            in.read(reinterpret_cast<char*>(&num_blocks), sizeof(uint64_t));
            logger() << "reading block stats (num_blocks = " << num_blocks << ")" << std::endl;
            blocks.resize(num_blocks);
            auto block_data = reinterpret_cast<char*>(blocks.data());
            in.read(block_data, num_blocks * sizeof(block_type));
        }

        template<class t_list>
        void process_list(std::unordered_map<uint32_t,uint64_t>& block_map,t_list& list,bool compute_gaps) {
            thread_local std::vector<uint32_t> buf(max_entry_width);
            size_t n = list.size();
            size_t full_strides = n / max_entry_width;
            auto lst_itr = list.begin();
            uint32_t prev = 0;
            for(size_t i=0;i<full_strides;i++) {
                // fill buf
                for(size_t j=0;j<max_entry_width;j++) {
                    buf[j] = *lst_itr - prev;
                    if(compute_gaps == true) prev = *lst_itr;
                    ++lst_itr;
                }
                // compute hashes (modified from murmur! hopefully still good...)
                const uint32_t m = 0x5bd1e995;
                const int r = 24;
                uint32_t hash = 0xDEADBEEF;
                size_t cur_len = 1;
                size_t cur_step = 0;
                for(size_t j=1;j<=max_entry_width;j++) {
                    uint32_t key = buf[j-1];
                    key *= m;
                    key ^= key >> r;
                    key *= m;
                    hash *= m;
                    hash ^= key;
                    if(j == cur_len) {
                        auto itr = block_map.find(hash);
                        if(itr != block_map.end()) {
                            auto block_idx = itr->second;
                            blocks[block_idx].freq += (max_entry_width >> cur_step);
                            blocks[block_idx].coverage = blocks[block_idx].freq * j;
                        } else {
                            // create a new block
                            block_type new_block;
                            new_block.freq = (max_entry_width >> cur_step);
                            new_block.coverage = new_block.freq * j;
                            new_block.entry_len = j;
                            for(size_t k=0;k<j;k++) {
                                new_block.entry[k] = buf[k];
                            }
                            block_map[hash] = blocks.size();
                            blocks.emplace_back(std::move(new_block));
                        }
                        cur_len *= 2;
                        cur_step++;
                    }
                }
            }
        }

        void try_to_store(std::string file_name) {
            logger() << "writing block stats" << std::endl;
            std::ofstream out(file_name.c_str());
            if(out) {
                uint64_t num_blocks = blocks.size();
                out.write(reinterpret_cast<char const*>(&num_blocks), sizeof(uint64_t));
                out.write(reinterpret_cast<char const*>(blocks.data()), num_blocks* sizeof(block_type));
            }
        }

        std::vector<block_type> blocks;
    };
}
