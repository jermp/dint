#pragma once

#include "binary_collection.hpp"
#include "hash_utils.hpp"
#include "util.hpp"

#include <boost/progress.hpp>

#include <unordered_map>

const size_t MIN_SIZE_THRES = 256;

namespace ds2i {
    namespace util {
        constexpr bool is_power_of_two(uint64_t x) {
            return (x & (x - 1)) == 0;
        }
    }

    #pragma pack(push, 1)
    template<uint32_t width>
    struct block_info {
        uint64_t hash;
        uint64_t freq;
        uint8_t entry_len;
        uint32_t entry[width];
    };
    #pragma pack(pop)

    template<class bm_type>
    void update_entry(const uint32_t* entry,size_t n,bm_type& block_map)
    {
    	using b_type = typename bm_type::value_type::second_type;
    	auto hash = hash_u32(entry,n);
		auto itr = block_map.find(hash);
		if(itr != block_map.end())
			itr->second.freq++;
		else {
			b_type new_block;
			new_block.hash = hash;
			new_block.freq = 1;
			new_block.entry_len = n;
			std::copy(entry,entry+n,new_block.entry);
			block_map[hash] = new_block;
		}
    }

    struct stats_geometric {
    	static std::string type() {
    		return "geometric";
    	}

    	template<class bm_type>
    	static void update_stats(const std::vector<uint32_t>& buf,bm_type& block_map) {
    		auto b = buf.data();
    		for(size_t size_u32=1;size_u32<=buf.size();size_u32*=2) {
    			update_entry(b,size_u32,block_map);
    		}
    	}
    };

    struct stats_adjusted {
    	static std::string type() {
    		return "adjusted";
    	}

    	template<class bm_type>
    	static void update_stats(std::vector<uint32_t>& buf,bm_type& block_map) {
    		auto b = buf.data();
    		for(size_t size_u32=1;size_u32<=buf.size();size_u32*=2) {
    			for(size_t pos=0;pos<buf.size();pos+=size_u32) {
    				update_entry(b+pos,size_u32,block_map);
    			}
    		}
    	}
    };

    template<uint32_t max_entry_width,class stats_type>
    struct block_statistics {
    	static_assert(util::is_power_of_two(max_entry_width));
    	using block_type = block_info<max_entry_width>;
    	using bm_type = std::unordered_map<uint64_t,block_type>;

        static std::string type() {
            return "block_stats-L" + std::to_string(max_entry_width) + "-" + stats_type::type();
        }

    	block_statistics(binary_collection& input,bool compute_gaps) {
            DS2I_LOG << "creating block stats (type = " << type() 
                     << ") - MIN_SIZE_THRES = " << MIN_SIZE_THRES;
            bm_type block_map;
            block_map.max_load_factor(0.01);
            boost::progress_display progress(input.data_size());
            for (auto const& list: input) {
                size_t n = list.size();
                progress += n+1;
                if (n < MIN_SIZE_THRES) continue;
                process_list(block_map,list,compute_gaps);
            }

            DS2I_LOG << "serializing stats";
            blocks.resize(block_map.size());
            auto block_selector = [](auto& pair){return pair.second;};
            std::transform(block_map.begin(), block_map.end(),blocks.begin(), block_selector);
            DS2I_LOG << "sort by freq";
            auto freq_cmp = [](const auto& a,const auto& b) {return a.freq > b.freq;};
            std::sort(blocks.begin(),blocks.end(),freq_cmp);
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

        template<class t_list>
    	void process_list(bm_type& block_map,t_list& list,bool compute_gaps) {
    		thread_local std::vector<uint32_t> buf(max_entry_width);
            size_t n = list.size();
            size_t full_strides = n / max_entry_width;
            auto lst_itr = list.begin();
            uint32_t prev = 0;
            for(size_t i=0;i<full_strides;i++) {
                // (1) fill buf
                for(size_t j=0;j<max_entry_width;j++) {
                    buf[j] = *lst_itr - prev;
                    if(compute_gaps == true) prev = *lst_itr;
                    ++lst_itr;
                }
                // (2) process buf and update freq stats
                stats_type::update_stats(buf,block_map);
            }
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

