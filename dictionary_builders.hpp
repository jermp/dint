#pragma once

#include "dictionary.hpp"

namespace ds2i {

    template<uint32_t t_dict_entries = 65536, uint32_t t_dict_entry_width = 16>
    struct dictionary_builder {
        static const uint32_t dict_entries = t_dict_entries;
        static const uint32_t dict_entry_width = t_dict_entry_width;

        static void build(dictionary& dict, std::vector<uint32_t>& block_stats)
        {
            dictionary::builder builder(dict_entries, dict_entry_width);

            // TODO: greedy algorithm to build the dictionary


            builder.build(dict);
        }

    };
}
