#pragma once

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>

#include "freq_index.hpp"
#include "positive_sequence.hpp"
#include "partitioned_sequence.hpp"
#include "uniform_partitioned_sequence.hpp"
#include "binary_freq_collection.hpp"
#include "block_freq_index.hpp"
#include "block_codecs.hpp"
#include "mixed_block.hpp"

#include "dint_configuration.hpp"
#include "dictionary_builders.hpp"
#include "block_statistics.hpp"
#include "dict_freq_index.hpp"

namespace ds2i {

    typedef freq_index<compact_elias_fano,
                       positive_sequence<strict_elias_fano>> ef_index;

    typedef freq_index<indexed_sequence,
                       positive_sequence<>> single_index;

    typedef freq_index<
        uniform_partitioned_sequence<>,
        positive_sequence<uniform_partitioned_sequence<strict_sequence>>
        > uniform_index;

    typedef freq_index<
        partitioned_sequence<>,
        positive_sequence<partitioned_sequence<strict_sequence>>
        > opt_index;

    typedef block_freq_index<optpfor_block> block_optpfor_index;
    typedef block_freq_index<varint_G8IU_block> block_varintg8iu_index;
    typedef block_freq_index<interpolative_block> block_interpolative_index;
    typedef block_freq_index<qmx_block> block_qmx_index;
    typedef block_freq_index<mixed_block> block_mixed_index;
    typedef block_freq_index<u32_block> block_u32_index;
    typedef block_freq_index<vbyte_block> block_vbyte_index;
    typedef block_freq_index<simple16_block> block_simple16_index;
    typedef block_freq_index<varintgb_block> block_varintgb_index;
    typedef block_freq_index<maskedvbyte_block> block_maskedvbyte_index;
    typedef block_freq_index<streamvbyte_block> block_streamvbyte_index;

    // DINT codecs

    using large_dictionary_type = dictionary
                                <constants::num_entries,
                                 constants::max_entry_size>;
    using small_dictionary_type = dictionary
                                <256,
                                 large_dictionary_type::max_entry_size>;

    using adjusted_collector_type = adjusted<large_dictionary_type::max_entry_size>;
    using     full_collector_type =     full<large_dictionary_type::max_entry_size>;
    using    fixed_collector_type =    fixed<large_dictionary_type::max_entry_size>;

    using adjusted_block_stats_type = block_statistics<adjusted_collector_type>;
    using     full_block_stats_type = block_statistics<full_collector_type>;
    using    fixed_block_stats_type = block_statistics<fixed_collector_type>;

    using DSF = decreasing_static_frequencies<large_dictionary_type, small_dictionary_type, adjusted_block_stats_type>;
    // using DSV = decreasing_static_volume<large_dictionary_type, adjusted_block_stats_type>;
    // using PDF = prefix_discounted_frequencies<large_dictionary_type, full_block_stats_type>;
    // using LSS =     longest_to_shortest_sweep<large_dictionary_type, adjusted_block_stats_type>;
    // using LSO =             long_strings_only<large_dictionary_type, fixed_block_stats_type>;

    using DSF_block_dint_index = dict_freq_index<DSF, dint_block>;
    // using DSV_block_dint_index = dict_freq_index<DSV, dint_block>;
    // using PDF_block_dint_index = dict_freq_index<PDF, dint_block>;
    // using LSS_block_dint_index = dict_freq_index<LSS, dint_block>;
    // using LSO_block_dint_index = dict_freq_index<LSO, dint_block>;
}

#define DS2I_INDEX_TYPES (ef)(single)(uniform)(opt)(block_optpfor)(block_varintg8iu)(block_interpolative)(block_qmx)(block_mixed)(block_u32)(block_vbyte)(block_simple16)(block_varintgb)(block_maskedvbyte)(block_streamvbyte)(DSF_block_dint)//(DSV_block_dint)(PDF_block_dint)(LSS_block_dint)(LSO_block_dint)
#define DS2I_BLOCK_INDEX_TYPES (block_optpfor)(block_varintg8iu)(block_interpolative)(block_qmx)(block_mixed)(block_u32)(block_vbyte)(block_simple16)(block_varintgb)(block_maskedvbyte)(block_streamvbyte)
