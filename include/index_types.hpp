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
#include "dint_codecs.hpp"
#include "dictionary_types.hpp"
#include "dictionary_builders.hpp"
#include "block_statistics.hpp"
#include "dict_freq_index.hpp"

namespace ds2i {

    typedef freq_index<
                compact_elias_fano,
                positive_sequence<strict_elias_fano>
            > ef_index;

    typedef freq_index<
                indexed_sequence,
                positive_sequence<>
            > single_index;

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

    // DINT indexes

    // collector type
    using adjusted_collector_type = adjusted<constants::max_entry_size>;

    // statistic types
    using adjusted_block_stats_type = block_statistics<adjusted_collector_type>;
    using adjusted_block_multi_stats_type = block_multi_statistics<adjusted_collector_type>;

    // dictionary_builders
    using single_rectangular_builder = decreasing_static_frequencies<
                                                single_dictionary_rectangular_type,
                                                adjusted_block_stats_type
                                            >;

    using single_packed_builder = decreasing_static_frequencies<
                                                single_dictionary_packed_type,
                                                adjusted_block_stats_type
                                            >;

    using multi_packed_builder = decreasing_static_frequencies<
                                                multi_dictionary_packed_type,
                                                adjusted_block_multi_stats_type
                                            >;

    // DINT configurations (all use optimal block parsing)
    using single_rect_dint_index    = dict_freq_index<single_rectangular_builder, opt_dint_single_dict_block>;
    using single_packed_dint_index  = dict_freq_index<single_packed_builder,      opt_dint_single_dict_block>;
    using multi_packed_dint_index   = dict_freq_index<multi_packed_builder,       opt_dint_multi_dict_block>;
}

#define DS2I_INDEX_TYPES (ef)(single)(uniform)(opt)(block_optpfor)(block_varintg8iu)(block_interpolative)(block_qmx)(block_mixed)(block_u32)(block_vbyte)(block_simple16)(block_varintgb)(block_maskedvbyte)(block_streamvbyte)(single_rect_dint)(single_packed_dint)(multi_packed_dint)
#define DS2I_BLOCK_INDEX_TYPES (block_optpfor)(block_varintg8iu)(block_interpolative)(block_qmx)(block_mixed)(block_u32)(block_vbyte)(block_simple16)(block_varintgb)(block_maskedvbyte)(block_streamvbyte)