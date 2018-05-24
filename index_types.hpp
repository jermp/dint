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
#include "dictionary_types.hpp"
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

    typedef block_freq_index<ds2i::optpfor_block> block_optpfor_index;
    typedef block_freq_index<ds2i::varint_G8IU_block> block_varint_index;
    typedef block_freq_index<ds2i::interpolative_block> block_interpolative_index;
    typedef block_freq_index<ds2i::qmx_block> block_qmx_index;
    typedef block_freq_index<ds2i::mixed_block> block_mixed_index;
    typedef block_freq_index<ds2i::u32_block> block_u32_index;
    typedef block_freq_index<ds2i::vbyte_block> block_vbyte_index;
    typedef block_freq_index<ds2i::simple16_block> block_simple16_index;

    // DINT codec
    const uint32_t max_entry_width = 16;
    const uint32_t dict_entries = 65536;
    using block_stats_type = ds2i::block_statistics<max_entry_width,ds2i::stats_geometric>;
    using dict_type = ds2i::dint_dict_strategy_PDF<block_stats_type,dict_entries, max_entry_width>;
    using block_dint_index = dict_freq_index<dict_type,ds2i::dint_block>;
}

#define DS2I_INDEX_TYPES (ef)(single)(uniform)(opt)(block_optpfor)(block_varint)(block_interpolative)(block_mixed)(block_qmx)(block_u32)(block_vbyte)(block_simple16)(block_dint)
#define DS2I_BLOCK_INDEX_TYPES (block_optpfor)(block_varint)(block_interpolative)(block_qmx)(block_mixed)
