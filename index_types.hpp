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

    // DINT codecs
    const uint32_t max_entry_width = 16;
    const uint32_t num_entries = 65536;

    using block_stats_type = block_statistics<adjusted, // adjusted, full
                                              max_entry_width>;
    using dictionary_type = dictionary // rectangular, packed
                                <num_entries, max_entry_width>;

    using DSF = decreasing_static_frequencies<dictionary_type, block_stats_type>;
    using PDF = prefix_discounted_frequencies<dictionary_type, block_stats_type>;
    using LSS =     longest_to_shortest_sweep<dictionary_type, block_stats_type>;

    using DSF_block_dint_index = dict_freq_index<DSF, ds2i::dint_block>;
    using PDF_block_dint_index = dict_freq_index<PDF, ds2i::dint_block>;
    using LSS_block_dint_index = dict_freq_index<LSS, ds2i::dint_block>;
}

#define DS2I_INDEX_TYPES (ef)(single)(uniform)(opt)(block_optpfor)(block_varint)(block_interpolative)(block_mixed)(block_qmx)(block_u32)(block_vbyte)(block_simple16)(DSF_block_dint)(PDF_block_dint)(LSS_block_dint)
#define DS2I_BLOCK_INDEX_TYPES (block_optpfor)(block_varint)(block_interpolative)(block_qmx)(block_mixed)
