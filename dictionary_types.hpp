#pragma once

#include "dint_configuration.hpp"
#include "dictionary_rectangular.hpp"
#include "dictionary_packed.hpp"
#include "dictionary_overlapped.hpp"

namespace ds2i {

    using large_dictionary_type = dictionary_overlapped // dictionary_rectangular dictionary_packed dictionary_overlapped
                                <constants::num_entries,
                                 constants::max_entry_size>;

    using small_dictionary_type = dictionary_overlapped
                                <256,
                                 large_dictionary_type::max_entry_size>;

}
