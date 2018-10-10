#pragma once

#include "dint_configuration.hpp"

namespace ds2i {

    double cost(uint32_t block_size, uint32_t block_frequency) {
        return block_frequency * (constants::initial_bpi * block_size
                               -  constants::codeword_bits);
    }

    double compute_saving(uint32_t block_size, uint32_t block_frequency, uint64_t total_integers) {
        return cost(block_size, block_frequency) / total_integers;
    };

}
