#pragma once

#include "dint_configuration.hpp"

namespace ds2i {

    double cost(uint32_t block_size, uint32_t block_frequency) {
        return block_frequency * (initial_bpis * block_size - codeword_bits);
    }

    double bpi(uint32_t block_size, uint32_t block_frequency, uint64_t total_integers) {
        return cost(block_size, block_frequency) / total_integers;
    };

}
