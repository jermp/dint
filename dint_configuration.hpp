#pragma once

namespace ds2i {
    namespace constants {
        static const uint64_t min_size = 0;
        static const uint64_t max_size = 50000000;

        static const uint32_t max_fractal_steps = 5;
        static const uint32_t max_entry_size = 16; // 8, 16
        static const uint32_t num_entries = 65536;

        static const double codeword_bits = 16.0;
        static const double initial_bpi = 48.0;
        static const double eps = 0.0001;
    }
}
