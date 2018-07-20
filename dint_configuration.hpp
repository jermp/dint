#pragma once

namespace ds2i {
    namespace constants {

        static const uint64_t min_size = 0; // 4096
        static const uint64_t max_size = 50000000;
        static const uint32_t block_size = 256;

        static const uint32_t top_k = 256; // 256

        // NOTE: if we use powers of 2 for the max_entry_size,
        // then num_target_sizes is log(max_entry_size) + 1.

        /********/

        // b = 16, l = 16
        static const uint32_t num_target_sizes = 5;
        static const uint32_t max_entry_size = 16;
        static const uint32_t target_sizes[] = {16, 8, 4, 2, 1};
        static const uint32_t num_entries = 65536;

        // b = 16, l = 8
        // static const uint32_t num_target_sizes = 4;
        // static const uint32_t max_entry_size = 8;
        // static const uint32_t target_sizes[] = {8, 4, 2, 1};
        // static const uint32_t num_entries = 65536;

        // b = 16, l = 4
        // static const uint32_t num_target_sizes = 3;
        // static const uint32_t max_entry_size = 4;
        // static const uint32_t target_sizes[] = {4, 2, 1};
        // static const uint32_t num_entries = 65536;

        /********/

        // b = 12, l = 16
        // static const uint32_t num_target_sizes = 5;
        // static const uint32_t max_entry_size = 16;
        // static const uint32_t target_sizes[] = {16, 8, 4, 2, 1};
        // static const uint32_t num_entries = 4096;

        // b = 12, l = 8
        // static const uint32_t num_target_sizes = 4;
        // static const uint32_t max_entry_size = 8;
        // static const uint32_t target_sizes[] = {8, 4, 2, 1};
        // static const uint32_t num_entries = 4096;

        // b = 12, l = 4
        // static const uint32_t num_target_sizes = 3;
        // static const uint32_t max_entry_size = 4;
        // static const uint32_t target_sizes[] = {4, 2, 1};
        // static const uint32_t num_entries = 4096;

        /********/

        // b = 8, l = 16
        // static const uint32_t num_target_sizes = 5;
        // static const uint32_t max_entry_size = 16;
        // static const uint32_t target_sizes[] = {16, 8, 4, 2, 1};
        // static const uint32_t num_entries = 256;

        // b = 8, l = 8
        // static const uint32_t num_target_sizes = 4;
        // static const uint32_t max_entry_size = 8;
        // static const uint32_t target_sizes[] = {8, 4, 2, 1};
        // static const uint32_t num_entries = 256;

        // b = 8, l = 4
        // static const uint32_t num_target_sizes = 3;
        // static const uint32_t max_entry_size = 4;
        // static const uint32_t target_sizes[] = {4, 2, 1};
        // static const uint32_t num_entries = 256;

        /********/

        static const double codeword_bits = 16.0;
        static const double initial_bpi = 48.0;
        static const double eps = 0.0001;
    }
}
