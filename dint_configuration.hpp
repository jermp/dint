#pragma once

#include <cmath>
#include <limits>

#define EXCEPTIONS 2
#define INF std::numeric_limits<uint32_t>::max()

namespace ds2i {

    namespace constants {

        static const uint64_t min_size = 0; // 4096
        static const uint64_t max_size = 50000000;
        static const uint32_t block_size = 256;

        enum block_selector {
            max = 0,
            median = 1,
            mode = 2
        };

        static const int context = block_selector::max;

        // static const uint32_t num_selectors = 16;
        // static const uint32_t selector_codes[] =
        //        {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 19, 22, 25, 29};
        // //                           8  10  12  14  16  19  22  25  29
        // //                           9  11  13  15  17  20  23  26  30
        // //                                          18  21  24  27  31
        // //                                                      28  32

        // static const uint32_t num_selectors = 5;
        // static const uint32_t selector_codes[] = {1, 2, 3, 4, 5};

        static const uint32_t num_selectors = 6;
        static const uint32_t selector_codes[] = {0, 1, 2, 3, 4, 5};

        // static const uint32_t num_selectors = 32;
        // static const uint32_t selector_codes[] =
        //        { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        //         16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        //         31, 32};
        // //                           8  10  12  14  16  19  22  25  29
        // //                           9  11  13  15  17  20  23  26  30
        // //                                          18  21  24  27  31
        // //                                                      28  32

        // NOTE: if we use powers of 2 for the max_entry_size,
        // then num_target_sizes is log(max_entry_size) + 1.

        /********/

        // b = 16, l = 16
        static const uint32_t max_entry_size = 16;
        static const uint32_t target_sizes[] = {16, 8, 4, 2, 1};
        static const uint32_t num_entries = 65536;
        static const uint32_t log2_num_entries = 16;
        static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        // b = 16, l = 8
        // static const uint32_t max_entry_size = 8;
        // static const uint32_t target_sizes[] = {8, 4, 2, 1};
        // static const uint32_t num_entries = 65536;
        // static const uint32_t log2_num_entries = 16;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        // b = 16, l = 4
        // static const uint32_t max_entry_size = 4;
        // static const uint32_t target_sizes[] = {4, 2, 1};
        // static const uint32_t num_entries = 65536;
        // static const uint32_t log2_num_entries = 16;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        /********/

        // b = 12, l = 16
        // static const uint32_t max_entry_size = 16;
        // static const uint32_t target_sizes[] = {16, 8, 4, 2, 1};
        // static const uint32_t num_entries = 4096;
        // static const uint32_t log2_num_entries = 16;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        // b = 12, l = 8
        // static const uint32_t max_entry_size = 8;
        // static const uint32_t target_sizes[] = {8, 4, 2, 1};
        // static const uint32_t num_entries = 4096;
        // static const uint32_t log2_num_entries = 16;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        // b = 12, l = 4
        // static const uint32_t max_entry_size = 4;
        // static const uint32_t target_sizes[] = {4, 2, 1};
        // static const uint32_t num_entries = 4096;
        // static const uint32_t log2_num_entries = 16;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        /********/

        // b = 8, l = 16
        // static const uint32_t max_entry_size = 16;
        // static const uint32_t target_sizes[] = {16, 8, 4, 2, 1};
        // static const uint32_t num_entries = 256;
        // static const uint32_t log2_num_entries = 8;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        // b = 8, l = 8
        // static const uint32_t max_entry_size = 8;
        // static const uint32_t target_sizes[] = {8, 4, 2, 1};
        // static const uint32_t num_entries = 256;
        // static const uint32_t log2_num_entries = 8;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        // b = 8, l = 4
        // static const uint32_t max_entry_size = 4;
        // static const uint32_t target_sizes[] = {4, 2, 1};
        // static const uint32_t num_entries = 256;
        // static const uint32_t log2_num_entries = 8;
        // static const uint32_t num_target_sizes = std::log2(max_entry_size) + 1;

        /********/

        static const double codeword_bits = 16.0;
        static const double initial_bpi = 48.0;
        static const double eps = 0.0001;
    }
}
