#pragma once

#include "dictionary_types.hpp"
#include "statistics_collectors.hpp"

namespace ds2i {

struct dint_statistics {
    dint_statistics()
        : ints_distr(7, 0)
        , codewords_distr(7, 0)
        , occs(constants::num_entries, 0) {}

    std::vector<uint64_t>
        ints_distr;  // 0:runs; 1:1; 2:2; 3:4; 4:8; 5:16; 6:exceptions
    std::vector<uint64_t>
        codewords_distr;  // 0:runs; 1:1; 2:2; 3:4; 4:8; 5:16; 6:exceptions

    std::vector<uint64_t> occs;

    uint64_t decoded_ints_from_dict = 0;
    uint64_t dict_codewords = 0;
    uint64_t total_ints = 0;

    std::unordered_map<uint32_t, uint32_t> exceptions;

    void eat(uint32_t exception) {
        auto it = exceptions.find(exception);
        if (it == exceptions.end()) {
            exceptions[exception] = 1;
        } else {
            ++exceptions[exception];
        }
    }
};

struct single_dint {
    template <typename Dictionary>
    static uint8_t const* decode(Dictionary const& dict, uint8_t const* in,
                                 uint32_t* out, size_t n
                                 // , dint_statistics& stats
    ) {
        uint16_t const* ptr = reinterpret_cast<uint16_t const*>(in);
        // uint8_t const* ptr = in; // if b = 8
        for (size_t i = 0; i != n; ++ptr) {
            uint32_t index = *ptr;
            uint32_t decoded_ints = 1;

            // ++stats.occs[index];

            if (DS2I_LIKELY(index > EXCEPTIONS - 1)) {
                decoded_ints = dict.copy(index, out);

                // if (decoded_ints == 1) {
                //     stats.ints_distr[1] += 1;
                //     stats.codewords_distr[1] += 1;
                // } else if (decoded_ints == 2) {
                //     stats.ints_distr[2] += 2;
                //     stats.codewords_distr[2] += 1;
                // } else if (decoded_ints == 4) {
                //     stats.ints_distr[3] += 4;
                //     stats.codewords_distr[3] += 1;
                // } else if (decoded_ints == 8) {
                //     stats.ints_distr[4] += 8;
                //     stats.codewords_distr[4] += 1;
                // } else if (decoded_ints == 16 and index >
                // Dictionary::reserved - 1) {
                //     stats.ints_distr[5] += 16;
                //     stats.codewords_distr[5] += 1;
                // } else if (decoded_ints >= 16) {
                //     stats.ints_distr[0] += decoded_ints;
                //     stats.codewords_distr[0] += 1;
                // }
                // stats.dict_codewords++;
                // stats.decoded_ints_from_dict += decoded_ints;
                // stats.total_ints += decoded_ints;

            } else {
                // stats.ints_distr[6] += 1;
                // ++stats.total_ints;

                if (index == 1) {  // 4-byte exception
                    *out = *(reinterpret_cast<uint32_t const*>(++ptr));
                    ++ptr;
                    // stats.codewords_distr[6] += 3;

                    // if b = 8
                    // *out = *(reinterpret_cast<uint32_t const*>(++ptr));
                    // ptr += 3;
                } else {  // 2-byte exception
                    *out = *(++ptr);
                    // stats.codewords_distr[6] += 2;

                    // if b = 8
                    // *out = *(reinterpret_cast<uint16_t const*>(++ptr));
                    // ptr += 1;
                }
            }
            out += decoded_ints;
            i += decoded_ints;
        }

        return reinterpret_cast<uint8_t const*>(ptr);
        // if b = 8
        // return ptr;
    }
};

struct single_greedy_dint {
    // specialized encoding
    template <typename Builder>
    static void encode(Builder& builder, uint32_t const* in,
                       uint32_t /*universe*/, uint32_t n,
                       std::vector<uint8_t>& out) {
        uint32_t const* begin = in;
        uint32_t const* end = begin + n;

        while (begin < end) {
            uint32_t longest_run_size = 0;
            uint32_t run_size = std::min<uint64_t>(256, end - begin);
            uint32_t index = EXCEPTIONS;

            for (uint32_t const* ptr = begin; ptr != begin + run_size; ++ptr) {
                if (*ptr == 0) {
                    ++longest_run_size;
                } else {
                    break;
                }
            }

            if (longest_run_size >= 16) {
                uint32_t k = 256;
                while (longest_run_size < k and k > 16) {
                    ++index;
                    k /= 2;
                }
                write_index(index, out);
                begin += k;
            } else {
                for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
                    uint32_t sub_block_size = constants::target_sizes[s];
                    uint32_t len =
                        std::min<uint32_t>(sub_block_size, end - begin);
                    index = builder.lookup(begin, len);
                    if (index != Builder::invalid_index) {
                        write_index(index, out);
                        begin += len;
                        break;
                    }
                }

                if (index == Builder::invalid_index) {
                    uint32_t exception = *begin;
                    auto ptr = reinterpret_cast<uint8_t const*>(&exception);
                    if (exception < 65536) {
                        out.insert(out.end(), 0);
                        out.insert(out.end(), 0);  // comment if b = 8
                        out.insert(out.end(), ptr, ptr + 2);

                    } else {
                        out.insert(out.end(), 1);
                        out.insert(out.end(), 0);  // comment if b = 8
                        out.insert(out.end(), ptr, ptr + 4);
                    }

                    begin += 1;
                }
            }
        }
    }

    // generic decoding
    template <typename Dictionary>
    static uint8_t const* decode(Dictionary const& dict, uint8_t const* in,
                                 uint32_t* out, uint32_t /*universe*/, size_t n
                                 // , dint_statistics& stats
    ) {
        return single_dint::decode(dict, in, out, n
                                   // , stats
        );
    }

    static void write_index(uint32_t index, std::vector<uint8_t>& out) {
        auto ptr = reinterpret_cast<uint8_t const*>(&index);
        out.insert(out.end(), ptr, ptr + 2);  // b = 16
        // out.insert(out.end(), ptr, ptr + 1); // b = 8
    }
};

struct single_opt_dint {
    // specialized encoding
    template <typename Builder>
    static void encode(Builder& builder, uint32_t const* begin, uint64_t n,
                       std::vector<uint8_t>& out, int b) {
        std::vector<node> path;
        path.resize(n + 2);
        path[0] = {0, 1, 0};  // dummy node
        for (uint32_t i = 1; i < n + 1; ++i) {
            path[i] = {i - 1, 1, 3 * i};
        }

        for (uint32_t i = 0; i != n; ++i) {
            uint32_t longest_run_size = 0;
            uint32_t run_size = std::min<uint64_t>(256, n - i);
            uint32_t index = EXCEPTIONS;

            for (uint32_t j = i; j != i + run_size; ++j) {
                if (begin[j] == 0) {
                    ++longest_run_size;
                } else {
                    break;
                }
            }

            if (longest_run_size >= 16) {
                uint32_t k = 256;
                while (longest_run_size < k and k > 16) {
                    k /= 2;
                    ++index;
                }
                while (k >= 16) {
                    uint32_t c = path[i].cost + 1;
                    if (path[i + k].cost > c) {
                        path[i + k] = {i, index, c};
                    }

                    k /= 2;
                    ++index;
                }
            }

            for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
                uint32_t sub_block_size = constants::target_sizes[s];
                uint32_t len = std::min<uint32_t>(sub_block_size, n - i);
                index = builder.lookup(begin + i, len);
                if (index != Builder::invalid_index) {
                    uint32_t c = path[i].cost + 1;
                    if (path[i + len].cost > c) {
                        path[i + len] = {i, index, c};
                    }
                } else {
                    if (sub_block_size == 1) {  // exceptions
                        uint32_t exception = begin[i];
                        uint32_t c = path[i].cost + 2;  // small exception cost
                        index = 0;

                        if (exception > 65536 - 1) {
                            c += 1;  // large exception cost
                            index = 1;
                        }

                        if (path[i + 1].cost > c) {
                            path[i + 1] = {i, index, c};
                        }
                    }
                }
            }
        }

        std::vector<node> encoding;
        uint32_t i = n;
        while (i != 0) {
            uint32_t parent = path[i].parent;
            encoding.push_back(path[i]);
            i = parent;
        }

        std::reverse(encoding.begin(), encoding.end());
        encoding.emplace_back(n, 1, -1);  // final dummy node

        for (uint32_t i = 0, pos = 0; i < encoding.size() - 1; ++i) {
            uint32_t index = encoding[i].codeword;
            uint32_t len = encoding[i + 1].parent - encoding[i].parent;

            if (index > 1) {
                // ++builder.codewords;
                write_index(index, out, b);
            } else {
                assert(len == 1);
                uint32_t exception = begin[pos];
                auto ptr = reinterpret_cast<uint8_t const*>(&exception);

                if (index == 0) {
                    // ++builder.small_exceptions;
                    out.insert(out.end(), 0);
                    if (b == 16) {
                        out.insert(out.end(), 0);
                    }
                    out.insert(out.end(), ptr, ptr + 2);
                } else {
                    // ++builder.large_exceptions;
                    out.insert(out.end(), 1);
                    if (b == 16) {
                        out.insert(out.end(), 0);
                    }
                    out.insert(out.end(), ptr, ptr + 4);
                }
            }

            pos += len;
            assert(pos <= n);
        }

        encoding.clear();
    }

    template <typename Builder>
    static void encode(Builder& builder, uint32_t const* in,
                       uint32_t /*universe*/, uint32_t n,
                       std::vector<uint8_t>& out) {
        encode(builder, in, n, out, 16);
    }

    // generic decoding
    template <typename Dictionary>
    static uint8_t const* decode(Dictionary const& dict, uint8_t const* in,
                                 uint32_t* out, uint32_t /*universe*/, size_t n
                                 // , dint_statistics& stats
    ) {
        return single_dint::decode(dict, in, out, n
                                   // , stats
        );
    }

    static void write_index(uint32_t index, std::vector<uint8_t>& out, int b) {
        auto ptr = reinterpret_cast<uint8_t const*>(&index);
        assert(b == 8 or b == 16);
        out.insert(out.end(), ptr, ptr + b / 8);
    }
};

struct multi_opt_dint {
    // specialized encoding
    template <typename Builder>
    static void encode(Builder& builder, uint32_t dictionary_id,
                       uint32_t const* begin, uint64_t n,
                       std::vector<uint8_t>& out, int b) {
        std::vector<node> path;
        path.resize(n + 2);
        path[0] = {0, 1, 0};  // dummy node
        for (uint32_t i = 1; i < n + 1; ++i) {
            path[i] = {i - 1, 1, 3 * i};
        }

        for (uint32_t i = 0; i != n; ++i) {
            uint32_t longest_run_size = 0;
            uint32_t run_size = std::min<uint64_t>(256, n - i);
            uint32_t index = EXCEPTIONS;

            for (uint32_t j = i; j != i + run_size; ++j) {
                if (begin[j] == 0) {
                    ++longest_run_size;
                } else {
                    break;
                }
            }

            if (longest_run_size >= 16) {
                uint32_t k = 256;
                while (longest_run_size < k and k > 16) {
                    k /= 2;
                    ++index;
                }
                while (k >= 16) {
                    uint32_t c = path[i].cost + 1;
                    if (path[i + k].cost > c) {
                        path[i + k] = {i, index, c};
                    }

                    k /= 2;
                    ++index;
                }
            }

            for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
                uint32_t sub_block_size = constants::target_sizes[s];
                uint32_t len = std::min<uint32_t>(sub_block_size, n - i);
                index = builder.lookup(dictionary_id, begin + i, len, b);
                if (index != Builder::invalid_index) {
                    uint32_t c = path[i].cost + 1;
                    if (path[i + len].cost > c) {
                        path[i + len] = {i, index, c};
                    }
                } else {
                    if (sub_block_size == 1) {  // exceptions
                        uint32_t exception = begin[i];
                        uint32_t c = path[i].cost + 2;  // small exception cost
                        index = 0;

                        if (exception > 65536 - 1) {
                            c += 1;  // large exception cost
                            index = 1;
                        }

                        if (path[i + 1].cost > c) {
                            path[i + 1] = {i, index, c};
                        }
                    }
                }
            }
        }

        std::vector<node> encoding;
        uint32_t i = n;
        while (i != 0) {
            uint32_t parent = path[i].parent;
            encoding.push_back(path[i]);
            i = parent;
        }

        std::reverse(encoding.begin(), encoding.end());
        encoding.emplace_back(n, 1, -1);  // final dummy node

        for (uint32_t i = 0, pos = 0; i < encoding.size() - 1; ++i) {
            uint32_t index = encoding[i].codeword;
            uint32_t len = encoding[i + 1].parent - encoding[i].parent;

            if (index > 1) {
                // ++builder.codewords;
                single_opt_dint::write_index(index, out, b);
            } else {
                assert(len == 1);
                uint32_t exception = begin[pos];
                auto ptr = reinterpret_cast<uint8_t const*>(&exception);

                if (index == 0) {
                    // ++builder.small_exceptions;
                    out.insert(out.end(), 0);
                    if (b == 16) {
                        out.insert(out.end(), 0);
                    }
                    out.insert(out.end(), ptr, ptr + 2);
                } else {
                    // ++builder.large_exceptions;
                    out.insert(out.end(), 1);
                    if (b == 16) {
                        out.insert(out.end(), 0);
                    }
                    out.insert(out.end(), ptr, ptr + 4);
                }
            }

            pos += len;
            assert(pos <= n);
        }

        encoding.clear();
    }

    template <typename Builder>
    static void encode(Builder& builder, uint32_t const* in,
                       uint32_t /*universe*/, uint32_t n,
                       std::vector<uint8_t>& out) {
        uint32_t const* begin = in;
        uint64_t num_blocks =
            succinct::util::ceil_div(n, constants::block_size);
        uint64_t tail = n - (n / constants::block_size * constants::block_size);

        for (uint64_t b = 0; b != num_blocks; ++b) {
            uint64_t size = constants::block_size;
            if (b == num_blocks - 1 and tail != 0) {
                size = tail;
            }

            std::vector<std::vector<uint8_t>> encoded;

            // option 1: choose the best dictionary (exhaustive search)
            size_t best_size = size_t(-1);
            uint32_t selector_code = 0;
            encoded.resize(2 * constants::num_selectors);
            for (uint32_t s = 0; s != constants::num_selectors; ++s) {
                encode(builder, s, begin, size, encoded[s], 16);
                encode(builder, s, begin, size,
                       encoded[s + constants::num_selectors], 8);

                size_t smallest_size = encoded[s].size();
                uint32_t sc = s;
                if (encoded[s + constants::num_selectors].size() <=
                    smallest_size) {
                    smallest_size =
                        encoded[s + constants::num_selectors].size();
                    sc += constants::num_selectors;
                }

                if (smallest_size < best_size) {
                    best_size = smallest_size;
                    selector_code = sc;
                }
            }

            out.push_back(selector_code);
            out.insert(out.end(), encoded[selector_code].begin(),
                       encoded[selector_code].end());
            for (auto& e : encoded) {
                e.clear();
            }

            // // option 2: select the dictionary based on the context
            // selector sct;
            // encoded.resize(2);
            // uint32_t selector_code = sct.get(begin, size);
            // encode(builder, selector_code, begin, size, encoded[0], 16);
            // encode(builder, selector_code, begin, size, encoded[1],  8);
            // size_t smallest_size = encoded[0].size();
            // if (encoded[1].size() <= smallest_size) {
            //     selector_code += constants::num_selectors;
            // }

            // // control byte
            // out.push_back(selector_code);
            // out.insert(out.end(), encoded[selector_code >=
            // constants::num_selectors].begin(),
            //                       encoded[selector_code >=
            //                       constants::num_selectors].end());

            begin += size;
        }
    }

    // specialized decoding
    template <typename Dictionary>
    static uint8_t const* decode(Dictionary const& dict, uint8_t const* in,
                                 uint32_t* out, uint32_t /*universe*/, size_t n
                                 // , dint_statistics& stats
    ) {
        uint64_t num_blocks =
            succinct::util::ceil_div(n, constants::block_size);
        size_t tail = n - (n / constants::block_size * constants::block_size);
        // uint64_t sum = 0;
        for (uint64_t b = 0; b != num_blocks; ++b) {
            size_t size = constants::block_size;
            if (b == num_blocks - 1 and tail != 0) {
                size = tail;
            }

            uint8_t selector_code = *in;
            if (selector_code < constants::num_selectors) {
                uint16_t const* ptr = reinterpret_cast<uint16_t const*>(in + 1);
                for (size_t i = 0; i != size; ++ptr) {
                    uint32_t index = *ptr;
                    uint32_t decoded_ints = 1;

                    // ++stats.occs[index];

                    if (DS2I_LIKELY(index > EXCEPTIONS - 1)) {
                        decoded_ints = dict.copy(selector_code, index, out);

                        // stats.total_ints += decoded_ints;

                    } else {
                        // ++stats.total_ints;

                        if (index == 1) {  // 4-byte exception
                            uint32_t exception =
                                *(reinterpret_cast<uint32_t const*>(++ptr));
                            *out = exception;
                            ++ptr;

                            // stats.eat(exception);

                        } else {  // 2-byte exception
                            uint32_t exception = *(++ptr);
                            *out = exception;

                            // stats.eat(exception);
                        }
                    }
                    out += decoded_ints;
                    i += decoded_ints;
                }
                in = reinterpret_cast<uint8_t const*>(ptr);
            } else {
                selector_code -= constants::num_selectors;
                uint8_t const* ptr = in + 1;
                for (size_t i = 0; i != size; ++ptr) {
                    uint32_t index = *ptr;
                    uint32_t decoded_ints = 1;

                    // ++stats.occs[index];

                    if (DS2I_LIKELY(index > EXCEPTIONS - 1)) {
                        decoded_ints = dict.copy(selector_code, index, out);

                        // stats.total_ints += decoded_ints;

                    } else {
                        // ++stats.total_ints;

                        if (index == 1) {  // 4-byte exception
                            uint32_t exception =
                                *(reinterpret_cast<uint32_t const*>(++ptr));
                            *out = exception;
                            ptr += 3;

                            // stats.eat(exception);

                        } else {  // 2-byte exception
                            uint32_t exception =
                                *(reinterpret_cast<uint16_t const*>(++ptr));
                            *out = exception;
                            ptr += 1;

                            // stats.eat(exception);
                        }
                    }
                    out += decoded_ints;
                    i += decoded_ints;
                }
                in = ptr;
            }

            // sum += size;
            // assert(sum <= n);
        }

        // assert(sum == n);

        return in;
    }
};
}  // namespace ds2i