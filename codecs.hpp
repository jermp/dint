#pragma once

#include "FastPFor/headers/optpfor.h"
#include "FastPFor/headers/variablebyte.h"
#include "streamvbyte/include/streamvbyte.h"
#include "MaskedVByte/include/varintencode.h"
#include "MaskedVByte/include/varintdecode.h"
#include "VarIntG8IU.h"
#include "varintgb.h"
#include "interpolative_coding.hpp"
#include "qmx.hpp"

#include "succinct/util.hpp"
#include "util.hpp"
#include "dictionary.hpp"

#include <unordered_map>

namespace ds2i {

    typedef dictionary<constants::num_entries,
                       constants::max_entry_size> dictionary_type;

    // workaround: VariableByte::decodeArray needs the buffer size, while we
    // only know the number of values. It also pads to 32 bits. We need to
    // rewrite
    struct TightVariableByte {

      template <uint32_t i> static uint8_t extract7bits(const uint32_t val) {
        return static_cast<uint8_t>((val >> (7 * i)) & ((1U << 7) - 1));
      }

      template <uint32_t i>
      static uint8_t extract7bitsmaskless(const uint32_t val) {
        return static_cast<uint8_t>((val >> (7 * i)));
      }

      static void encode(const uint32_t *in, const size_t length, uint8_t *out,
                         size_t &nvalue) {
        uint8_t *bout = out;
        for (size_t k = 0; k < length; ++k) {
          const uint32_t val(in[k]);
          /**
           * Code below could be shorter. Whether it could be faster
           * depends on your compiler and machine.
           */
          if (val < (1U << 7)) {
            *bout = static_cast<uint8_t>(val | (1U << 7));
            ++bout;
          } else if (val < (1U << 14)) {
            *bout = extract7bits<0>(val);
            ++bout;
            *bout = extract7bitsmaskless<1>(val) | (1U << 7);
            ++bout;
          } else if (val < (1U << 21)) {
            *bout = extract7bits<0>(val);
            ++bout;
            *bout = extract7bits<1>(val);
            ++bout;
            *bout = extract7bitsmaskless<2>(val) | (1U << 7);
            ++bout;
          } else if (val < (1U << 28)) {
            *bout = extract7bits<0>(val);
            ++bout;
            *bout = extract7bits<1>(val);
            ++bout;
            *bout = extract7bits<2>(val);
            ++bout;
            *bout = extract7bitsmaskless<3>(val) | (1U << 7);
            ++bout;
          } else {
            *bout = extract7bits<0>(val);
            ++bout;
            *bout = extract7bits<1>(val);
            ++bout;
            *bout = extract7bits<2>(val);
            ++bout;
            *bout = extract7bits<3>(val);
            ++bout;
            *bout = extract7bitsmaskless<4>(val) | (1U << 7);
            ++bout;
          }
        }
        nvalue = bout - out;
      }

      static void encode_single(uint32_t val, std::vector<uint8_t> &out) {
        uint8_t buf[5];
        size_t nvalue;
        encode(&val, 1, buf, nvalue);
        out.insert(out.end(), buf, buf + nvalue);
      }

      static uint8_t const *decode(const uint8_t *in, uint32_t *out, size_t n) {
        const uint8_t *inbyte = in;
        for (size_t i = 0; i < n; ++i) {
          unsigned int shift = 0;
          for (uint32_t v = 0;; shift += 7) {
            uint8_t c = *inbyte++;
            v += ((c & 127) << shift);
            if ((c & 128)) {
              *out++ = v;
              break;
            }
          }
        }
        return inbyte;
      }
    };

    struct header {
        static void write(uint32_t n, uint32_t universe,
                          std::vector<uint8_t>& out)
        {
            TightVariableByte::encode_single(n, out);
            TightVariableByte::encode_single(universe, out);
        }

        static uint8_t const* read(uint8_t const* in,
                                   uint32_t* n, uint32_t* universe)
        {
            uint8_t const* inbuf = in;
            inbuf = TightVariableByte::decode(inbuf, n, 1);
            inbuf = TightVariableByte::decode(inbuf, universe, 1);
            return inbuf;
        }
    };

    struct interpolative {

        static void encode(uint32_t const* in,
                           uint32_t universe, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* builder = nullptr)
        {
            (void) builder;
            std::vector<uint32_t> inbuf(n);
            std::vector<uint32_t> outbuf;
            inbuf[0] = *in;
            for (size_t i = 1; i < n; ++i) {
                inbuf[i] = inbuf[i - 1] + in[i];
            }

            bit_writer bw(outbuf);
            bw.write_interpolative(inbuf.data(), n - 1, 0, universe);
            uint8_t const *bufptr = (uint8_t const *)outbuf.data();
            out.insert(out.end(), bufptr,
                       bufptr + succinct::util::ceil_div(bw.size(), 8));
        }

        static uint8_t const* DS2I_NOINLINE decode(uint8_t const* in,
                                                   uint32_t* out,
                                                   uint32_t universe, uint32_t n,
                                                   dictionary_type const* dict = nullptr)
        {
            (void) dict;
            uint8_t const* inbuf = in;
            out[n - 1] = universe;
            size_t read_interpolative = 0;
            if (n > 1) {
                bit_reader br((uint32_t const *)inbuf);
                br.read_interpolative(out, n - 1, 0, universe);
                for (size_t i = n - 1; i > 0; --i) {
                    out[i] -= out[i - 1];
                }
                read_interpolative = succinct::util::ceil_div(br.position(), 8);
            }

            return inbuf + read_interpolative;
        }
    };

    struct optpfor {

        struct codec_type : FastPFor::OPTPFor<4, FastPFor::Simple16<false>> {

            uint8_t const* force_b;

            uint32_t findBestB(const uint32_t *in, uint32_t len) {
                // trick to force the choice of b from a parameter
                if (force_b) {
                    return *force_b;
                }

                // this is mostly a cut&paste from FastPFor, but we stop the
                // optimization early as the b to test becomes larger than maxb
                uint32_t b = 0;
                uint32_t bsize = std::numeric_limits<uint32_t>::max();
                const uint32_t mb = FastPFor::maxbits(in, in + len);
                uint32_t i = 0;
                while (mb > 28 + possLogs[i])
                    ++i; // some schemes such as Simple16 don't code numbers greater than 28

                for (; i < possLogs.size(); i++) {
                    if (possLogs[i] > mb && possLogs[i] >= mb)
                        break;
                    const uint32_t csize = tryB(possLogs[i], in, len);

                    if (csize <= bsize) {
                        b = possLogs[i];
                        bsize = csize;
                    }
                }

                return b;
            }

            const uint32_t* decodeArray(const uint32_t *in, const size_t, uint32_t *out, size_t & nvalue)
            {
                const uint32_t * const initout(out);
                const uint32_t numBlocks = *in++;
                size_t out_len = 0;
                for (uint32_t i = 0; i < numBlocks; i++) {
                    in = decodeBlock(in, out, out_len);
                    out += out_len;
                }
                nvalue = out - initout;
                assert(nvalue == numBlocks * BlockSize);
                return in;
            }
        };

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/,
                           uint8_t const* b = nullptr) // if non-null forces b
        {
            if (n == 1) {
                TightVariableByte::encode_single(in[0], out);
                return;
            }

            thread_local codec_type optpfor_codec;
            thread_local std::vector<uint8_t> buf(2 * 4 * n);
            size_t out_len = buf.size();
            optpfor_codec.force_b = b;
            if (n % 128) {
                n = (n / 128 + 1) * 128;
            }
            optpfor_codec.encodeArray(in, n, reinterpret_cast<uint32_t *>(buf.data()), out_len); // encodeBlock
            out_len *= 4;
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* DS2I_NOINLINE decode(uint8_t const* in,
                                                   uint32_t* out,
                                                   uint64_t /*universe*/, uint32_t n,
                                                   dictionary_type const* /*dict*/)
        {
            if (DS2I_UNLIKELY(n == 1)) {
                return TightVariableByte::decode(in, out, 1);
            }

            thread_local codec_type optpfor_codec; // pfor decoding is *not* thread-safe
            size_t out_len = n;
            uint8_t const *ret;
            if (n % 128) {
                n = (n / 128 + 1) * 128;
            }
            ret = reinterpret_cast<uint8_t const *>(optpfor_codec.decodeArray(
                        reinterpret_cast<uint32_t const *>(in), n, out, out_len
                    ));
            assert(out_len == n);
            return ret;
        }
    };

    struct varintg8iu {

        struct codec_type : VarIntG8IU {

            // rewritten version of decodeBlock optimized for when the output
            // size is known rather than the input
            // the buffers pointed by src and dst must be respectively at least
            // 9 and 8 elements large
            uint32_t decodeBlock(uint8_t const *&src, uint32_t *dst) const
            {
                uint8_t desc = *src;
                src += 1;
                const __m128i data =
                    _mm_lddqu_si128(reinterpret_cast<__m128i const *>(src));
                src += 8;
                const __m128i result = _mm_shuffle_epi8(data, vecmask[desc][0]);
                _mm_storeu_si128(reinterpret_cast<__m128i *>(dst), result);
                int readSize = maskOutputSize[desc];

                if (readSize > 4) {
                    const __m128i result2 = _mm_shuffle_epi8(
                        data, vecmask[desc][1]);
                    _mm_storeu_si128(
                        reinterpret_cast<__m128i *>(dst + 4),
                        result2);
                }

                return readSize;
            }
        };

        static void encode(uint32_t const* in,
                           uint32_t universe, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            // if (n == 1) {
            //     TightVariableByte::encode_single(in[0], out);
            //     return;
            // }

            if (n < 8) {
                interpolative::encode(in, universe, n, out);
                return;
            }

            thread_local codec_type varint_codec;
            thread_local std::vector<uint8_t> buf(2 * 4 * n);
            size_t out_len = buf.size();
            const uint32_t *src = in;
            unsigned char *dst = buf.data();
            size_t srclen = n * 4;
            size_t dstlen = out_len;
            out_len = 0;
            while (srclen > 0 && dstlen >= 9) {
                out_len += varint_codec.encodeBlock(src, srclen, dst, dstlen);
            }
            assert(srclen == 0);
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        // we only allow varint to be inlined (others have DS2I_NOILINE)
        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t universe, uint32_t n,
                                     dictionary_type const* /*dict*/)
        {
            // if (DS2I_UNLIKELY(n == 1)) {
            //     return TightVariableByte::decode(in, out, 1);
            // }

            if (DS2I_UNLIKELY(n < 8)) {
                return interpolative::decode(in, out, universe, n);
            }

            static codec_type varint_codec; // decodeBlock is thread-safe

            size_t out_len = 0;
            uint8_t const *src = in;
            uint32_t *dst = out;
            while (out_len <= (n - 8)) {
                out_len += varint_codec.decodeBlock(src, dst + out_len);
            }

            // decodeBlock can overshoot, so we decode the last blocks in a
            // local buffer
            while (out_len < n) {
                uint32_t buf[8];
                size_t read = varint_codec.decodeBlock(src, buf);
                size_t needed = std::min(read, n - out_len);
                memcpy(dst + out_len, buf, needed * 4);
                out_len += needed;
            }
            assert(out_len == n);
            return src;
        }
    };

    struct qmx {

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            if (n == 1) {
                TightVariableByte::encode_single(in[0], out);
                return;
            }

            static const uint64_t overflow = 512;
            thread_local std::vector<uint8_t> buf( (overflow * 4) + 2 * 4 * n);
            QMX::codec qmx_codec(n);
            size_t compressed_bytes = qmx_codec.encode(buf.data(), in);
            TightVariableByte::encode_single(compressed_bytes, out);
            out.insert(out.end(), buf.data(), buf.data() + compressed_bytes);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, uint32_t n,
                                     dictionary_type const* /*dict*/)
        {
            if (DS2I_UNLIKELY(n == 1)) {
                return TightVariableByte::decode(in, out, 1);
            }

            uint32_t compressed_bytes = 0;
            in = TightVariableByte::decode(in, &compressed_bytes, 1);
            QMX::codec qmx_codec(n);
            qmx_codec.decode(out, in, compressed_bytes);
            return in + compressed_bytes;
        }
    };

    struct vbyte {

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            std::vector<uint8_t> buf(2 * 4 * n);
            size_t out_len = buf.size();
            TightVariableByte::encode(in, n, buf.data(), out_len);
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n,
                                     dictionary_type const* /*dict*/)
        {
            return TightVariableByte::decode(in, out, n);
        }
    };

    struct u32 {

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            size_t srclen = n * sizeof(uint32_t);
            const uint8_t* src = (const uint8_t*)in;
            out.insert(out.end(), src, src + srclen);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n,
                                     dictionary_type const* /*dict*/)
        {
            const uint8_t* src = (const uint8_t*)in;
            uint8_t* dst = (uint8_t*)out;
            size_t n4 = n * 4;
            for (size_t i = 0; i < n4; i++) {
                *dst++ = *src++;
            }
            return src;
        }
    };

    struct simple16 {

        using codec_type = FastPFor::Simple16<false>;

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            thread_local codec_type simple16_codec;
            std::vector<uint8_t> buf(2 * 8 * n);
            size_t out_len = buf.size();
            simple16_codec.encodeArray(in, n, reinterpret_cast<uint32_t*>(buf.data()), out_len);
            out_len *= 4;
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n,
                                     dictionary_type const* /*dict*/)
        {
            uint8_t const* ret;
            thread_local codec_type simple16_codec;
            ret = reinterpret_cast<uint8_t const*>(simple16_codec.decodeArray(reinterpret_cast<uint32_t const*>(in), 1,
                out, n));
            return ret;
        }
    };

    struct streamvbyte {

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            uint32_t *src = const_cast<uint32_t*>(in);
            std::vector<uint8_t> buf(streamvbyte_max_compressedbytes(n));
            size_t out_len = streamvbyte_encode(src, n, buf.data());
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n,
                                     dictionary_type const* /*dict*/)
        {
            auto read = streamvbyte_decode(in, out, n);
            return in + read;
        }
    };

    struct maskedvbyte {

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            uint32_t* src = const_cast<uint32_t*>(in);
            std::vector<uint8_t> buf(2 * n * sizeof(uint32_t));
            size_t out_len = vbyte_encode(src, n, buf.data());
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n,
                                     dictionary_type const* /*dict*/)
        {
            auto read = masked_vbyte_decode(in, out, n);
            return in + read;
        }
    };

    struct varintgb {

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out,
                           dictionary_type::builder const* /*builder*/)
        {
            thread_local VarIntGB<false> varintgb_codec;
            thread_local std::vector<uint8_t> buf(2 * n * sizeof(uint32_t));
            size_t out_len = varintgb_codec.encodeArray(in, n, buf.data());
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n,
                                     dictionary_type const* /*dict*/)
        {
            thread_local VarIntGB<false> varintgb_codec;
            auto read = varintgb_codec.decodeArray(in, n, out);
            return read + in;
        }
    };

    struct dint_statistics
    {
        dint_statistics(uint32_t /*num_entries*/)
            : ints_distr(7, 0)
            , codewords_distr(7, 0)
        {}

        std::vector<uint64_t> ints_distr;       // 0:runs; 1:1; 2:2; 3:4; 4:8; 5:16; 6:exceptions
        std::vector<uint64_t> codewords_distr;  // 0:runs; 1:1; 2:2; 3:4; 4:8; 5:16; 6:exceptions

        uint64_t decoded_ints_from_dict = 0;
        uint64_t dict_codewords = 0;
    };

    struct node {
        node()
        {}

        node(uint32_t p, uint32_t w, uint32_t c)
            : parent(p), codeword(w), cost(c)
        {}

        uint32_t parent;
        uint32_t codeword;
        uint32_t cost;
    };

    struct dint {

        // NOTE: greedy parsing
        static void
        encode(uint32_t const* in,
               uint32_t /*universe*/, uint32_t n,
               std::vector<uint8_t>& out,
               dictionary_type::builder* builder)
        {
            uint32_t const* begin = in;
            uint32_t const* end = begin + n;
            // uint64_t written_ints = 0;

            while (begin < end)
            {
                uint32_t longest_run_size = 0;
                uint32_t run_size = std::min<uint64_t>(256, end - begin);
                uint32_t index = EXCEPTIONS;
                uint32_t len = 1;

                for (uint32_t const* ptr  = begin;
                                     ptr != begin + run_size;
                                   ++ptr)
                {
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
                    len = k;
                } else {
                    for (uint32_t s = 0; s < constants::num_target_sizes; ++s)
                    {
                        uint32_t sub_block_size = constants::target_sizes[s];
                        len = std::min<uint32_t>(sub_block_size, end - begin);
                        index = builder->lookup(begin, len);
                        if (index != dictionary_type::invalid_index) {
                            write_index(index, out);
                            break;
                        }
                    }

                    if (index == dictionary_type::invalid_index)
                    {
                        len = 1;
                        uint32_t exception = *begin;
                        auto ptr = reinterpret_cast<uint8_t const*>(&exception);

                        // USED WITH EXCEPTIONS = 1
                        // out.insert(out.end(), 0);
                        // out.insert(out.end(), 0); // comment if b = 8
                        // out.insert(out.end(), ptr, ptr + 4);

                        if (exception < 65536) {
                            out.insert(out.end(), 0);
                            out.insert(out.end(), 0); // comment if b = 8
                            out.insert(out.end(), ptr, ptr + 2);
                        } else {
                            out.insert(out.end(), 1);
                            out.insert(out.end(), 0); // comment if b = 8
                            out.insert(out.end(), ptr, ptr + 4);
                        }
                    }
                }

                begin += len;
                // written_ints += len;
            }

            // return written_ints;
        }

        // NOTE: optimal parsing
        // static void encode(uint32_t const* in,
        //                    uint32_t /*universe*/, uint32_t n,
        //                    std::vector<uint8_t>& out,
        //                    dictionary_type::builder* builder)
        // {
        //     // NOTE: everything at the beginning is a large exception
        //     // costs are in shorts! (1 short = 16 bits)
        //     std::vector<node> path(n + 2);
        //     path[0] = {0, 1, 0}; // dummy node
        //     for (uint32_t i = 1; i < n + 1; ++i) {
        //         path[i] = {i - 1, 1, 3 * i};
        //     }

        //     // logger() << "finding shortest path" << std::endl;

        //     for (uint32_t i = 0; i < n; ++i)
        //     {
        //         // std::cout << "current node " << i << ":\n";
        //         // std::cout << "parent " << path[i].parent << "\n";
        //         // std::cout << "codeword " << path[i].codeword << "\n";
        //         // std::cout << "cost " << path[i].cost << "\n\n";

        //         uint32_t longest_run_size = 0;
        //         uint32_t run_size = std::min<uint64_t>(256, n - i);
        //         uint32_t index = EXCEPTIONS;

        //         for (uint32_t j = i; j != i + run_size; ++j) {
        //             if (in[j] == 0) {
        //                 ++longest_run_size;
        //             } else {
        //                 break;
        //             }
        //         }

        //         // std::cout << "longest_run_size " << longest_run_size << std::endl;

        //         if (longest_run_size >= 16) {
        //             uint32_t k = 256;
        //             while (longest_run_size < k and k > 16) {
        //                 k /= 2;
        //                 ++index;
        //             }
        //             while (k >= 16) {
        //                 // std::cout << "k " << k << "; index " << index << std::endl;

        //                 uint32_t c = path[i].cost + 1;
        //                 if (path[i + k].cost > c) {
        //                     path[i + k] = {i, index, c};
        //                 }

        //                 k /= 2;
        //                 ++index;
        //             }

        //             // std::cout << std::endl;
        //         }

        //         for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
        //             uint32_t sub_block_size = constants::target_sizes[s];
        //             uint32_t len = std::min<uint32_t>(sub_block_size, n - i);
        //             index = builder->lookup(in + i, len);
        //             if (index != dictionary_type::invalid_index) {
        //                 uint32_t c = path[i].cost + 1;
        //                 if (path[i + len].cost > c) {
        //                     path[i + len] = {i, index, c};
        //                 }
        //             } else {
        //                 if (sub_block_size == 1) { // exceptions
        //                     uint32_t exception = in[i];
        //                     uint32_t c = path[i].cost + 2; // small exception cost
        //                     index = 0;

        //                     if (exception > 65536) {
        //                         c += 1; // large exception cost
        //                         index = 1;
        //                     }

        //                     if (path[i + 1].cost > c) {
        //                         path[i + 1] = {i, index, c};
        //                     }
        //                 }
        //             }
        //         }
        //     }

        //     // std::cout << "min_cost " << path.back().cost << std::endl;

        //     std::vector<node> encoding;
        //     uint32_t i = n;
        //     while (i != 0) {
        //         uint32_t parent = path[i].parent;
        //         encoding.push_back(path[i]);
        //         i = parent;
        //     }

        //     // std::cout << std::endl;
        //     std::reverse(encoding.begin(), encoding.end());
        //     encoding.emplace_back(n, 1, -1); // final dummy node

        //     // logger() << "encoding" << std::endl;
        //     uint32_t pos = 0;
        //     uint32_t cost = 0;
        //     for (uint32_t i = 0; i < encoding.size() - 1; ++i) {
        //         uint32_t index = encoding[i].codeword;
        //         uint32_t len = encoding[i + 1].parent - encoding[i].parent;

        //         assert(len == builder->size(index));

        //         cost += 1;

        //         if (index > 1) {
        //             write_index(index, out);
        //         } else {

        //             assert(len == 1);
        //             uint32_t exception = in[pos];
        //             auto ptr = reinterpret_cast<uint8_t const*>(&exception);
        //             cost += 1;

        //             if (index == 0) {
        //                 out.insert(out.end(), 0);
        //                 out.insert(out.end(), 0); // comment if b = 8
        //                 out.insert(out.end(), ptr, ptr + 2);
        //             } else {
        //                 cost += 1;
        //                 out.insert(out.end(), 1);
        //                 out.insert(out.end(), 0); // comment if b = 8
        //                 out.insert(out.end(), ptr, ptr + 4);
        //             }
        //         }

        //         // std::cout << "cost " << cost << "; pos " << pos << "; len " << len << "; codeword: " << index << std::endl;

        //         pos += len;
        //     }

        //     assert(pos == n);

        //     // std::cout << "pos " << pos << "/" << block_size << std::endl;
        //     // std::cout << "cost = " << cost << std::endl;
        // }

        // static uint8_t const* decode(uint8_t const* in,
        //                              uint32_t* out,
        //                              uint32_t /*universe*/, size_t n,
        //                              dictionary_type const* dict
        //                              , dint_statistics& stats
        //                              )
        // {
        //     uint16_t const* ptr = reinterpret_cast<uint16_t const*>(in);
        //     // uint8_t const* ptr = in;

        //     for (size_t i = 0; i != n; ++ptr) {
        //         uint32_t index = *ptr;
        //         uint32_t decoded_ints = 1;

        //         if (DS2I_LIKELY(index > dictionary_type::reserved - 1))
        //         {
        //             // NOTE1: on Gov2 and decoding the docIDs,
        //             // this IF if executed for 90.35% of the codewords

        //             // NOTE2: on Gov2 and decoding the docIDSs,
        //             // if we count the groups of 8 consecutive codewords that requires
        //             // only a copy, i.e., this IF, we cover 82.35% of the codewords

        //             decoded_ints = dict->copy(index, out);

        //             // if (decoded_ints == 1) {
        //             //     stats.ints_distr[1] += 1;
        //             //     stats.codewords_distr[1] += 1;
        //             // } else if (decoded_ints == 16) {
        //             //     stats.ints_distr[5] += 16;
        //             //     stats.codewords_distr[5] += 1;
        //             // }

        //             // std::cout << "0\n";

        //             // stats.dict_codewords++;
        //             // stats.decoded_ints_from_dict += decoded_ints;

        //         } else {

        //             // NOTE: on Gov2 and decoding the docIDs,
        //             // this IF if executed for 9.64% of the codewords

        //             static const uint32_t run_lengths[] = {0, 1, // exceptions
        //                                                    256, 128, 64, 32, 16};
        //             decoded_ints = run_lengths[index];

        //             if (DS2I_UNLIKELY(decoded_ints == 1)) {
        //                 *out = *(reinterpret_cast<uint32_t const*>(++ptr));
        //                 ++ptr;

        //                 // stats.ints_distr[6] += 1;
        //                 // stats.codewords_distr[6] += 3;

        //                 // needed when b = 8
        //                 // ptr += 2;
        //             }

        //             if (DS2I_UNLIKELY(decoded_ints == 0)) { // 2-byte exception
        //                 // *out = *(reinterpret_cast<uint16_t const*>(++ptr)); // when b = 8
        //                 *out = *(++ptr);
        //                 decoded_ints = 1;

        //                 // needed when b = 8
        //                 // ptr += 1;
        //             }
        //         }

        //         out += decoded_ints;
        //         i += decoded_ints;

        //         // if (decoded_ints >= 16) {
        //         //     stats.ints_distr[0] += decoded_ints;
        //         //     stats.codewords_distr[0] += 1;
        //         // } else if (decoded_ints == 8) {
        //         //     stats.ints_distr[4] += 8;
        //         //     stats.codewords_distr[4] += 1;
        //         // } else if (decoded_ints == 4) {
        //         //     stats.ints_distr[3] += 4;
        //         //     stats.codewords_distr[3] += 1;
        //         // } else if (decoded_ints == 2) {
        //         //     stats.ints_distr[2] += 2;
        //         //     stats.codewords_distr[2] += 1;
        //         // }
        //     }

        //     return reinterpret_cast<uint8_t const*>(ptr);
        // }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n,
                                     dictionary_type const* dict
                                     , dint_statistics& stats
                                     )
        {
            // uint16_t const* ptr = reinterpret_cast<uint16_t const*>(in);
            // for (size_t i = 0; i != n; ++ptr) {
            //     uint32_t index = *ptr;
            //     uint32_t decoded_ints = 1;
            //     if (DS2I_LIKELY(index > EXCEPTIONS - 1)) {
            //         decoded_ints = dict->copy(index, out);
            //     } else {
            //         if (index == 1) { // 4-byte exception
            //             *out = *(reinterpret_cast<uint32_t const*>(++ptr));
            //             ++ptr;
            //         } else { // 2-byte exception
            //             *out = *(++ptr);
            //         }
            //     }
            //     out += decoded_ints;
            //     i += decoded_ints;
            // }
            // return reinterpret_cast<uint8_t const*>(ptr);


            uint16_t const* ptr = reinterpret_cast<uint16_t const*>(in);

            uint32_t j = 0;
            uint32_t sum = 0;

            static const int len = 4;
            static uint32_t S[] = {0, 0, 0, 0
                                 // , 0, 0, 0, 0
                            };

            static uint32_t I[] = {0, 0, 0, 0
                                 // , 0, 0, 0, 0
                            };

            // static uint32_t tmp[][16] = {
            //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            //   , {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            //   , {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            //   , {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            //   // , {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            //   // , {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            //   // , {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            //   // , {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            // };

            // static uint32_t tmp[] = {
            //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            //   , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            //   , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            //   , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            //   // , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            //   // , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            //   // , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            //   // , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            // };

            // std::cout << "n " << n << std::endl;

            for (size_t i = 0; i < n;)
            {
                // std::cout << i << std::endl;

                while (j != len) {
                    uint32_t index = *ptr;
                    if (DS2I_LIKELY(index > EXCEPTIONS - 1)) {

                        uint32_t size = dict->size(index);
                        sum += size;
                        S[j] = sum;
                        I[j] = index;
                        ++j;

                    } else {
                        if (index == 1) { // 4-byte exception
                            *(out + sum) = *(reinterpret_cast<uint32_t const*>(++ptr));
                            ++ptr;
                        } else { // 2-byte exception
                            *(out + sum) = *(++ptr);
                        }
                        ++sum;
                    }

                    ++ptr;

                    if (i + sum >= n) {
                        break;
                    }
                }

                dict->copy(I[0], out +   0 );
                dict->copy(I[1], out + S[0]);
                dict->copy(I[2], out + S[1]);
                dict->copy(I[3], out + S[2]);
                // dict->copy(I[4], out + S[3]);
                // dict->copy(I[5], out + S[4]);
                // dict->copy(I[6], out + S[5]);
                // dict->copy(I[7], out + S[6]);

                // dict->copy(I[0], tmp[0]);
                // dict->copy(I[1], tmp[1]);
                // dict->copy(I[2], tmp[2]);
                // dict->copy(I[3], tmp[3]);
                // dict->copy(I[4], tmp[4]);
                // dict->copy(I[5], tmp[5]);
                // dict->copy(I[6], tmp[6]);
                // dict->copy(I[7], tmp[7]);

                // dict->copy(I[0], tmp +   0 );
                // dict->copy(I[1], tmp + S[0]);
                // dict->copy(I[2], tmp + S[1]);
                // dict->copy(I[3], tmp + S[2]);
                // memcpy(out, tmp, len * 16 * sizeof(uint32_t));

                out += sum;
                i += sum;

                j = 0;
                sum = 0;
            }

            return reinterpret_cast<uint8_t const*>(ptr);
        }

    private:
        static void write_index(uint32_t index, std::vector<uint8_t>& out) {
            auto ptr = reinterpret_cast<uint8_t const*>(&index);
            out.insert(out.end(), ptr, ptr + 2); // b = 16
            // out.insert(out.end(), ptr, ptr + 1); // b = 8
        }
    };

    // struct dint { // NOTE: no exception code

    //     static uint64_t encode(uint32_t const* in,
    //                            uint32_t /*universe*/, uint32_t n,
    //                            std::vector<uint8_t>& out,
    //                            dictionary_type::builder* builder
    //                           )
    //     {
    //         uint32_t const* begin = in;
    //         uint32_t const* end = begin + n;

    //         uint64_t written = 0;
    //         while (begin < end)
    //         {
    //             uint32_t index = EXCEPTIONS;
    //             for (uint32_t s = 0; s < constants::num_target_sizes; ++s)
    //             {
    //                 uint32_t sub_block_size = constants::target_sizes[s];
    //                 uint32_t len = std::min<uint32_t>(sub_block_size, end - begin);
    //                 index = builder->lookup(begin, len);
    //                 if (index != dictionary_type::invalid_index) {
    //                     write_index(index, out);
    //                     begin += len;
    //                     written += len;
    //                     break;
    //                 }
    //             }

    //             if (index == dictionary_type::invalid_index) {
    //                 begin += 1;
    //             }
    //         }

    //         return written;
    //     }

    //     static uint8_t const* decode(uint8_t const* in,
    //                                  uint32_t* out,
    //                                  uint32_t /*universe*/, size_t n,
    //                                  dictionary_type const* dict
    //                                  , dint_statistics& stats
    //                                  )
    //     {
    //         uint16_t const* ptr = reinterpret_cast<uint16_t const*>(in);
    //         for (size_t i = 0; i != n; ++ptr) {
    //             uint32_t index = *ptr;
    //             uint32_t decoded_ints = 1;
    //             if (DS2I_LIKELY(index > EXCEPTIONS - 1)) {
    //                 decoded_ints = dict->copy(index, out);
    //             } else {
    //                 if (index == 1) { // 4-byte exception
    //                     *out = *(reinterpret_cast<uint32_t const*>(++ptr));
    //                     ++ptr;

    //                     // stats.ints_distr[6] += 1;
    //                     // stats.codewords_distr[6] += 3;

    //                     // needed when b = 8
    //                     // ptr += 2;
    //                 } else { // 2-byte exception
    //                     // *out = *(reinterpret_cast<uint16_t const*>(++ptr)); // when b = 8
    //                     *out = *(++ptr);

    //                     // needed when b = 8
    //                     // ptr += 1;
    //                 }
    //             }
    //             out += decoded_ints;
    //             i += decoded_ints;
    //         }
    //         return reinterpret_cast<uint8_t const*>(ptr);
    //     }

    // private:
    //     static void write_index(uint32_t index, std::vector<uint8_t>& out) {
    //         auto ptr = reinterpret_cast<uint8_t const*>(&index);
    //         out.insert(out.end(), ptr, ptr + 2);
    //     }
    // };

    #define CODECS (interpolative)(optpfor)(varintg8iu)(qmx)(vbyte)(u32)(simple16)(streamvbyte)(maskedvbyte)(varintgb)(dint)
}
