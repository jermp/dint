#pragma once

#include "FastPFor/headers/VarIntG8IU.h"
#include "FastPFor/headers/optpfor.h"
#include "FastPFor/headers/variablebyte.h"
#include "streamvbyte/include/streamvbyte.h"
#include "MaskedVByte/include/varintencode.h"
#include "MaskedVByte/include/varintdecode.h"
#include "varintgb.h"
#include "interpolative_coding.hpp"
#include "qmx_codec.hpp"
#include "succinct/util.hpp"
#include "util.hpp"
#include "dictionary.hpp"

#include "dint_configuration.hpp"

#include <limits>

namespace ds2i {

    // workaround: VariableByte::decodeArray needs the buffer size, while we
    // only know the number of values. It also pads to 32 bits. We need to
    // rewrite
    class TightVariableByte {
    public:
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

    struct interpolative_block {
      static const uint64_t block_size = constants::block_size;
      static const uint64_t overflow = 0;

      static void encode(uint32_t const *in, uint32_t sum_of_values, size_t n,
                         std::vector<uint8_t> &out) {
        assert(n <= block_size);
        thread_local std::vector<uint32_t> inbuf(block_size);
        thread_local std::vector<uint32_t> outbuf;
        inbuf[0] = *in;
        for (size_t i = 1; i < n; ++i) {
          inbuf[i] = inbuf[i - 1] + in[i];
        }

        if (sum_of_values == uint32_t(-1)) {
          sum_of_values = inbuf[n - 1];
          TightVariableByte::encode_single(sum_of_values, out);
        }

        bit_writer bw(outbuf);
        bw.write_interpolative(inbuf.data(), n - 1, 0, sum_of_values);
        uint8_t const *bufptr = (uint8_t const *)outbuf.data();
        out.insert(out.end(), bufptr,
                   bufptr + succinct::util::ceil_div(bw.size(), 8));
      }

      static uint8_t const *DS2I_NOINLINE decode(uint8_t const *in, uint32_t *out,
                                                 uint32_t sum_of_values, size_t n) {
        assert(n <= block_size);
        uint8_t const *inbuf = in;
        if (sum_of_values == uint32_t(-1)) {
          inbuf = TightVariableByte::decode(inbuf, &sum_of_values, 1);
        }

        out[n - 1] = sum_of_values;
        size_t read_interpolative = 0;
        if (n > 1) {
          bit_reader br((uint32_t const *)inbuf);
          br.read_interpolative(out, n - 1, 0, sum_of_values);
          for (size_t i = n - 1; i > 0; --i) {
            out[i] -= out[i - 1];
          }
          read_interpolative = succinct::util::ceil_div(br.position(), 8);
        }

        return inbuf + read_interpolative;
      }
    };

    struct optpfor_block {

      struct codec_type : FastPFor::OPTPFor<4, FastPFor::Simple16<false>> {

        uint8_t const *force_b;

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
      };

      static const uint64_t block_size = codec_type::BlockSize;
      static const uint64_t overflow = 0;

      static void encode(uint32_t const *in, uint32_t sum_of_values, size_t n,
                         std::vector<uint8_t> &out,
                         uint8_t const *b = nullptr) // if non-null forces b
      {
        thread_local codec_type optpfor_codec;
        thread_local std::vector<uint8_t> buf(2 * 4 * block_size);
        assert(n <= block_size);

        if (n < block_size) {
          interpolative_block::encode(in, sum_of_values, n, out);
          return;
        }

        size_t out_len = buf.size();

        optpfor_codec.force_b = b;
        optpfor_codec.encodeBlock(in, reinterpret_cast<uint32_t *>(buf.data()),
                                  out_len);
        out_len *= 4;
        out.insert(out.end(), buf.data(), buf.data() + out_len);
      }

      static uint8_t const *DS2I_NOINLINE decode(uint8_t const *in, uint32_t *out,
                                                 uint32_t sum_of_values, size_t n) {
        thread_local codec_type optpfor_codec; // pfor decoding is *not* thread-safe
        assert(n <= block_size);

        if (DS2I_UNLIKELY(n < block_size)) {
          return interpolative_block::decode(in, out, sum_of_values, n);
        }

        size_t out_len = block_size;
        uint8_t const *ret;

        ret = reinterpret_cast<uint8_t const *>(optpfor_codec.decodeBlock(
            reinterpret_cast<uint32_t const *>(in), out, out_len));
        assert(out_len == n);
        return ret;
      }
    };

    struct varint_G8IU_block {
      static const uint64_t block_size = constants::block_size;
      static const uint64_t overflow = 0;

      struct codec_type : FastPFor::VarIntG8IU {

        // rewritten version of decodeBlock optimized for when the output
        // size is known rather than the input
        // the buffers pointed by src and dst must be respectively at least
        // 9 and 8 elements large
        uint32_t decodeBlock(uint8_t const *&src, uint32_t *dst) const {
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
                data, vecmask[desc][1]); //__builtin_ia32_pshufb128(data, shf2);
            _mm_storeu_si128(
                reinterpret_cast<__m128i *>(dst + 4),
                result2); //__builtin_ia32_storedqu(dst + (16), result2);
          }

          return readSize;
        }
      };

      static void encode(uint32_t const *in, uint32_t sum_of_values, size_t n,
                         std::vector<uint8_t> &out) {
        thread_local codec_type varint_codec;
        thread_local std::vector<uint8_t> buf(2 * 4 * block_size);
        assert(n <= block_size);

        if (n < block_size) {
          interpolative_block::encode(in, sum_of_values, n, out);
          return;
        }

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
      static uint8_t const *decode(uint8_t const *in, uint32_t *out,
                                   uint32_t sum_of_values, size_t n) {
        static codec_type varint_codec; // decodeBlock is thread-safe
        assert(n <= block_size);

        if (DS2I_UNLIKELY(n < block_size)) {
          return interpolative_block::decode(in, out, sum_of_values, n);
        }

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

    struct qmx_block
    {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 512; // qmx can potentially overshoot...

        static void encode(uint32_t const *in, uint32_t sum_of_values, size_t n,
                           std::vector<uint8_t> &out) {
            assert(n <= block_size);
            if (n < block_size) {
                interpolative_block::encode(in, sum_of_values, n, out);
                return;
            }
            thread_local QMX::codec<block_size> qmx_codec;
            thread_local std::vector<uint8_t> buf( (overflow*4) + 2 * 4 * block_size);
            size_t out_len = buf.size();
            out_len = qmx_codec.encode(buf.data(),in);
            TightVariableByte::encode_single(out_len, out);
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        // we only allow varint to be inlined (others have DS2I_NOILINE)
        static uint8_t const *decode(uint8_t const *in, uint32_t *out,
                                   uint32_t sum_of_values, size_t n) {
            static QMX::codec<block_size> qmx_codec; // decodeBlock is thread-safe
            assert(n <= block_size);

            if (DS2I_UNLIKELY(n < block_size)) {
                return interpolative_block::decode(in, out, sum_of_values, n);
            }

            uint32_t enc_len = 0;
            in = TightVariableByte::decode(in, &enc_len, 1);
            qmx_codec.decode(out,in, enc_len);
            return in + enc_len;
        }
    };

    struct vbyte_block {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 0;

        static void encode(uint32_t const* in, uint32_t /* sum_of_values */,
            size_t n, std::vector<uint8_t>& out)
        {
            std::vector<uint8_t> buf(2 * 4 * block_size);
            size_t out_len = buf.size();
            TightVariableByte::encode(in, n, buf.data(), out_len);
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in, uint32_t* out,
            uint32_t /* sum_of_values */, size_t n)
        {
            return TightVariableByte::decode(in, out, n);
        }
    };

    struct u32_block {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 0;

        static void encode(uint32_t const* in, uint32_t /* sum_of_values */,
            size_t n, std::vector<uint8_t>& out)
        {
            assert(n <= block_size);
            size_t srclen = n * sizeof(uint32_t);
            const uint8_t* src = (const uint8_t*)in;
            out.insert(out.end(), src, src + srclen);
        }

        static uint8_t const* decode(uint8_t const* in, uint32_t* out,
            uint32_t /* sum_of_values */, size_t n)
        {
            assert(n <= block_size);
            const uint8_t* src = (const uint8_t*)in;
            uint8_t* dst = (uint8_t*)out;
            size_t n4 = n * 4;
            for (size_t i = 0; i < n4; i++) {
                *dst++ = *src++;
            }
            return src;
        }
    };

    struct simple16_block
    {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 512;
        using codec_type = FastPFor::Simple16<false>;

        static void encode(uint32_t const* in, uint32_t /* sum_of_values */,
            size_t n, std::vector<uint8_t>& out)
        {
            thread_local codec_type simple16_codec;
            assert(n <= block_size);
            // XXX this could be threadlocal static
            std::vector<uint8_t> buf(2 * 8 * block_size);
            size_t out_len = buf.size();
            simple16_codec.encodeArray(in, n, reinterpret_cast<uint32_t*>(buf.data()), out_len);
            out_len *= 4;
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in, uint32_t* out,
            uint32_t /* sum_of_values */, size_t n)
        {
            assert(n <= block_size);
            uint8_t const* ret;
            thread_local codec_type simple16_codec;
            ret = reinterpret_cast<uint8_t const*>(simple16_codec.decodeArray(reinterpret_cast<uint32_t const*>(in), 1,
                out, n));
            return ret;
        }
    };

    struct streamvbyte_block
    {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 512;

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out)
        {
            uint32_t *src = const_cast<uint32_t*>(in);
            std::vector<uint8_t> buf(streamvbyte_max_compressedbytes(n));
            size_t out_len = streamvbyte_encode(src, n, buf.data());
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n)
        {
            auto read = streamvbyte_decode(in, out, n);
            return in + read;
        }
    };

    struct maskedvbyte_block
    {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 512;

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out)
        {
            uint32_t* src = const_cast<uint32_t*>(in);
            std::vector<uint8_t> buf(2 * n * sizeof(uint32_t));
            size_t out_len = vbyte_encode(src, n, buf.data());
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n)
        {
            auto read = masked_vbyte_decode(in, out, n);
            return in + read;
        }
    };

    struct varintgb_block
    {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 512;

        static void encode(uint32_t const* in,
                           uint32_t /*universe*/, uint32_t n,
                           std::vector<uint8_t>& out)
        {
            thread_local VarIntGB<false> varintgb_codec;
            thread_local std::vector<uint8_t> buf(2 * n * sizeof(uint32_t));
            size_t out_len = varintgb_codec.encodeArray(in, n, buf.data());
            out.insert(out.end(), buf.data(), buf.data() + out_len);
        }

        static uint8_t const* decode(uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t /*universe*/, size_t n)
        {
            thread_local VarIntGB<false> varintgb_codec;
            auto read = varintgb_codec.decodeArray(in, n, out);
            return read + in;
        }
    };

    struct node
    {
        node()
        {}

        node(uint32_t p, uint32_t w, uint32_t c)
            : parent(p), codeword(w), cost(c)
        {}

        uint32_t parent;
        uint32_t codeword;
        uint32_t cost;
    };

    struct dint_block
    {
        static const uint64_t block_size = constants::block_size;
        static const uint64_t overflow = 512;

        // NOTE: greedy parsing
        // template<typename Builder>
        // static void encode(Builder& builder, // TODO: restore constness
        //                    uint32_t const* in,
        //                    uint32_t sum_of_values,
        //                    uint32_t n,
        //                    std::vector<uint8_t>& out)
        // {
        //     if (n < block_size) {
        //         interpolative_block::encode(in, sum_of_values, n, out);
        //         return;
        //     }

        //     uint32_t const* begin = in;
        //     uint32_t const* end = begin + n;

        //     uint32_t i = 0;
        //     uint32_t pos = 0;
        //     uint32_t cost = 0;

        //     while (begin < end)
        //     {
        //         uint32_t longest_run_size = 0;
        //         uint32_t run_size = std::min<uint64_t>(256, end - begin);
        //         uint32_t index = EXCEPTIONS;

        //         for (uint32_t const* ptr  = begin;
        //                              ptr != begin + run_size;
        //                            ++ptr)
        //         {
        //             if (*ptr == 0) {
        //                 ++longest_run_size;
        //             } else {
        //                 break;
        //             }
        //         }

        //         if (longest_run_size >= 16) {
        //             uint32_t k = 256;
        //             while (longest_run_size < k and k > 16) {
        //                 ++index;
        //                 k /= 2;
        //             }

        //             cost += 1;
        //             std::cout << "cost " << cost << "; pos " << pos << "; len " << k << "; codeword: " << index << std::endl;
        //             write_index(index, out);
        //             ++builder.codewords;
        //             begin += k;
        //             pos += k;

        //         } else {
        //             for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
        //                 uint32_t sub_block_size = constants::target_sizes[s];
        //                 uint32_t len = std::min<uint32_t>(sub_block_size, end - begin);
        //                 index = builder.lookup(begin, len);
        //                 if (index != Builder::invalid_index) {

        //                     cost += 1;
        //                     std::cout << "cost " << cost << "; pos " << pos << "; len " << len << "; codeword: " << index << std::endl;
        //                     write_index(index, out);
        //                     ++builder.codewords;
        //                     begin += len;
        //                     pos += len;

        //                     break;
        //                 }
        //             }

        //             if (index == Builder::invalid_index)
        //             {
        //                 uint32_t exception = *begin;
        //                 auto ptr = reinterpret_cast<uint8_t const*>(&exception);

        //                 // USED WITH EXCEPTIONS = 1
        //                 // out.insert(out.end(), 0);
        //                 // out.insert(out.end(), 0); // comment if b = 8
        //                 // out.insert(out.end(), ptr, ptr + 4);

        //                 if (exception < 65536) {

        //                     ++builder.small_exceptions;
        //                     out.insert(out.end(), 0);
        //                     out.insert(out.end(), 0); // comment if b = 8
        //                     out.insert(out.end(), ptr, ptr + 2);
        //                     cost += 2;

        //                     std::cout << "cost " << cost << "; pos " << pos << "; len 1; codeword: 0" << std::endl;

        //                 } else {

        //                     ++builder.large_exceptions;
        //                     out.insert(out.end(), 1);
        //                     out.insert(out.end(), 0); // comment if b = 8
        //                     out.insert(out.end(), ptr, ptr + 4);
        //                     cost += 3;

        //                     std::cout << "cost " << cost << "; pos " << pos << "; len 1; codeword: 1" << std::endl;
        //                 }

        //                 pos += 1;
        //                 begin += 1;
        //             }
        //         }

        //         ++i;
        //     }

        //     std::cout << "cost = " << cost << std::endl;
        // }

        // NOTE: optimal parsing
        template<typename Builder>
        static void encode(Builder& builder, // TODO: restore constness
                           uint32_t const* in,
                           uint32_t sum_of_values,
                           uint32_t n,
                           std::vector<uint8_t>& out)
        {
            if (n < block_size) {
                interpolative_block::encode(in, sum_of_values, n, out);
                return;
            }

            // NOTE: everything at the beginning is a large exception
            // costs are in shorts! (1 short = 16 bits)
            std::vector<node> path(n + 2);
            path[0] = {0, 1, 0}; // dummy node
            for (uint32_t i = 1; i < n + 1; ++i) {
                path[i] = {i - 1, 1, 3 * i};
            }

            // logger() << "finding shortest path" << std::endl;

            for (uint32_t i = 0; i < n; ++i)
            {
                // std::cout << "current node " << i << ":\n";
                // std::cout << "parent " << path[i].parent << "\n";
                // std::cout << "codeword " << path[i].codeword << "\n";
                // std::cout << "cost " << path[i].cost << "\n\n";

                uint32_t longest_run_size = 0;
                uint32_t run_size = std::min<uint64_t>(256, n - i);
                uint32_t index = EXCEPTIONS;

                for (uint32_t j = i; j != i + run_size; ++j) {
                    if (in[j] == 0) {
                        ++longest_run_size;
                    } else {
                        break;
                    }
                }

                // std::cout << "longest_run_size " << longest_run_size << std::endl;

                if (longest_run_size >= 16) {
                    uint32_t k = 256;
                    while (longest_run_size < k and k > 16) {
                        k /= 2;
                        ++index;
                    }
                    while (k >= 16) {
                        // std::cout << "k " << k << "; index " << index << std::endl;

                        uint32_t c = path[i].cost + 1;
                        if (path[i + k].cost > c) {
                            path[i + k] = {i, index, c};
                        }

                        k /= 2;
                        ++index;
                    }

                    // std::cout << std::endl;
                }

                for (uint32_t s = 0; s < constants::num_target_sizes; ++s) {
                    uint32_t sub_block_size = constants::target_sizes[s];
                    uint32_t len = std::min<uint32_t>(sub_block_size, n - i);
                    index = builder.lookup(in + i, len);
                    if (index != Builder::invalid_index) {
                        uint32_t c = path[i].cost + 1;
                        if (path[i + len].cost > c) {
                            path[i + len] = {i, index, c};
                        }
                    } else {
                        if (sub_block_size == 1) { // exceptions
                            uint32_t exception = in[i];
                            uint32_t c = path[i].cost + 2; // small exception cost
                            index = 0;

                            if (exception > 65536) {
                                c += 1; // large exception cost
                                index = 1;
                            }

                            if (path[i + 1].cost > c) {
                                path[i + 1] = {i, index, c};
                            }
                        }
                    }
                }
            }

            // std::cout << "min_cost " << path.back().cost << std::endl;

            std::vector<node> encoding;
            uint32_t i = n;
            while (i != 0) {
                uint32_t parent = path[i].parent;
                encoding.push_back(path[i]);
                i = parent;
            }

            // std::cout << std::endl;
            std::reverse(encoding.begin(), encoding.end());
            encoding.emplace_back(n, 1, inf); // final dummy node

            // logger() << "encoding" << std::endl;
            uint32_t pos = 0;
            uint32_t cost = 0;
            for (uint32_t i = 0; i < encoding.size() - 1; ++i) {
                uint32_t index = encoding[i].codeword;
                uint32_t len = encoding[i + 1].parent - encoding[i].parent;

                assert(len == builder.size(index));

                cost += 1;

                if (index > 1) {
                    ++builder.codewords;
                    write_index(index, out);
                } else {

                    assert(len == 1);
                    uint32_t exception = in[pos];
                    auto ptr = reinterpret_cast<uint8_t const*>(&exception);
                    cost += 1;

                    if (index == 0) {
                        ++builder.small_exceptions;
                        out.insert(out.end(), 0);
                        out.insert(out.end(), 0); // comment if b = 8
                        out.insert(out.end(), ptr, ptr + 2);
                    } else {
                        cost += 1;
                        ++builder.large_exceptions;
                        out.insert(out.end(), 1);
                        out.insert(out.end(), 0); // comment if b = 8
                        out.insert(out.end(), ptr, ptr + 4);
                    }
                }

                // std::cout << "cost " << cost << "; pos " << pos << "; len " << len << "; codeword: " << index << std::endl;

                pos += len;
            }

            assert(pos == n);

            // std::cout << "pos " << pos << "/" << block_size << std::endl;
            // std::cout << "cost = " << cost << std::endl;
        }

        template<typename Dictionary>
        static uint8_t const* decode(Dictionary const& dict,
                                     uint8_t const* in,
                                     uint32_t* out,
                                     uint32_t sum_of_values,
                                     size_t n)
        {
            if (DS2I_UNLIKELY(n < block_size)) {
                return interpolative_block::decode(in, out, sum_of_values, n);
            }

            uint16_t const* ptr = reinterpret_cast<uint16_t const*>(in); // b = 16
            // uint8_t const* ptr = in; // b = 8
            for (size_t i = 0; i != n; ++ptr)
            {
                uint32_t index = *ptr;
                uint32_t decoded_ints = 1;
                if (DS2I_LIKELY(index > Dictionary::reserved - 1)) {
                    decoded_ints = dict.copy(index, out);
                } else {
                    static const uint32_t run_lengths[] = {0, 1, // exceptions
                                                           256, 128, 64, 32, 16};
                    decoded_ints = run_lengths[index];

                    if (DS2I_UNLIKELY(decoded_ints == 1)) { // 4-byte exception
                        *out = *(reinterpret_cast<uint32_t const*>(++ptr));
                        ++ptr;

                        // needed when b = 8
                        // ptr += 2;
                    }

                    if (DS2I_UNLIKELY(decoded_ints == 0)) { // 2-byte exception
                        // *out = *(reinterpret_cast<uint16_t const*>(++ptr)); // when b = 8
                        *out = *(++ptr);
                        decoded_ints = 1;
                        // needed when b = 8
                        // ptr += 1;
                    }
                }

                out += decoded_ints;
                i += decoded_ints;
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
}
