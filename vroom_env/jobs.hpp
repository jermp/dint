#pragma once

#include "semiasync_queue.hpp"

namespace ds2i {

    // for generic block-codec
    template<typename Iterator, typename Encoder>
    struct enc_sequence_adder : semiasync_queue::job {
        enc_sequence_adder(
            Iterator begin,
            uint64_t n,
            boost::progress_display& progress,
            std::vector<uint8_t>& output,
            bool docs,
            uint64_t& num_processed_lists,
            uint64_t& num_total_ints
        )
            : begin(begin)
            , n(n)
            , universe(0)
            , progress(progress)
            , output(output)
            , docs(docs)
            , num_processed_lists(num_processed_lists)
            , num_total_ints(num_total_ints)
        {}

        virtual void prepare()
        {
            std::vector<uint32_t> buf;
            buf.reserve(n);
            uint32_t prev = docs ? -1 : 0;
            for (uint64_t i = 0; i != n; ++i, ++begin) {
                buf.push_back(*begin - prev - 1);
                if (docs) prev = *begin;
                universe += buf.back();
            }
            assert(buf.size() == n);

            Encoder::encode(
                buf.data(), universe, buf.size(), tmp
            );
        }

        virtual void commit() {
            header::write(n, universe, output);
            output.insert(output.end(), tmp.begin(), tmp.end());
            progress += n + 1;
            ++num_processed_lists;
            num_total_ints += n;
        }

        Iterator begin;
        uint64_t n;
        uint32_t universe;
        boost::progress_display& progress;
        std::vector<uint8_t> tmp;
        std::vector<uint8_t>& output;
        bool docs;
        uint64_t& num_processed_lists;
        uint64_t& num_total_ints;
    };

    // for DINT
    template<typename Iterator, typename Encoder, typename Builder>
    struct dint_sequence_adder : semiasync_queue::job {
        dint_sequence_adder(
            Iterator begin,
            uint64_t n,
            Builder& builder,
            boost::progress_display& progress,
            std::vector<uint8_t>& output,
            bool docs,
            uint64_t& num_processed_lists,
            uint64_t& num_total_ints
        )
            : begin(begin)
            , n(n)
            , universe(0)
            , builder(builder)
            , progress(progress)
            , output(output)
            , docs(docs)
            , num_processed_lists(num_processed_lists)
            , num_total_ints(num_total_ints)
        {}

        virtual void prepare()
        {
            std::vector<uint32_t> buf;
            buf.reserve(n);
            uint32_t prev = docs ? -1 : 0;
            for (uint64_t i = 0; i != n; ++i, ++begin) {
                buf.push_back(*begin - prev - 1);
                if (docs) prev = *begin;
                universe += buf.back();
            }
            assert(buf.size() == n);

            Encoder::encode(
                builder, buf.data(), universe, buf.size(), tmp
            );
        }

        virtual void commit() {
            header::write(n, universe, output);
            output.insert(output.end(), tmp.begin(), tmp.end());
            progress += n + 1;
            ++num_processed_lists;
            num_total_ints += n;
        }

        Iterator begin;
        uint64_t n;
        uint32_t universe;
        Builder& builder;
        boost::progress_display& progress;
        std::vector<uint8_t> tmp;
        std::vector<uint8_t>& output;
        bool docs;
        uint64_t& num_processed_lists;
        uint64_t& num_total_ints;
    };

    // for PEF
    template<typename Iterator>
    struct pef_sequence_adder : semiasync_queue::job {
        pef_sequence_adder(
            Iterator begin,
            uint64_t n, uint64_t universe,
            succinct::bit_vector_builder& bvb,
            boost::progress_display& progress,
            bool docs,
            uint64_t& num_processed_lists,
            uint64_t& num_total_ints
        )
            : begin(begin)
            , n(n)
            , universe(universe)
            , bvb(bvb)
            , progress(progress)
            , docs(docs)
            , num_processed_lists(num_processed_lists)
            , num_total_ints(num_total_ints)
        {}

        virtual void prepare()
        {
            if (not docs) { // on freqs, PEF needs the perfix sums
                universe = 0;
                auto in = begin;
                for (uint64_t i = 0; i != n; ++i, ++in) {
                    universe += *in;
                }
                universe += 1;
            }
            pef::encode(begin, universe, n, tmp, not docs);
        }

        virtual void commit()
        {
            progress += n + 1;
            ++num_processed_lists;
            num_total_ints += n;
            uint64_t offset = bvb.size() + tmp.size();
            bvb.append_bits(offset + 64 + 32 + 32, 64);
            bvb.append_bits(universe, 32);
            bvb.append_bits(n, 32);
            bvb.append(tmp);
        }

        Iterator begin;
        uint64_t n;
        uint64_t universe;
        succinct::bit_vector_builder& bvb;
        succinct::bit_vector_builder tmp;
        boost::progress_display& progress;
        bool docs;
        uint64_t& num_processed_lists;
        uint64_t& num_total_ints;
    };

}
