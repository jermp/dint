#pragma once

#include "util.hpp"

namespace ds2i {

void print_statistics(std::string type, char const* encoded_data_filename,
                      std::vector<double> const& timings,
                      uint64_t num_decoded_ints, uint64_t num_decoded_lists
                      // , dint_statistics const& stats
) {
    static const uint64_t billion = 1000000000;
    double tot_elapsed =
        std::accumulate(timings.begin(), timings.end(), double(0.0));
    double ns_x_int = tot_elapsed * billion / num_decoded_ints;
    uint64_t ints_x_sec = uint64_t(1 / ns_x_int * billion);

    logger() << "elapsed time " << tot_elapsed << " [sec]" << std::endl;
    logger() << ns_x_int << " [ns] x int" << std::endl;
    logger() << ints_x_sec << " ints x [sec]" << std::endl;

    // logger() << "avg. # of decoded integers x codeword: " <<
    // double(stats.decoded_ints_from_dict) / stats.dict_codewords << std::endl;

    // stats to std output
    std::cout << "{";
    std::cout << "\"filename\": \"" << encoded_data_filename << "\", ";
    std::cout << "\"num_sequences\": \"" << num_decoded_lists << "\", ";
    std::cout << "\"num_integers\": \"" << num_decoded_ints << "\", ";
    std::cout << "\"type\": \"" << type << "\", ";
    std::cout << "\"tot_elapsed_time\": \"" << tot_elapsed << "\", ";
    std::cout << "\"ns_x_int\": \"" << ns_x_int << "\", ";
    std::cout << "\"ints_x_sec\": \"" << ints_x_sec << "\"";
    std::cout << "}" << std::endl;

    // // binary entropy
    // {
    //     uint64_t total_exceptions = 0;

    //     uint64_t total_codewords = 0;
    //     for (uint64_t i = 0; i != constants::num_entries; ++i) {
    //         total_codewords += stats.occs[i];
    //     }

    //     for (auto const& pair: stats.exceptions) {
    //         total_exceptions += pair.second;
    //     }

    //     logger() << "total_exceptions " << total_exceptions << std::endl;

    //     double exceptions_entropy = 0.0;
    //     for (auto const& pair: stats.exceptions) {
    //         double p = double(pair.second) / total_exceptions;
    //         exceptions_entropy += p * std::log2(1.0 / p);
    //     }

    //     double entropy = 0.0;
    //     for (uint64_t i = 0; i != constants::num_entries; ++i) {
    //         double x = stats.occs[i];
    //         if (x != 0) {
    //             double p = x / total_codewords;
    //             entropy += p * std::log2(1.0 / p);
    //         }
    //     }

    //     logger() << "total_ints " << stats.total_ints << std::endl;
    //     logger() << "binary entropy of the source: " << (entropy *
    //     total_codewords + exceptions_entropy * total_exceptions) / 8.0 /
    //     constants::GiB << " [GiB]" << std::endl; logger() << "binary entropy
    //     x codeword: " << entropy << " bits" << std::endl; logger() << "binary
    //     entropy x exception: " << exceptions_entropy << " bits" << std::endl;
    //     logger() << "binary entropy x integer: " << (entropy *
    //     total_codewords + exceptions_entropy * total_exceptions) /
    //     stats.total_ints << " bits" << std::endl;
    // }

    // uint64_t total_codewords = 0;
    // uint64_t total_decoded_ints = 0;
    // for (uint64_t i = 0; i < stats.codewords_distr.size(); ++i) {
    //     total_codewords += stats.codewords_distr[i];
    //     total_decoded_ints += stats.ints_distr[i];
    // }

    // std::cout << "total_codewords " << total_codewords << std::endl;
    // std::cout << "total_decoded_ints " << total_decoded_ints << std::endl;

    // for (uint64_t i = 0; i < stats.codewords_distr.size(); ++i) {
    //     if (i == 0) std::cout << "freq:\n";
    //     else if (i == stats.codewords_distr.size() - 1) std::cout <<
    //     "rare:\n"; else std::cout << (uint32_t(1) << (i - 1)) << ":\n";
    //     std::cout << "\t codewords: " << stats.codewords_distr[i] * 100.0 /
    //     total_codewords << "%" << std::endl; std::cout << "\t integers: " <<
    //     stats.ints_distr[i] * 100.0 / total_decoded_ints << "%" << std::endl;
    // }
}

}  // namespace ds2i