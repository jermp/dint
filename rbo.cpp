#include <iostream>
#include <vector>
#include <fstream>

#include <boost/format.hpp>

#include <unordered_set>

#include "util.hpp"

template<typename RandomAccessIterator>
void rbo(RandomAccessIterator A, RandomAccessIterator B,
         uint64_t n, double p)
{
    double weight = 1.0 - p;
    double rbo_min = 0.0;
    double contrib = 0.0;
    uint64_t overlap = 0;

    std::unordered_set<uint32_t> seen_A;
    std::unordered_set<uint32_t> seen_B;
    size_t i = 1;

    for (; i <= n; ++i) {
        if (seen_A.count(A[i - 1]) == 1 or seen_B.count(B[i - 1]) == 1) {
            ds2i::logger() << "ERROR: Duplicate value at line " << i << std::endl;
        }

        if (A[i - 1] == B[i - 1]) {
            overlap++;
        } else {
            overlap += seen_B.count(A[i - 1]);
            overlap += seen_A.count(B[i - 1]);
        }
        seen_A.insert(A[i - 1]);
        seen_B.insert(B[i - 1]);
        contrib = weight * overlap / i;
        rbo_min += contrib;
        weight *= p;
    }

    auto max_overlap = overlap;
    auto rbo_max = rbo_min;
    const double EPSILON = 1e-15;
    while (weight > EPSILON) {
        i++;
        contrib = weight * overlap / i;
        rbo_min += contrib;
        if (max_overlap == i - 1) {
            // both new elements must be the same novel one
            max_overlap += 1;
        } else {
            // two new elements can be assumed, both ones that
            // appeared already in the other list
            max_overlap += 2;
        }
        rbo_max += weight * max_overlap / i;
        // prepare for the next pair of imaginary values
        weight *= p;
    }
    boost::format rbofmt("\t rbo(p = %1$.5f) = %2$.6f + %3$.6f (n=%4$7d, d=%5$7d)");
    ds2i::logger() << rbofmt % p % rbo_min % (rbo_max - rbo_min) % n % (i - 1) << std::endl;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t<filename1> <filename2>" << std::endl;
        std::cerr << "The files should contain one integer per line." << std::endl;
        return 1;
    }

    std::string filename_A = argv[1];
    std::string filename_B = argv[2];
    std::ifstream in_A(filename_A.c_str());
    std::ifstream in_B(filename_B.c_str());
    std::vector<uint32_t> A;
    std::vector<uint32_t> B;

    uint32_t x;
    while (in_A >> x) A.push_back(x);
    while (in_B >> x) B.push_back(x);

    if (A.size() != B.size()) {
        throw std::runtime_error("input files should contain the same number of lines");
    }

    auto P = {0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 0.999, 0.9999, 0.99999};
    for (auto p: P) {
        rbo(A.begin(), B.begin(), A.size(), p);
    }

    return 0;
}
