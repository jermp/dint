#include <iostream>
#include <vector>
#include <cmath>

// n = # of possible target sizes
// k = # of words read at a time
// generate a complete n-ary tree of height k

// void generate(int d, int n, int k)
// {
//     if (d == k + 1) return;

//     for (int i = 1; i <= n; ++i) {
//         for (int p = 1; p < d; ++p) {
//             std::cout << "\t";
//         }
//         std::cout << i << "\n";
//         generate(d + 1, n, k);
//     }
// }

void generate(std::vector<int>& prefix,
              std::vector<int> const& targets,
              int& i, int N, int l,
              int d, int n, int k)
{
    if (d == k + 1) return;

    for (int j = 1; j <= n; ++j)
    {
        prefix[d - 1] = targets[j - 1];
        if (d == k)
        {
            // for (int p = 0; p < k; ++p) {
            //     std::cout << prefix[p] << " ";
            // }
            // std::cout << "\n";

            // std::cout << "\n";

            int sum = 0;
            for (int p = 0; p < k; ++p) {
                sum += prefix[p];
            }
            std::vector<std::vector<int>> matrix(l, std::vector<int>(k, sum));

            int index = 0;
            for (int r = 0; r < k; ++r) {
                int codeword_length = prefix[r];
                for (int p = 0; p < codeword_length; ++p) {
                    auto& v = matrix[p];
                    v[r] = index++;
                }
            }

            for (int t = 0; t < l; ++t) {
                std::cout << "\t{";
                auto const& v = matrix[t];
                for (int p = 0; p < k; ++p) {
                    std::cout << v[p];
                    if (p == k - 1) {
                        std::cout << "}";
                        if (i != N - 1) {
                            std::cout << ",";
                        }
                    } else {
                        std::cout << ", ";
                    }
                    ++i;
                }
                std::cout << "\n";
            }

            std::cout << "\n";
        }

        generate(prefix, targets, i, N, l, d + 1, n, k);
    }
}

int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << ":\n"
                  << "\t<n> <k>"
                  << std::endl;
        return 1;
    }

    int n = std::atoi(argv[1]);
    int k = std::atoi(argv[2]);
    std::vector<int> targets(n);
    std::cout << "enter targets' sizes:\n";

    int l = 0; // max. target size
    for (int i = 0; i < n; ++i) {
        int x;
        std::cin >> x;
        targets[i] = x;
        if (x > l) l = x;
    }

    int N = std::pow(n, k) * l;
    std::cout << "static uint32_t indices[" << N << "][" << k << "] = {" << std::endl;

    // for (auto x: targets) {
    //     std::cout << x << " ";
    // }
    // std::cout << std::endl;

    // generate(1, n, k);
    std::vector<int> prefix(k);
    std::cout << std::endl;
    int i = 0;
    generate(prefix, targets, i, N * l, l, 1, n, k);

    std::cout << "};" << std::endl;

    return 0;
}

