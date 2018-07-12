#include <iostream>
#include <vector>

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
              int d, int n, int k)
{
    if (d == k + 1) return;

    for (int i = 1; i <= n; ++i)
    {
        prefix[d - 1] = targets[i - 1];
        if (d == k)
        {
            for (int p = 0; p < k; ++p) {
                std::cout << prefix[p] << " ";
            }
            std::cout << "\n";

            std::cout << std::endl;

            int l = 0;
            int sum = 0;
            for (int p = 0; p < n; ++p) {
                if (targets[p] > l) {
                    l = targets[p];
                }
            }
            for (int p = 0; p < k; ++p) {
                sum += prefix[p];
            }
            std::vector<std::vector<int>> matrix(l, std::vector<int>(k, sum));
            l = 0;
            for (int r = 0; r < k; ++r) {
                int codeword_length = prefix[r];
                for (int p = 0; p < codeword_length; ++p) {
                    auto& v = matrix[p];
                    v[r] = l++;
                }
            }

            for (auto const& v: matrix) {
                for (auto x: v) {
                    std::cout << x << " ";
                }
                std::cout << std::endl;
            }

            std::cout << std::endl;
            std::cout << std::endl;
        }

        generate(prefix, targets, d + 1, n, k);
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
    for (int i = 0; i < n; ++i) {
        int x;
        std::cin >> x;
        targets[i] = x;
    }

    // for (auto x: targets) {
    //     std::cout << x << " ";
    // }
    // std::cout << std::endl;

    // generate(1, n, k);
    std::vector<int> prefix(k);
    std::cout << std::endl;
    generate(prefix, targets, 1, n, k);

    return 0;
}

