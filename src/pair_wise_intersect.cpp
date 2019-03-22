#include <iostream>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <succinct/mapper.hpp>

#include "index_types.hpp"
#include "util.hpp"

typedef uint32_t term_id_type;
typedef std::vector<term_id_type> term_id_vec;

bool read_query(term_id_vec& ret, std::istream& is = std::cin) {
    ret.clear();
    std::string line;
    if (!std::getline(is, line))
        return false;
    std::istringstream iline(line);
    term_id_type term_id;
    while (iline >> term_id) {
        ret.push_back(term_id);
    }

    return true;
}

template <typename Enum>
static uint64_t intersect(uint64_t num_docs, std::vector<Enum>& enums,
                          std::vector<uint32_t>& out) {
    // increasing frequency
    if (enums[0].size() > enums[1].size()) {
        std::swap(enums[0], enums[1]);
    }

    uint64_t results = 0;
    uint64_t candidate = enums[0].docid();
    size_t i = 1;
    while (candidate < num_docs) {
        for (; i < 2; ++i) {
            enums[i].next_geq(candidate);
            if (enums[i].docid() != candidate) {
                candidate = enums[i].docid();
                i = 0;
                break;
            }
        }

        if (i == 2) {
            out[results] = candidate;
            ++results;
            enums[0].next();
            candidate = enums[0].docid();
            i = 1;
        }
    }

    return results;
}

template <typename Index>
void perftest(const char* index_filename) {
    using namespace ds2i;

    Index index;
    logger() << "Loading index from " << index_filename << std::endl;
    boost::iostreams::mapped_file_source m(index_filename);
    succinct::mapper::map(index, m);

    std::vector<term_id_vec> queries;
    term_id_vec q;
    while (read_query(q)) {
        assert(q.size() == 2);
        queries.push_back(q);
    }

    uint32_t num_queries = queries.size();
    logger() << "Executing " << num_queries << " pair-wise intersections..."
             << std::endl;

    uint64_t num_docs = index.num_docs();
    std::vector<uint32_t> out(num_docs);

    double total_usecs = 0.0;
    // first run if for warming up
    static const int runs = 10 + 1;
    size_t total = 0;

    typedef typename Index::document_enumerator enum_type;
    std::vector<enum_type> qq;
    qq.reserve(2);
    for (int run = 0; run != runs; ++run) {
        double start = get_time_usecs();
        for (uint32_t i = 0; i != num_queries; ++i) {
            qq.clear();
            for (auto term : queries[i]) {
                qq.push_back(index[term]);
            }
            uint64_t size = intersect(num_docs, qq, out);
            total += size;
        }
        double end = get_time_usecs();
        double elapsed = end - start;
        if (run) {
            total_usecs += elapsed;
        }
    }

    // for debug
    std::cout << total << std::endl;

    printf(
        "\t %d intersections took %lf [musecs] (avg. among %d "
        "runs)\n",
        num_queries, total_usecs / (runs - 1), runs - 1);
    printf(
        "\t %lf [musecs] per intersection (avg. among %d "
        "queries)\n",
        total_usecs / (runs - 1) / num_queries, num_queries);
}

int main(int argc, const char** argv) {
    using namespace ds2i;

    int mandatory = 3;
    if (argc < mandatory) {
        std::cerr << argv[0] << " <index_type> <index_filename> < query_log"
                  << std::endl;
        return 1;
    }

    std::string index_type = argv[1];
    const char* index_filename = argv[2];

    if (false) {
#define LOOP_BODY(R, DATA, T)                              \
    }                                                      \
    else if (index_type == BOOST_PP_STRINGIZE(T)) {        \
        perftest<BOOST_PP_CAT(T, _index)>(index_filename); \
        /**/

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, DS2I_INDEX_TYPES);
#undef LOOP_BODY
    } else {
        logger() << "ERROR: Unknown index type " << index_type << std::endl;
    }

    return 0;
}
