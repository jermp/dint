#pragma once

#include <succinct/mapper.hpp>
#include "util.hpp"

using ds2i::logger;

template <typename InputCollection, typename Collection>
void verify_collection(InputCollection const& input, const char* filename)
{
    Collection coll;
    boost::iostreams::mapped_file_source m(filename);
    succinct::mapper::map(coll, m);

    DS2I_LOG << "Checking the written data, just to be extra safe...";
    size_t s = 0;
    for (auto seq: input) {
        auto e = coll[s];
        if (e.size() != seq.docs.size()) {
            DS2I_LOG << "sequence " << s
                     << " has wrong length! ("
                     << e.size() << " != " << seq.docs.size() << ")"
                    ;
            exit(1);
        }

        for (size_t i = 0; i < e.size(); ++i, e.next()) {
            uint64_t docid = *(seq.docs.begin() + i);
            uint64_t freq = *(seq.freqs.begin() + i);

            if (docid != e.docid()) {
                DS2I_LOG << "docid in sequence " << s
                         << " differs at position " << i << "!";
                DS2I_LOG << e.docid() << " != " << docid;
                DS2I_LOG << "sequence length: " << seq.docs.size();

                exit(1);
            }

            if (freq != e.freq()) {
                DS2I_LOG << "freq in sequence " << s
                         << " differs at position " << i << "!";
                DS2I_LOG << e.freq() << " != " << freq;
                DS2I_LOG << "sequence length: " << seq.docs.size();

                exit(1);
            }
        }

        s += 1;
    }
    DS2I_LOG << "Everything is OK!";
}

