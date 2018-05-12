#pragma once

#include <vector>
#include <algorithm>

namespace ds2i {

    template<typename T, typename C>
    struct heap {
        heap()
        {}

        heap(uint64_t size, C const& comparator)
            : heap()
        {
            init(size, comparator);
        }

        void init(uint64_t size, C const& comparator) {
            m_q.reserve(size);
            m_comparator = comparator;
        }

        void push(T const& x) {
            m_q.push_back(x);
            std::push_heap(m_q.begin(), m_q.end(), m_comparator);
        }

        T const& top() {
            return m_q.front();
        }

        void pop() {
            std::pop_heap(m_q.begin(), m_q.end(), m_comparator);
            m_q.pop_back();
        }

        void heapify() {
            sink(0);
        }

        void clear() {
            m_q.clear();
        }

        bool empty() const {
            return m_q.empty();
        }

        inline uint64_t size() const {
            return m_q.size();
        }

    private:
        uint64_t m_k;
        std::vector<T> m_q;
        C m_comparator;

        void sink(uint64_t pos) {
            assert(pos <= size());
            while (2 * pos + 1 < size()) {
                uint64_t i = 2 * pos + 1;
                if (i + 1 < size() and m_comparator(m_q[i], m_q[i + 1])) {
                    ++i;
                }
                if (not m_comparator(m_q[pos], m_q[i])) break;
                std::swap(m_q[pos], m_q[i]);
                pos = i;
            }
        }
    };

}
