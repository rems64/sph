#pragma once

#include <cstddef>
#include <vector>
template <typename T> class RingBuffer {
public:
    RingBuffer(size_t n) : m_count(n), m_values(n), m_last_index(n){};

    const size_t count() const { return m_count; }
    void push_back(T v) {
        m_last_index = (m_last_index + 1) % m_count;
        m_values[m_last_index] = v;
    }
    const T *ptr() const { return m_values.data(); }
    size_t offset() const { return m_last_index; }

private:
    std::vector<T> m_values;
    size_t m_count;
    size_t m_last_index;
};
