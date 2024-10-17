// nei_utility.h
#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <thread>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

const float PI = glm::pi<float>();
const float EPSILON = 1e-4f;
const int MAX_DEPTH = 6;

// XorShift Random Number Generator
struct XorShift {
    uint32_t state;

    XorShift(uint32_t seed = 0) : state(seed) {
        if (state == 0) {
            state = 1; 
        }
    }

    inline uint32_t next() {
        uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
        return x;
    }

    inline float next_float() {
        return (next() & 0xFFFFFF) / 16777216.0f; 
    }
};

// Sobol sequence generator class
class Sobol {
public:
    Sobol(int dimensions = 2) : m_Dimension(dimensions), m_Index(0) {
        if (m_Dimension > max_dimensions) {
            std::cerr << "Sobol sequence supports up to " << max_dimensions << " dimensions." << std::endl;
            m_Dimension = max_dimensions;
        }

        initialize_direction_numbers();

        for (int dim = 0; dim < m_Dimension; ++dim) {
            X[dim] = 0;
        }
    }

    glm::dvec2 next() {
        glm::dvec2 sample(0.0);

        for (int dim = 0; dim < m_Dimension; ++dim) {
            int c = count_trailing_zeros(m_Index + 1);
            if (c >= direction_numbers[dim].size()) {
                c = direction_numbers[dim].size() - 1;
            }
            X[dim] ^= direction_numbers[dim][c];
            sample[dim] = static_cast<double>(X[dim]) / static_cast<double>(1ULL << max_bits);
        }

        ++m_Index;
        return sample;
    }

private:
    int m_Dimension;
    uint32_t m_Index;
    static const int max_dimensions = 2; 
    static const int max_bits = 32; 

    std::vector<std::vector<uint64_t>> direction_numbers = std::vector<std::vector<uint64_t>>(max_dimensions, std::vector<uint64_t>());

    uint64_t X[max_dimensions];

    int count_trailing_zeros(uint32_t n) const {
        if (n == 0) return max_bits;
        int count = 0;
        while ((n & 1) == 0) {
            n >>= 1;
            ++count;
        }
        return count;
    }

    void initialize_direction_numbers() 
    {
        for (int i = 1; i <= max_bits; ++i) {
            direction_numbers[0].push_back(1ULL << (max_bits - i));
        }

        direction_numbers[1].resize(max_bits);
        uint64_t v = 1ULL << (max_bits - 1);
        for (int i = 0; i < max_bits; ++i) {
            direction_numbers[1][i] = v;
            v ^= v >> 1;
        }
    }
};

inline float random_float() {
    thread_local static XorShift gen_local(std::random_device{}());
    return gen_local.next_float();
}

inline glm::vec3 random_in_unit_sphere() {
    while (true) {
        glm::vec3 p(random_float() * 2.0f - 1.0f,
                    random_float() * 2.0f - 1.0f,
                    random_float() * 2.0f - 1.0f);
        if (glm::dot(p, p) < 1.0f)
            return p;
    }
}

// Utility transformations
inline glm::mat4 translation(const glm::vec3& offset) {
    return glm::translate(glm::mat4(1.0f), offset);
}

inline glm::mat4 rotation(float angle, const glm::vec3& axis) {
    return glm::rotate(glm::mat4(1.0f), angle, axis);
}

inline glm::mat4 scale(const glm::vec3& factors) {
    return glm::scale(glm::mat4(1.0f), factors);
}
