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

extern std::mt19937 gen;
extern std::uniform_real_distribution<float> dis;

// Define a more optimized XorShift class
struct XorShift {
    uint32_t state;

    // Constructor with an initial seed
    XorShift(uint32_t seed = 0) : state(seed) {
        if (state == 0) {
            state = 1; // Avoid a zero initial state
        }
    }

    // Generate the next random number
    inline uint32_t next() {
        uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
        return x;
    }

    // Generate a floating-point number in [0, 1)
    inline float next_float() {
        return (next() & 0xFFFFFF) / 16777216.0f; // 2^24 = 16777216
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
            return p; // Remove glm::normalize(p)
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
