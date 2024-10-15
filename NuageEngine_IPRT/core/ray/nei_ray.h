// nei_ray.h
#pragma once

#include <glm/glm.hpp>

// Structure representing a ray with an origin and direction
struct Ray
{
    glm::vec3 m_Origin;      // Origin of the ray
    glm::vec3 m_Direction;   // Direction of the ray

    // Default constructor initializing the ray with default values
    Ray() : m_Origin(glm::vec3(0.0f)), m_Direction(glm::vec3(0.0f, 0.0f, -1.0f)) {}

    // Constructor initializing the ray with given origin and direction
    Ray(const glm::vec3& origin, const glm::vec3& direction) : m_Origin(origin), m_Direction(glm::normalize(direction)) {}
};
