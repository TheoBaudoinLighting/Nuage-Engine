// nei_aabb.h
#pragma once

#include <glm/glm.hpp>
#include "nei_ray.h"

class AABB
{
public:
    glm::vec3 m_Min;
    glm::vec3 m_Max;

    // Default constructor initializing min and max to zero vectors
    AABB() : m_Min(glm::vec3(0.0f)), m_Max(glm::vec3(0.0f)) {}

    // Constructor initializing min and max with given vectors
    AABB(const glm::vec3& min, const glm::vec3& max) : m_Min(min), m_Max(max) {}

    // Static method to create a surrounding box from two AABBs
    static AABB surrounding_box(const AABB& box0, const AABB& box1)
    {
        glm::vec3 small = glm::min(box0.m_Min, box1.m_Min);
        glm::vec3 big = glm::max(box0.m_Max, box1.m_Max);
        return AABB(small, big);
    }

    // Method to check if a ray hits the AABB
    bool Hit(const Ray& ray, float t_min, float t_max) const
    {
        for (int a = 0; a < 3; a++) {
            float invD = 1.0f / ray.m_Direction[a];
            float t0 = (m_Min[a] - ray.m_Origin[a]) * invD;
            float t1 = (m_Max[a] - ray.m_Origin[a]) * invD;
            if (invD < 0.0f)
                std::swap(t0, t1);
            t_min = t0 > t_min ? t0 : t_min;
            t_max = t1 < t_max ? t1 : t_max;
            if (t_max <= t_min)
                return false;
        }
        return true;
    }

    // New method to calculate the centroid
    glm::vec3 Centroid() const {
        return 0.5f * (m_Min + m_Max);
    }
};
