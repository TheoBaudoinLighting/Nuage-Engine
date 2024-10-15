// nei_aabb.h
#pragma once

#include <glm/glm.hpp>
#include "nei_ray.h"

class AABB {
public:
    AABB() {}
    AABB(const glm::vec3& min, const glm::vec3& max) : m_Min(min), m_Max(max) {}

    const glm::vec3& Min() const { return m_Min; }
    const glm::vec3& Max() const { return m_Max; }

    bool Hit(const Ray& ray, float t_min, float t_max) const {
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

private:
    glm::vec3 m_Min;
    glm::vec3 m_Max;
};

inline AABB SurroundingBox(const AABB& box0, const AABB& box1) {
    glm::vec3 small(fmin(box0.Min().x, box1.Min().x),
                    fmin(box0.Min().y, box1.Min().y),
                    fmin(box0.Min().z, box1.Min().z));

    glm::vec3 big(fmax(box0.Max().x, box1.Max().x),
                  fmax(box0.Max().y, box1.Max().y),
                  fmax(box0.Max().z, box1.Max().z));

    return AABB(small, big);
}
