// nei_plane.h
#pragma once

#include <glm/glm.hpp>
#include "nei_hittable.h"

struct Plane : public Hittable
{
    glm::vec3 m_Point;       // A point on the plane
    glm::vec3 m_Normal;      // The normal vector of the plane
    MaterialPtr m_Material;  // The material of the plane

    // Constructor initializing the plane with a point, normal, and material
    Plane(const glm::vec3& point, const glm::vec3& normal, const MaterialPtr& material)
        : m_Point(point), m_Normal(glm::normalize(normal)), m_Material(material) {}

    // Function to check if a ray hits the plane
    bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const override
    {
        float denom = glm::dot(m_Normal, ray.m_Direction);
        if (fabs(denom) > 1e-6f)
        {
            float t = glm::dot(m_Point - ray.m_Origin, m_Normal) / denom;
            if (t >= t_min && t <= t_max)
            {
                record.m_T = t;
                record.m_Point = ray.m_Origin + t * ray.m_Direction;
                record.m_Normal = denom < 0.0f ? m_Normal : -m_Normal;
                record.m_Material = m_Material;
                return true;
            }
        }
        return false;
    }

    // Function to get the bounding box of the plane
    bool BoundingBox(AABB& output_box) const override
    {
        // Planes are infinite, so no finite bounding box
        return false;
    }

    // Function to get the material of the plane
    MaterialPtr GetMaterial() const
    {
        return m_Material;
    }
};

// Type definition for a shared pointer to a Plane object
typedef std::shared_ptr<Plane> PlanePtr;
