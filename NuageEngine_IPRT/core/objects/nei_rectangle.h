// nei_rectangle.h
#pragma once

#include <glm/glm.hpp>
#include "nei_hittable.h"
#include "nei_material.h"
#include "nei_aabb.h"
#include "nei_utility.h" // Include this for XorShift

class Rectangle : public Hittable
{
public:
    Rectangle(const glm::vec3& corner, const glm::vec3& edge1, const glm::vec3& edge2, const MaterialPtr& material)
        : m_Corner(corner), m_Edge1(edge1), m_Edge2(edge2), m_Material(material)
    {
       m_Normal = -glm::normalize(glm::cross(m_Edge1, m_Edge2));
    }

    bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const override
    {
        float denom = glm::dot(m_Normal, ray.m_Direction);
        if (fabs(denom) > 1e-6f)
        {
            float t = glm::dot(m_Corner - ray.m_Origin, m_Normal) / denom;
            if (t >= t_min && t <= t_max)
            {
                glm::vec3 p = ray.m_Origin + t * ray.m_Direction;
                glm::vec3 d = p - m_Corner;

                float dot11 = glm::dot(m_Edge1, m_Edge1);
                float dot12 = glm::dot(m_Edge1, m_Edge2);
                float dot22 = glm::dot(m_Edge2, m_Edge2);
                float dot1p = glm::dot(m_Edge1, d);
                float dot2p = glm::dot(m_Edge2, d);

                float denom_bary = dot11 * dot22 - dot12 * dot12;
                float u = (dot22 * dot1p - dot12 * dot2p) / denom_bary;
                float v = (dot11 * dot2p - dot12 * dot1p) / denom_bary;

                if (u >= 0.0f && u <= 1.0f && v >= 0.0f && v <= 1.0f)
                {
                    record.m_T = t;
                    record.m_Point = p;
                    record.m_Normal = denom < 0.0f ? m_Normal : -m_Normal;
                    record.m_Material = m_Material;

                    return true;
                }
            }
        }
        return false;
    }

    bool BoundingBox(AABB& output_box) const override
    {
        glm::vec3 corner0 = m_Corner;
        glm::vec3 corner1 = m_Corner + m_Edge1;
        glm::vec3 corner2 = m_Corner + m_Edge2;
        glm::vec3 corner3 = m_Corner + m_Edge1 + m_Edge2;

        glm::vec3 min_corner = glm::min(glm::min(corner0, corner1), glm::min(corner2, corner3));
        glm::vec3 max_corner = glm::max(glm::max(corner0, corner1), glm::max(corner2, corner3));

        output_box = AABB(min_corner, max_corner);
        return true;
    }

    // Method to sample a point on the rectangle
    glm::vec3 SamplePoint(XorShift& gen_local) const {
        float u = gen_local.next_float();
        float v = gen_local.next_float();
        return m_Corner + u * m_Edge1 + v * m_Edge2;
    }

    // Area of the rectangle
    float Area() const {
        return glm::length(glm::cross(m_Edge1, m_Edge2));
    }

    // Method to get the normal of the rectangle
    glm::vec3 GetNormal() const {
        return m_Normal;
    }

    // Method to get the material of the rectangle
    MaterialPtr GetMaterial() const {
        return m_Material;
    }

private:
    glm::vec3 m_Corner;      // One corner of the rectangle
    glm::vec3 m_Edge1;       // One edge vector
    glm::vec3 m_Edge2;       // Another edge vector
    glm::vec3 m_Normal;      // Normal of the rectangle
    MaterialPtr m_Material;  // Material of the rectangle
};

typedef std::shared_ptr<Rectangle> RectanglePtr;
