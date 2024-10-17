// nei_transform.h
#pragma once

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include "nei_hittable.h"
#include "nei_aabb.h"

class Transform : public Hittable {
public:
    Transform(const HittablePtr& object, const glm::mat4& transform)
        : m_Object(object), m_Transform(transform)
    {
        m_InverseTransform = glm::inverse(m_Transform);
        m_InverseTranspose = glm::transpose(m_InverseTransform);
    }

    virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const override 
    {
        glm::vec4 origin4 = glm::vec4(ray.m_Origin, 1.0f);
        glm::vec4 direction4 = glm::vec4(ray.m_Direction, 0.0f);

        glm::vec3 transformed_origin = glm::vec3(m_InverseTransform * origin4);
        glm::vec3 transformed_direction = glm::normalize(glm::vec3(m_InverseTransform * direction4));

        Ray transformed_ray(transformed_origin, transformed_direction);

        if (m_Object->Hit(transformed_ray, t_min, t_max, record))
        {
            // Transform the intersection point and normal back into world space
            glm::vec4 world_point = m_Transform * glm::vec4(record.m_Point, 1.0f);
            record.m_Point = glm::vec3(world_point) / world_point.w;

            glm::vec4 world_normal = m_InverseTranspose * glm::vec4(record.m_Normal, 0.0f);
            record.m_Normal = glm::normalize(glm::vec3(world_normal));

            return true;
        }

        return false;
    }

    virtual bool BoundingBox(AABB& output_box) const override 
    {
        AABB child_box;
        if (!m_Object->BoundingBox(child_box))
            return false;

        glm::vec3 min(FLT_MAX);
        glm::vec3 max(-FLT_MAX);

        for (int i = 0; i < 8; i++) {
            glm::vec3 corner(
                (i & 1) ? child_box.Max().x : child_box.Min().x,
                (i & 2) ? child_box.Max().y : child_box.Min().y,
                (i & 4) ? child_box.Max().z : child_box.Min().z
            );

            glm::vec3 transformed_corner = glm::vec3(m_Transform * glm::vec4(corner, 1.0f));

            min = glm::min(min, transformed_corner);
            max = glm::max(max, transformed_corner);
        }

        output_box = AABB(min, max);
        return true;
    }

    HittablePtr GetInnerObject() const { return m_Object; }

    glm::mat4 GetTransform() const { return m_Transform; }

    glm::mat4 GetInverseTransform() const { return m_InverseTransform; }

    glm::mat4 GetInverseTranspose() const { return m_InverseTranspose; }

    MaterialPtr GetMaterial() const override { return m_Object->GetMaterial(); }

private:
    HittablePtr m_Object;           // Inner object
    glm::mat4 m_Transform;          // Transformation matrix
    glm::mat4 m_InverseTransform;   // Inverse transformation matrix
    glm::mat4 m_InverseTranspose;   // Inverse transpose matrix
};

typedef std::shared_ptr<Transform> TransformPtr;
typedef std::shared_ptr<Transform> TransformPtr;
