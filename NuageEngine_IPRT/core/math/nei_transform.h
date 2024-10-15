// nei_transform.h
#pragma once

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include "nei_hittable.h"
#include "nei_aabb.h"

class Transform : public Hittable
{
public:
    Transform(const HittablePtr& object, const glm::mat4& transform)
        : m_Object(object), m_Transform(transform)
    {
        m_InverseTransform = glm::inverse(m_Transform);
        m_InverseTranspose = glm::transpose(m_InverseTransform);
    }

    // Method to access the inner object
    HittablePtr GetInnerObject() const { return m_Object; }

    virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const override
    {
        // Transform the ray into the local space of the object
        glm::vec4 origin4 = glm::vec4(ray.m_Origin, 1.0f);
        glm::vec4 direction4 = glm::vec4(ray.m_Direction, 0.0f);

        glm::vec3 transformed_origin = glm::vec3(m_InverseTransform * origin4);
        glm::vec3 transformed_direction = glm::normalize(glm::vec3(m_InverseTransform * direction4));

        Ray transformed_ray(transformed_origin, transformed_direction);

        if (m_Object->Hit(transformed_ray, t_min, t_max, record))
        {
            // Transform the hit point and normal back to world space
            glm::vec4 world_point = m_Transform * glm::vec4(record.m_Point, 1.0f);
            record.m_Point = glm::vec3(world_point) / world_point.w;

            glm::vec4 world_normal = m_InverseTranspose * glm::vec4(record.m_Normal, 0.0f);
            record.m_Normal = glm::normalize(glm::vec3(world_normal));

            return true;
        }

        return false;
    }

    bool BoundingBox(AABB& output_box) const override
    {
        AABB local_box;
        if (!m_Object->BoundingBox(local_box))
            return false;

        // Transform the corners of the bounding box
        glm::vec3 corners[8];
        corners[0] = glm::vec3(local_box.m_Min.x, local_box.m_Min.y, local_box.m_Min.z);
        corners[1] = glm::vec3(local_box.m_Max.x, local_box.m_Min.y, local_box.m_Min.z);
        corners[2] = glm::vec3(local_box.m_Min.x, local_box.m_Max.y, local_box.m_Min.z);
        corners[3] = glm::vec3(local_box.m_Max.x, local_box.m_Max.y, local_box.m_Min.z);
        corners[4] = glm::vec3(local_box.m_Min.x, local_box.m_Min.y, local_box.m_Max.z);
        corners[5] = glm::vec3(local_box.m_Max.x, local_box.m_Min.y, local_box.m_Max.z);
        corners[6] = glm::vec3(local_box.m_Min.x, local_box.m_Max.y, local_box.m_Max.z);
        corners[7] = glm::vec3(local_box.m_Max.x, local_box.m_Max.y, local_box.m_Max.z);

        glm::vec3 transformed_corners[8];
        for (int i = 0; i < 8; ++i)
        {
            glm::vec4 transformed = m_Transform * glm::vec4(corners[i], 1.0f);
            transformed_corners[i] = glm::vec3(transformed) / transformed.w;
        }

        glm::vec3 new_min = transformed_corners[0];
        glm::vec3 new_max = transformed_corners[0];
        for (int i = 1; i < 8; ++i)
        {
            new_min = glm::min(new_min, transformed_corners[i]);
            new_max = glm::max(new_max, transformed_corners[i]);
        }

        output_box = AABB(new_min, new_max);
        return true;
    }

    glm::mat4 GetTransform() const { return m_Transform; }
    glm::mat4 GetInverseTransform() const { return m_InverseTransform; }
    glm::mat4 GetInverseTranspose() const { return m_InverseTranspose; }

private:
    HittablePtr m_Object;
    glm::mat4 m_Transform;
    glm::mat4 m_InverseTransform;
    glm::mat4 m_InverseTranspose;
};

typedef std::shared_ptr<Transform> TransformPtr;
