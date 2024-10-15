// nei_sphere.h
#pragma once

#include <glm/glm.hpp>
#include <memory>
#include "nei_ray.h"
#include "nei_hittable.h"

struct Sphere : public Hittable {
    glm::vec3 m_Center;       // Center of the sphere
    float m_Radius;           // Radius of the sphere
    MaterialPtr m_Material;   // Material of the sphere

    Sphere(const glm::vec3& center, float radius, const MaterialPtr& material)
        : m_Center(center), m_Radius(radius), m_Material(material) {}

    virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const override {
        glm::vec3 oc = ray.m_Origin - m_Center;
        float a = glm::dot(ray.m_Direction, ray.m_Direction);
        float b = glm::dot(oc, ray.m_Direction);
        float c = glm::dot(oc, oc) - m_Radius * m_Radius;
        float discriminant = b * b - a * c;

        if (discriminant > 0.0f) {
            float sqrt_disc = std::sqrt(discriminant);
            float temp = (-b - sqrt_disc) / a;
            if (temp < t_max && temp > t_min) {
                record.m_T = temp;
                record.m_Point = ray.m_Origin + record.m_T * ray.m_Direction;
                record.m_Normal = glm::normalize(record.m_Point - m_Center);
                record.m_Material = m_Material;
                return true;
            }
            temp = (-b + sqrt_disc) / a;
            if (temp < t_max && temp > t_min) {
                record.m_T = temp;
                record.m_Point = ray.m_Origin + record.m_T * ray.m_Direction;
                record.m_Normal = glm::normalize(record.m_Point - m_Center);
                record.m_Material = m_Material;
                return true;
            }
        }
        return false;
    }

    virtual bool BoundingBox(AABB& output_box) const override {
        output_box = AABB(
            m_Center - glm::vec3(m_Radius + 0.0001f),
            m_Center + glm::vec3(m_Radius + 0.0001f)
        );
        return true;
    }
};

typedef std::shared_ptr<Sphere> SpherePtr;
