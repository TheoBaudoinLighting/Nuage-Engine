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

    float Area() const {
        return 4 * glm::pi<float>() * m_Radius * m_Radius;
    }

    float Pdf(const glm::vec3& origin, const glm::vec3& direction) const {
        HitRecord rec;
        Ray ray(origin, direction);
        if (!this->Hit(ray, EPSILON, std::numeric_limits<float>::max(), rec)) {
            return 0.0f;
        }

        float area = this->Area();
        float distance_squared = glm::length(rec.m_Point - origin);
        float cosine = glm::dot(rec.m_Normal, -direction);

        return distance_squared / (cosine * area);
    }

    virtual bool BoundingBox(AABB& output_box) const override {
        output_box = AABB(
            m_Center - glm::vec3(m_Radius + 0.0001f),
            m_Center + glm::vec3(m_Radius + 0.0001f)
        );
        return true;
    }

    MaterialPtr GetMaterial() const override { return m_Material; }
};

typedef std::shared_ptr<Sphere> SpherePtr;
