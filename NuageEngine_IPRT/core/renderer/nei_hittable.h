// nei_hittable.h
#pragma once

#include <glm/glm.hpp>

#include "nei_aabb.h"
#include <memory>
#include "nei_material.h"

// Structure to store hit record details
struct HitRecord {
    float m_T; // Time of hit
    glm::vec3 m_Point; // Point of intersection
    glm::vec3 m_Normal; // Normal at the intersection
    MaterialPtr m_Material; // Material of the intersected object
};

struct Ray;

// Abstract class for hittable objects
struct Hittable {
    virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const = 0;
    virtual bool BoundingBox(AABB& output_box) const = 0;
};

// Type definition for a shared pointer to a Hittable object
typedef std::shared_ptr<Hittable> HittablePtr;
