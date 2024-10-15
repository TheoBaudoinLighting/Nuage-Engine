// nei_hittable.h
#pragma once

#include <glm/glm.hpp>
#include <memory>

#include "nei_aabb.h"
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
    // Pure virtual function to check if a ray hits the object
    virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const = 0;
    // Pure virtual function to get the bounding box of the object
    virtual bool BoundingBox(AABB& box) const = 0;
};

// Type definition for a shared pointer to a Hittable object
typedef std::shared_ptr<Hittable> HittablePtr;
