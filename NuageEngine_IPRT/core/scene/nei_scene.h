// nei_scene.h
#pragma once

#include <vector>
#include <memory>
#include "nei_hittable.h"
#include "nei_bvh.h"

// Class representing the scene containing all objects
class Scene {
public:
    Scene() {}

    // Add an object to the scene
    void AddObject(const HittablePtr& object) {
        m_Objects.push_back(object);
    }

    // Build the BVH from the objects in the scene
    void BuildBVH() {
        if (m_Objects.empty()) {
            m_BVH = nullptr;
            return;
        }
        m_BVH = std::make_shared<BVHNode>(m_Objects, 0, m_Objects.size());
    }

    // Test the intersection of a ray with the scene
    bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const {
        if (!m_BVH)
            return false;
        return m_BVH->Hit(ray, t_min, t_max, record);
    }

    // Get the list of objects in the scene
    const std::vector<HittablePtr>& GetObjects() const {
        return m_Objects;
    }

private:
    std::vector<HittablePtr> m_Objects; // List of objects in the scene
    BVHNodePtr m_BVH;                    // Root of the BVH
};
