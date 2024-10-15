// nei_scene.h
#pragma once

#include <vector>
#include <memory>
#include "nei_hittable.h"
#include "nei_bvh.h"

class Scene
{
public:
    Scene() {}

    void AddObject(const HittablePtr& object)
    {
        AABB box;
        if (object->BoundingBox(box)) {
            m_FiniteObjects.push_back(object);
        } else {
            m_InfiniteObjects.push_back(object);
        }
    }

    const std::vector<HittablePtr>& GetFiniteObjects() const {
        return m_FiniteObjects;
    }

    const std::vector<HittablePtr>& GetInfiniteObjects() const {
        return m_InfiniteObjects;
    }

    void BuildBVH()
    {
        if (!m_FiniteObjects.empty())
            m_BVH = std::make_shared<BVHNode>(m_FiniteObjects, 0, m_FiniteObjects.size());
    }

    bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const
    {
        HitRecord temp_rec;
        bool hit_anything = false;
        float closest_so_far = t_max;

        // Test finite objects directly
        for (const auto& object : m_FiniteObjects) {
            if (object->Hit(ray, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.m_T;
                record = temp_rec;
            }
        }

        // Test infinite objects
        for (const auto& object : m_InfiniteObjects) {
            if (object->Hit(ray, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.m_T;
                record = temp_rec;
            }
        }

        return hit_anything;
    }

private:
    std::vector<HittablePtr> m_FiniteObjects;
    std::vector<HittablePtr> m_InfiniteObjects;
    std::shared_ptr<BVHNode> m_BVH;
};
