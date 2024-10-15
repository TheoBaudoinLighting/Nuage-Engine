// nei_bvh.h
#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include "nei_hittable.h"
#include "nei_aabb.h"

class BVHNode : public Hittable {
public:
    BVHNode();

    BVHNode(const std::vector<HittablePtr>& objects, size_t start, size_t end);

    virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const override;

    virtual bool BoundingBox(AABB& output_box) const override;

private:
    HittablePtr m_Left;
    HittablePtr m_Right;
    AABB m_Box;
};

typedef std::shared_ptr<BVHNode> BVHNodePtr;
