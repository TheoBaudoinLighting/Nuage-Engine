#pragma once

#include "nei_hittable.h"
#include "nei_aabb.h"
#include <memory>
#include <vector>

// Class representing a node in the Bounding Volume Hierarchy (BVH)
class BVHNode : public Hittable
{
public:
    // Default constructor
    BVHNode() {}

    // Constructor to build a BVH node from a list of objects
    BVHNode(std::vector<HittablePtr>& objects, size_t start, size_t end);

    // Function to check if a ray hits the BVH node
    bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const override;

    // Function to get the bounding box of the BVH node
    bool BoundingBox(AABB& output_box) const override;

private:
    HittablePtr left;  // Left child node
    HittablePtr right; // Right child node
    AABB box;          // Bounding box of the node
};

// Type definition for a shared pointer to a BVHNode object
typedef std::shared_ptr<BVHNode> BVHNodePtr;
