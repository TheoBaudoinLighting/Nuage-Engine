// nei_bvh.cpp
#include "nei_bvh.h"
#include <iostream>
#include <random>

// Default constructor
BVHNode::BVHNode() {}

// Constructor that builds the BVH from a list of objects
BVHNode::BVHNode(const std::vector<HittablePtr>& objects, size_t start, size_t end) {
    auto objs = objects; // Copy for manipulation

    // Select a random axis (0 = x, 1 = y, 2 = z)
    int axis = std::random_device{}() % 3;

    // Comparator to sort objects along the selected axis
    auto comparator = [axis](const HittablePtr& a, const HittablePtr& b) -> bool {
        AABB box_a, box_b;
        if (!a->BoundingBox(box_a) || !b->BoundingBox(box_b)) {
            std::cerr << "Error: No bounding box in BVH construction.\n";
        }
        return box_a.Min()[axis] < box_b.Min()[axis];
    };

    size_t object_span = end - start;

    if (object_span == 1) {
        m_Left = m_Right = objs[start];
    }
    else if (object_span == 2) {
        if (comparator(objs[start], objs[start + 1])) {
            m_Left = objs[start];
            m_Right = objs[start + 1];
        }
        else {
            m_Left = objs[start + 1];
            m_Right = objs[start];
        }
    }
    else {
        std::sort(objs.begin() + start, objs.begin() + end, comparator);

        size_t mid = start + object_span / 2;
        m_Left = std::make_shared<BVHNode>(objs, start, mid);
        m_Right = std::make_shared<BVHNode>(objs, mid, end);
    }

    AABB box_left, box_right;

    if (!m_Left->BoundingBox(box_left) || !m_Right->BoundingBox(box_right)) {
        std::cerr << "Error: No bounding box in BVH construction.\n";
    }

    m_Box = SurroundingBox(box_left, box_right);
}

// Method to test intersection with a ray
bool BVHNode::Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const {
    if (!m_Box.Hit(ray, t_min, t_max))
        return false;

    bool hit_left = m_Left->Hit(ray, t_min, t_max, record);
    bool hit_right = m_Right->Hit(ray, t_min, hit_left ? record.m_T : t_max, record);

    return hit_left || hit_right;
}

// Method to get the bounding box of the node
bool BVHNode::BoundingBox(AABB& output_box) const {
    output_box = m_Box;
    return true;
}
