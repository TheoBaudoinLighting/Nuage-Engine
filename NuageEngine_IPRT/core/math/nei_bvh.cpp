// nei_bvh.cpp
#include "nei_bvh.h"
#include <iostream>
#include <algorithm>

bool box_compare(const HittablePtr& a, const HittablePtr& b, int axis) {
    AABB box_a;
    AABB box_b;

    if (!a->BoundingBox(box_a) || !b->BoundingBox(box_b)) {
        std::cerr << "No bounding box in bvh_node constructor.\n";
        return false;
    }

    float a_centroid = box_a.Centroid()[axis];
    float b_centroid = box_b.Centroid()[axis];

    return a_centroid < b_centroid;
}

bool box_x_compare(const HittablePtr& a, const HittablePtr& b) {
    return box_compare(a, b, 0);
}

bool box_y_compare(const HittablePtr& a, const HittablePtr& b) {
    return box_compare(a, b, 1);
}

bool box_z_compare(const HittablePtr& a, const HittablePtr& b) {
    return box_compare(a, b, 2);
}

BVHNode::BVHNode(std::vector<HittablePtr>& objects, size_t start, size_t end) {
    auto& objs = objects;

    // Initialize centroid_bounds with the centroid of the first object
    AABB temp_box;
    if (!objs[start]->BoundingBox(temp_box)) {
        std::cerr << "No bounding box in bvh_node constructor.\n";
    }
    glm::vec3 centroid = temp_box.Centroid();
    AABB centroid_bounds(centroid, centroid);

    // Expand centroid_bounds to include the centroids of all objects
    for (size_t i = start + 1; i < end; ++i) {
        if (!objs[i]->BoundingBox(temp_box)) {
            std::cerr << "No bounding box in bvh_node constructor.\n";
        }
        centroid = temp_box.Centroid();
        centroid_bounds = AABB::surrounding_box(centroid_bounds, AABB(centroid, centroid));
    }

    glm::vec3 extent = centroid_bounds.m_Max - centroid_bounds.m_Min;
    int axis = extent.x > extent.y ? (extent.x > extent.z ? 0 : 2) : (extent.y > extent.z ? 1 : 2);

    auto comparator = (axis == 0) ? box_x_compare
                    : (axis == 1) ? box_y_compare
                                  : box_z_compare;

    size_t object_span = end - start;

    if (object_span == 1) {
        left = right = objs[start];
    } else if (object_span == 2) {
        if (comparator(objs[start], objs[start + 1])) {
            left = objs[start];
            right = objs[start + 1];
        } else {
            left = objs[start + 1];
            right = objs[start];
        }
    } else {
        std::sort(objs.begin() + start, objs.begin() + end, comparator);

        size_t mid = start + object_span / 2;
        left = std::make_shared<BVHNode>(objs, start, mid);
        right = std::make_shared<BVHNode>(objs, mid, end);
    }

    AABB box_left, box_right;

    if (!left->BoundingBox(box_left) || !right->BoundingBox(box_right)) {
        std::cerr << "No bounding box in bvh_node constructor.\n";
    }

    box = AABB::surrounding_box(box_left, box_right);
}

bool BVHNode::Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const {
    if (!box.Hit(ray, t_min, t_max))
        return false;

    HitRecord left_rec, right_rec;
    bool hit_left = left->Hit(ray, t_min, t_max, left_rec);
    bool hit_right = right->Hit(ray, t_min, t_max, right_rec);

    if (hit_left && hit_right) {
        if (left_rec.m_T < right_rec.m_T) {
            record = left_rec;
        } else {
            record = right_rec;
        }
        return true;
    } else if (hit_left) {
        record = left_rec;
        return true;
    } else if (hit_right) {
        record = right_rec;
        return true;
    } else {
        return false;
    }
}

bool BVHNode::BoundingBox(AABB& output_box) const {
    output_box = box;
    return true;
}
