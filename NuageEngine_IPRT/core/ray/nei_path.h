// nei_path.h
#pragma once

#include <glm/glm.hpp>
#include <memory>

#include "nei_hittable.h"
#include "nei_material.h"

struct PathVertex {
    glm::vec3 position;       // Position of the vertex
    glm::vec3 normal;         // Normal at the vertex
    MaterialPtr material;     // Material of the vertex
    glm::vec3 emission;       // Emission at the vertex
    glm::vec3 throughput;     // Throughput at the vertex
    glm::vec3 direct_light;   // Direct light at the vertex
    bool is_light;            // Flag to check if the vertex is a light source
    HittablePtr geometry;     // Geometry of the vertex
};

typedef std::shared_ptr<PathVertex> PathVertexPtr;
