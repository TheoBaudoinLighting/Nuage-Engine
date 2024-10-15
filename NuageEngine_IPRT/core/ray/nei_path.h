// nei_path.h
#pragma once

#include <glm/glm.hpp>
#include <memory>
#include "nei_material.h"

// Structure representing a vertex in a path
struct PathVertex {
    glm::vec3 position;    // Position of the vertex
    glm::vec3 normal;      // Normal at the vertex
    MaterialPtr material;  // Material at the vertex
    glm::vec3 emission;    // Emission at the vertex
    glm::vec3 throughput;  // Throughput at the vertex
};

// Typedef for a shared pointer to a PathVertex
typedef std::shared_ptr<PathVertex> PathVertexPtr;
