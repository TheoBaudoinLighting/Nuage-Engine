// nei_material.h
#pragma once

#include <glm/glm.hpp>
#include <memory>

class Material {
public:
    glm::vec3 m_Albedo;     // Reflection color
    glm::vec3 m_Emission;   // Emission color
    bool m_IsEmissive;      // Indicates if the material is emissive

    // Constructor
    Material(const glm::vec3& albedo = glm::vec3(0.0f),
             const glm::vec3& emission = glm::vec3(0.0f))
        : m_Albedo(albedo), m_Emission(emission), m_IsEmissive(emission != glm::vec3(0.0f)) 
        {
            m_IsEmissive = glm::length(m_Emission) > 1e-6f;
        }
};

typedef std::shared_ptr<Material> MaterialPtr;
