// nei_camera.h
#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "nei_ray.h"

class Camera
{
public:
    Camera(
        glm::vec3 position = glm::vec3(0.0f, 0.0f, 1.5f),
        glm::vec3 look_at = glm::vec3(0.0f, 0.0f, -1.0f),
        glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f),
        float fov = 60.0f,
        float aspect_ratio = 1.0f)
    {
        m_Position = position;
        m_Fov = fov;
        m_AspectRatio = aspect_ratio;

        // Calculate the viewport dimensions
        float theta = glm::radians(m_Fov);
        float h = tan(theta / 2.0f);
        float viewport_height = 2.0f * h;
        float viewport_width = m_AspectRatio * viewport_height;

        // Calculate camera basis vectors
        m_W = glm::normalize(m_Position - look_at);
        m_U = glm::normalize(glm::cross(up, m_W));
        m_V = glm::cross(m_W, m_U);

        // Calculate the horizontal and vertical vectors of the viewport
        m_Horizontal = m_U * viewport_width;
        m_Vertical = m_V * viewport_height;

        // Calculate the lower left corner of the viewport
        m_LowerLeftCorner = m_Position - m_Horizontal / 2.0f - m_Vertical / 2.0f - m_W;
    }

    // Generate a ray from the camera through the viewport
    Ray GetRay(float u, float v) const
    {
        glm::vec3 direction = m_LowerLeftCorner + u * m_Horizontal + v * m_Vertical - m_Position;
        return Ray(m_Position, glm::normalize(direction));
    }

private:
    glm::vec3 m_Position;         // Camera position
    glm::vec3 m_LowerLeftCorner;  // Lower left corner of the viewport
    glm::vec3 m_Horizontal;       // Horizontal vector of the viewport
    glm::vec3 m_Vertical;         // Vertical vector of the viewport
    float m_Fov;                  // Field of view
    float m_AspectRatio;          // Aspect ratio
    glm::vec3 m_W;                // Camera basis vector W
    glm::vec3 m_U;                // Camera basis vector U
    glm::vec3 m_V;                // Camera basis vector V
};