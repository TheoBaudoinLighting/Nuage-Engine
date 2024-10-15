#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <string>

#include "nei_path.h"
#include "nei_scene.h"
#include "nei_camera.h"
#include "nei_renderer_utility.h"
#include "nei_utility.h"

class Renderer
{
public:
    Renderer(const RendererParameters& params);

    void Render(const Scene& scene, const Camera& camera);

private:
    const RendererParameters& m_Params;

    // Private rendering functions
    glm::vec3 double_trace(const Ray& ray, const Scene& scene, XorShift& gen_local);
    std::vector<PathVertex> trace_eye_path(const Ray& ray, const Scene& scene, XorShift& gen_local);
    std::vector<PathVertex> trace_light_path(const Scene& scene, XorShift& gen_local);
    glm::vec3 connect_paths(const PathVertex& eye_vertex, const PathVertex& light_vertex, const Scene& scene);

    // Utility functions for rendering
    void update_statistics(glm::vec3& mean, glm::vec3& M2, int n, const glm::vec3& new_sample);
    std::string generate_unique_filename(const std::string& base, const std::string& ext);

    // Storage for rendered pixels
    std::vector<glm::vec3> m_Pixels;
};
