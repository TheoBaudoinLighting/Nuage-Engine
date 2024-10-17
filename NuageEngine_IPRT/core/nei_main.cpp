// MIT License
// 
// Copyright (c) 2024 Théo Baudoin-Malnoë
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// nei_main.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <atomic>
#include <mutex>
#include <omp.h>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
//#include <glm/gtx/string_cast.hpp>

#include "nei_renderer.h"
#include "nei_material.h"
#include "nei_rectangle.h"
#include "nei_sphere.h"
#include "nei_plane.h"
#include "nei_transform.h"
#include "nei_camera.h"
#include "nei_scene.h"

int main()
{
    RendererParameters params;
    params.image_width = 640;
    params.image_height = 360;
    params.samples_per_pixel = 128;  // Min = 4
    params.max_depth = 3;
    params.russian_roulette = 0.5f;

    // Create a scene
    Scene scene;

    // Create materials
    auto red = std::make_shared<Material>(glm::vec3(0.65f, 0.05f, 0.05f));
    auto green = std::make_shared<Material>(glm::vec3(0.12f, 0.45f, 0.15f));
    auto white = std::make_shared<Material>(glm::vec3(0.73f, 0.73f, 0.73f));
    auto light = std::make_shared<Material>(glm::vec3(1.0f), glm::vec3(10.0f, 10.0f, 10.0f));
    light->m_IsEmissive = true;

    // Create planes
    // Left wall (green)
    auto left_wall = std::make_shared<Rectangle>(
        glm::vec3(-1.0f, -1.0f, -1.0f),
        glm::vec3(0.0f, 2.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 2.0f),
        green
    );
    scene.AddObject(left_wall);

    // Right wall (red)
    auto right_wall = std::make_shared<Rectangle>(
        glm::vec3(1.0f, -1.0f, -1.0f),
        glm::vec3(0.0f, 2.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 2.0f),
        red
    );
    scene.AddObject(right_wall);

    // Floor
    auto floor = std::make_shared<Rectangle>(
        glm::vec3(-1.0f, -1.0f, -1.0f),
        glm::vec3(2.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 2.0f),
        white
    );
    scene.AddObject(floor);

    // Ceiling
    auto ceiling = std::make_shared<Rectangle>(
        glm::vec3(-1.0f, 1.0f, -1.0f),
        glm::vec3(2.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 2.0f),
        white
    );
    scene.AddObject(ceiling);

    // Back wall
    auto back_wall = std::make_shared<Rectangle>(
        glm::vec3(-1.0f, -1.0f, -1.0f),
        glm::vec3(2.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 2.0f, 0.0f),
        white
    );
    scene.AddObject(back_wall);

    // Add rectangular light
    glm::vec3 light_corner(-0.5f, 1.0f, -0.5f); 
    glm::vec3 light_edge1(1.0f, 0.0f, 0.0f);     
    glm::vec3 light_edge2(0.0f, 0.0f, 1.0f);     
    auto rectangle_light = std::make_shared<Rectangle>(light_corner, light_edge1, light_edge2, light);
    scene.AddObject(rectangle_light);

    // Add a non-emissive white sphere (optional)
    auto sphere = std::make_shared<Sphere>(
        glm::vec3(0.0f, 0.0f, -0.5f),
        0.3f,
        white
    );
    glm::mat4 transform = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f, 1.0f, 1.0f));
    auto transformed_sphere = std::make_shared<Transform>(sphere, transform);
    scene.AddObject(transformed_sphere);

    // Create camera
    Camera cam(
        glm::vec3(0.0f, 0.0f, 3.0f),
        glm::vec3(0.0f, 0.0f, -0.5f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        40.0f,
        16.0f / 9.0f
    );

    // End of scene setup

    // Build BVH
    scene.BuildBVH();

    // Create a Renderer 
    Renderer renderer(params);

    // Launch the renderer
    renderer.Render(scene, cam);

    return 0;
}
