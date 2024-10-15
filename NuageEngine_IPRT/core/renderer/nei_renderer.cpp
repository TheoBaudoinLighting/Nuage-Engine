// nei_renderer.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <omp.h>

#include "nei_renderer.h"
#include "nei_utility.h"
#include "nei_ray.h"
#include "nei_path.h"
#include "nei_scene.h"
#include "nei_camera.h"
#include "nei_hittable.h"
#include "nei_material.h"
#include "nei_plane.h"
#include "nei_rectangle.h"
#include "nei_sphere.h"
#include "nei_transform.h"

// Constructor for the Renderer class
Renderer::Renderer(const RendererParameters& params)
    : m_Params(params), m_Pixels(params.image_width * params.image_height, glm::vec3(0.0f))
{
}

void Renderer::Render(const Scene& scene, const Camera& camera)
{
    // Rendering preparation
    int image_width = m_Params.image_width;
    int image_height = m_Params.image_height;
    int samples_per_pixel = m_Params.samples_per_pixel;
    int max_depth = m_Params.max_depth;

    // Multithreading management with OpenMP
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    omp_set_nested(0);
    std::cerr << "Using " << num_threads << " threads with nested parallelism disabled.\n";

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cerr << "Starting parallel rendering...\n";

    // Parallel rendering loop
#pragma omp parallel for schedule(guided)
    for (int j = 0; j < image_height; ++j) {
        XorShift gen_local(std::random_device{}() + omp_get_thread_num());

        for (int i = 0; i < image_width; ++i) {
            glm::vec3 mean_color(0.0f);
            glm::vec3 M2(0.0f);
            int current_samples = 0;
            const int max_samples = samples_per_pixel;

            while (current_samples < max_samples) {
                float u = (i + gen_local.next_float()) / static_cast<float>(image_width - 1);
                float v = (j + gen_local.next_float()) / static_cast<float>(image_height - 1);

                Ray r = camera.GetRay(u, v);
                glm::vec3 sample_color = double_trace(r, scene, gen_local);

                current_samples++;
                update_statistics(mean_color, M2, current_samples, sample_color);
            }

            // Use mean_color as the final pixel color
            glm::vec3 pixel_color = mean_color;

            // Apply gamma correction (gamma 2.2)
            float r_val = std::sqrt(pixel_color.r);
            float g_val = std::sqrt(pixel_color.g);
            float b_val = std::sqrt(pixel_color.b);
            m_Pixels[j * image_width + i] = glm::vec3(r_val, g_val, b_val);
        }

        // Display progress
        int completed = j + 1;
        if (completed % (image_height / 100) == 0) {
#pragma omp critical
            {
                std::cerr << "Progress: " << (100.0 * completed / image_height) << "%\n" << std::flush;
            }
        }
    }

    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> render_duration = end_time - start_time;
    double render_seconds = render_duration.count();
    std::cerr << "Render time: " << render_seconds << " seconds\n";

    // Save the rendered image
    std::cerr << "Writing pixels to file...\n";
    std::string filename = generate_unique_filename("render", "ppm");
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "Error: Unable to open output file.\n";
        return;
    }
    ofs << "P6\n" << image_width << " " << image_height << "\n255\n";

    std::vector<unsigned char> binary_pixels;
    binary_pixels.reserve(image_width * image_height * 3);
    for (int j = image_height - 1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            glm::vec3 color = m_Pixels[j * image_width + i];
            unsigned char ir = static_cast<unsigned char>(std::clamp(color.r, 0.0f, 0.999f) * 255.0f);
            unsigned char ig = static_cast<unsigned char>(std::clamp(color.g, 0.0f, 0.999f) * 255.0f);
            unsigned char ib = static_cast<unsigned char>(std::clamp(color.b, 0.0f, 0.999f) * 255.0f);
            binary_pixels.push_back(ir);
            binary_pixels.push_back(ig);
            binary_pixels.push_back(ib);
        }
    }
    ofs.write(reinterpret_cast<char*>(binary_pixels.data()), binary_pixels.size());
    ofs.close();
    std::cerr << "Image saved as '" << filename << "'.\n";
}
// Function to generate a cosine-weighted sample direction
inline glm::vec3 cosine_weighted_sample(const glm::vec3& normal, const glm::vec3& u, const glm::vec3& v, XorShift& gen_local) {
    float r1 = gen_local.next_float();
    float r2 = gen_local.next_float();
    float phi = 2.0f * PI * r1;
    float sqrt_r2 = std::sqrt(r2);
    float x = std::cos(phi) * sqrt_r2;
    float y = std::sin(phi) * sqrt_r2;
    float z = std::sqrt(1.0f - r2);

    return u * x + v * y + normal * z;
}

// Function to multiply a vector by a scalar
inline glm::vec3 multiply(const glm::vec3& vec, int scalar) {
    return vec * static_cast<float>(scalar);
}

// Function to compute the variance of a set of samples
inline float compute_variance(const std::vector<glm::vec3>& samples, glm::vec3 mean) {
    float variance_r = 0.0f;
    float variance_g = 0.0f;
    float variance_b = 0.0f;
    for (const auto& sample : samples) {
        float diff_r = sample.r - mean.r;
        float diff_g = sample.g - mean.g;
        float diff_b = sample.b - mean.b;
        variance_r += diff_r * diff_r;
        variance_g += diff_g * diff_g;
        variance_b += diff_b * diff_b;
    }
    variance_r /= samples.size();
    variance_g /= samples.size();
    variance_b /= samples.size();
    return (variance_r + variance_g + variance_b) / 3.0f;
}

// Function to update mean and variance using Welford's algorithm
void Renderer::update_statistics(glm::vec3& mean, glm::vec3& M2, int n, const glm::vec3& new_sample) {
    glm::vec3 delta = new_sample - mean; // Calculate the difference between the new sample and the current mean
    mean += delta / static_cast<float>(n); // Update the mean
    glm::vec3 delta2 = new_sample - mean; // Calculate the difference between the new sample and the updated mean
    M2 += delta * delta2; // Update the M2 value
}

// Trace a path from the camera
std::vector<PathVertex> Renderer::trace_eye_path(const Ray& ray, const Scene& scene, XorShift& gen_local) {
    std::vector<PathVertex> path;
    Ray current_ray = ray;
    glm::vec3 throughput(1.0f);
    int depth = 0;
    int max_depth = m_Params.max_depth;

    while (depth < max_depth) {
        HitRecord rec;
        if (scene.Hit(current_ray, EPSILON, std::numeric_limits<float>::max(), rec)) 
        {
            PathVertex vertex;
            vertex.position = rec.m_Point;
            vertex.normal = rec.m_Normal;
            vertex.material = rec.m_Material;
            vertex.emission = rec.m_Material->m_Emission;
            vertex.throughput = throughput;
            path.push_back(vertex);

            if (rec.m_Material->m_IsEmissive) {
                break;
            }

            glm::vec3 w = rec.m_Normal;
            glm::vec3 a = (fabs(w.x) > 0.1f) ? glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
            glm::vec3 v = glm::normalize(glm::cross(a, w));
            glm::vec3 u = glm::cross(w, v);

            glm::vec3 new_direction = cosine_weighted_sample(w, u, v, gen_local);
            float cos_theta = glm::dot(rec.m_Normal, new_direction);
            glm::vec3 brdf = rec.m_Material->m_Albedo / glm::pi<float>();

            throughput *= brdf * cos_theta;

            current_ray = Ray(rec.m_Point + EPSILON * new_direction, new_direction);

            depth++;
        } else {
            break;
        }
    }

    return path;
}

// Trace a path from the light
std::vector<PathVertex> Renderer::trace_light_path(const Scene& scene, XorShift& gen_local)
{
    std::vector<std::shared_ptr<Hittable>> emissive_objects;

    // Retrieve emissive objects from finite objects
    for (const auto& object : scene.GetFiniteObjects()) {
        // Rectangle
        auto rect = std::dynamic_pointer_cast<Rectangle>(object);
        if (rect && rect->GetMaterial()->m_IsEmissive) {
            emissive_objects.push_back(rect);
            continue;
        }

        // Sphere
        auto sphere = std::dynamic_pointer_cast<Sphere>(object);
        if (sphere && sphere->m_Material->m_IsEmissive) {
            emissive_objects.push_back(sphere);
            continue;
        }

        // Transform
        auto transform = std::dynamic_pointer_cast<Transform>(object);
        if (transform) {
            auto child = transform->GetInnerObject();
            if (child) {
                auto child_rect = std::dynamic_pointer_cast<Rectangle>(child);
                if (child_rect && child_rect->GetMaterial()->m_IsEmissive) {
                    emissive_objects.push_back(child_rect);
                    continue;
                }

                auto child_sphere = std::dynamic_pointer_cast<Sphere>(child);
                if (child_sphere && child_sphere->m_Material->m_IsEmissive) {
                    emissive_objects.push_back(child_sphere);
                    continue;
                }

                // Add other types of emissive objects if necessary
            }
        }
    }

    // Retrieve emissive objects from infinite objects
    for (const auto& object : scene.GetInfiniteObjects()) {
        // Plane
        auto plane = std::dynamic_pointer_cast<Plane>(object);
        if (plane && plane->GetMaterial()->m_IsEmissive) {
            emissive_objects.push_back(plane);
            continue;
        }

        // Add other types of infinite emissive objects if necessary
    }

    if (emissive_objects.empty()) {
        return {};
    }

    // Select a random emissive object
    size_t light_index = gen_local.next() % emissive_objects.size();
    auto selected_light = emissive_objects[light_index];

    glm::vec3 sampled_point;
    glm::vec3 sampled_normal;
    glm::vec3 emission;

    if (auto rect = std::dynamic_pointer_cast<Rectangle>(selected_light)) {
        sampled_point = rect->SamplePoint(gen_local);
        sampled_normal = rect->GetNormal();
        emission = rect->GetMaterial()->m_Emission;
    }
    else if (auto sphere = std::dynamic_pointer_cast<Sphere>(selected_light)) {
        float u = gen_local.next_float();
        float v = gen_local.next_float();
        float theta = 2.0f * PI * u;
        float phi = std::acos(1.0f - 2.0f * v);
        float x = std::sin(phi) * std::cos(theta);
        float y = std::sin(phi) * std::sin(theta);
        float z = std::cos(phi);
        sampled_normal = glm::vec3(x, y, z);
        sampled_point = sphere->m_Center + sampled_normal * sphere->m_Radius;
        emission = sphere->m_Material->m_Emission;
    }
    else if (auto plane = std::dynamic_pointer_cast<Plane>(selected_light)) {
        // Sample a random point on the plane
        // Since the plane is infinite, we need to define a sampling method
        // For example, sample within a limited area to avoid issues
        float range = 10.0f; // Define a sampling range
        float u = (gen_local.next_float() * 2.0f - 1.0f) * range;
        float v = (gen_local.next_float() * 2.0f - 1.0f) * range;
        // Find a point on the plane using two orthogonal vectors
        glm::vec3 a = (fabs(plane->m_Normal.x) > 0.1f) ? glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
        glm::vec3 u_dir = glm::normalize(glm::cross(a, plane->m_Normal));
        glm::vec3 v_dir = glm::normalize(glm::cross(plane->m_Normal, u_dir));
        sampled_point = plane->m_Point + u * u_dir + v * v_dir;
        sampled_normal = plane->m_Normal;
        emission = plane->m_Material->m_Emission;
    }
    else {
        return {};
    }

    glm::vec3 w = sampled_normal;
    glm::vec3 a = (fabs(w.x) > 0.1f) ? glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
    glm::vec3 v = glm::normalize(glm::cross(a, w));
    glm::vec3 u = glm::cross(w, v);

    glm::vec3 direction = cosine_weighted_sample(w, u, v, gen_local);

    glm::vec3 throughput = emission;

    Ray current_ray(sampled_point, direction);

    std::vector<PathVertex> path;
    int depth = 0;

    while (depth < m_Params.max_depth) {
        HitRecord rec;
        if (scene.Hit(current_ray, EPSILON, std::numeric_limits<float>::max(), rec)) {
            PathVertex vertex;
            vertex.position = rec.m_Point;
            vertex.normal = rec.m_Normal;
            vertex.material = rec.m_Material;
            vertex.emission = rec.m_Material->m_Emission;
            vertex.throughput = throughput;
            path.push_back(vertex);

            if (rec.m_Material->m_IsEmissive && depth > 0) {
                break;
            }

            glm::vec3 w_dir = rec.m_Normal;
            glm::vec3 a_dir = (fabs(w_dir.x) > 0.1f) ? glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
            glm::vec3 v_dir = glm::normalize(glm::cross(a_dir, w_dir));
            glm::vec3 u_dir = glm::cross(w_dir, v_dir);

            glm::vec3 new_direction = cosine_weighted_sample(w_dir, u_dir, v_dir, gen_local);
            float cos_theta = glm::dot(rec.m_Normal, new_direction);
            glm::vec3 brdf = rec.m_Material->m_Albedo; // Remove division by π

            throughput *= brdf * cos_theta;
            current_ray = Ray(rec.m_Point + EPSILON * new_direction, new_direction);

            depth++;
        }
        else {
            break;
        }
    }

    return path;
}

// Connect paths between eye and light vertices
glm::vec3 Renderer::connect_paths(const PathVertex& eye_vertex, const PathVertex& light_vertex, const Scene& scene) {
    glm::vec3 direction = light_vertex.position - eye_vertex.position;
    float distance_squared = glm::dot(direction, direction);
    float distance = std::sqrt(distance_squared);
    direction = glm::normalize(direction);

    Ray connecting_ray(eye_vertex.position, direction);

    HitRecord rec;
    if (scene.Hit(connecting_ray, EPSILON, distance, rec)) {
        return glm::vec3(0.0f);
    }

    // Calculate BRDF for eye vertex
    glm::vec3 f_r = eye_vertex.material->m_Albedo / glm::pi<float>();

    // Calculate BRDF for light vertex
    glm::vec3 f_l = light_vertex.material->m_Albedo / glm::pi<float>();

    float cos_theta_eye = glm::dot(eye_vertex.normal, direction);
    float cos_theta_light = glm::dot(light_vertex.normal, -direction);

    // Check if the angles are valid
    if (cos_theta_eye <= 0.0f || cos_theta_light <= 0.0f) {
        return glm::vec3(0.0f);
    }

    // Calculate the contribution
    glm::vec3 contribution = f_r * f_l * cos_theta_eye * cos_theta_light / distance_squared;

    return contribution;
}

// Double trace function
glm::vec3 Renderer::double_trace(const Ray& ray, const Scene& scene, XorShift& gen_local)
{
    glm::vec3 radiance(0.0f);

    // Trace a path from the camera
    std::vector<PathVertex> eye_path = trace_eye_path(ray, scene, gen_local);

    // Add emission from emissive vertices in the eye path
    for (const auto& vertex : eye_path) {
        if (vertex.material->m_IsEmissive) {
            radiance += vertex.throughput * vertex.emission;
        }
    }

    // Trace a path from the light
    std::vector<PathVertex> light_path = trace_light_path(scene, gen_local);

    // Connect the paths
    for (size_t i = 0; i < eye_path.size(); ++i) {
        for (size_t j = 0; j < light_path.size(); ++j) {
            glm::vec3 contrib = connect_paths(eye_path[i], light_path[j], scene);

            // Multiply by the throughput of both paths
            contrib *= eye_path[i].throughput;
            contrib *= light_path[j].throughput;

            radiance += contrib;
        }
    }

    return radiance;
}

// Compute the squared distance between two colors
inline float compute_color_distance_sq(const glm::vec3& a, const glm::vec3& b) {
    glm::vec3 diff = a - b;
    return glm::dot(diff, diff);
}

/*void bilateral_denoise(...) { ... }*/

// Generate a unique filename
std::string Renderer::generate_unique_filename(const std::string& base, const std::string& ext) {
    std::string filename;
    int file_index = 1;
    do {
        std::ostringstream oss;
        oss << base << "_" << std::setw(3) << std::setfill('0') << file_index++ << "." << ext;
        filename = oss.str();
    } while (std::filesystem::exists(filename)); // Check if the file already exists
    return filename;
}
