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

// Constructor of the Renderer class
Renderer::Renderer(const RendererParameters& params)
    : m_Params(params), m_Pixels(params.image_width * params.image_height, glm::vec3(0.0f)) {}

void Renderer::Render(const Scene& scene, const Camera& camera)
{
    // Prepare for rendering
    int image_width = m_Params.image_width;
    int image_height = m_Params.image_height;
    int samples_per_pixel = m_Params.samples_per_pixel;

    Scene scene_copy = scene; 
    scene_copy.BuildBVH();

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

            // Gamma correction
            float r_val = std::sqrt(pixel_color.r);
            float g_val = std::sqrt(pixel_color.g);
            float b_val = std::sqrt(pixel_color.b);
            m_Pixels[j * image_width + i] = glm::vec3(r_val, g_val, b_val);
        }
        // Display progress
#pragma omp critical
        {
            std::cerr << "\rProgress: " << std::fixed << std::setprecision(2) << (100.0 * (j + 1) / image_height) << "%" << std::flush;
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

// Function to generate a cosine-weighted sample
inline glm::vec3 cosine_weighted_sample(const glm::vec3& normal, XorShift& gen_local) 
{
    float r1 = gen_local.next_float();
    float r2 = gen_local.next_float();
    float phi = 2.0f * PI * r1;
    float sqrt_r2 = std::sqrt(r2);
    float x = std::cos(phi) * sqrt_r2;
    float y = std::sin(phi) * sqrt_r2;
    float z = std::sqrt(1.0f - r2);

    // Orthonormal basis
    glm::vec3 u = glm::normalize(glm::cross((fabs(normal.x) > 0.1f ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0)), normal));
    glm::vec3 v = glm::cross(normal, u);

    return u * x + v * y + normal * z;
}

// MIS weighting function
float Renderer::mis_weight(const PathVertex& eye_vertex, const PathVertex& light_vertex, const Scene& scene) 
{
    float pdf_light = compute_pdf_light(eye_vertex, light_vertex, scene);
    float pdf_bsdf = compute_pdf_bsdf(eye_vertex, light_vertex);

    // Balance heuristic for MIS, used to combine the light and BSDF PDFs
    float weight = 0.0f;
    if (pdf_light + pdf_bsdf > 0.0f) {
        weight = pdf_bsdf / (pdf_light + pdf_bsdf);
    }

    return weight;
}

// Trace a path from the camera
std::vector<PathVertex> Renderer::trace_eye_path(const Ray& ray, const Scene& scene, XorShift& gen_local) 
{
    // Trace a path from the camera
    std::vector<PathVertex> path;
    Ray current_ray = ray;
    glm::vec3 throughput(1.0f);
    int depth = 0;
    int max_depth = m_Params.max_depth;

    while (depth < max_depth) 
    {
        HitRecord rec;
        if (scene.Hit(current_ray, EPSILON, std::numeric_limits<float>::max(), rec)) {
            PathVertex vertex;
            vertex.position = rec.m_Point;
            vertex.normal = rec.m_Normal;
            vertex.material = rec.m_Material;
            vertex.emission = rec.m_Material->m_Emission;
            vertex.throughput = throughput;
            vertex.geometry = rec.m_Object;

            if (rec.m_Material->m_IsEmissive) {
                vertex.is_light = true;
                vertex.direct_light = glm::vec3(0.0f);
                path.push_back(vertex);
                break;
            }

            // Calculate direct light via NEE
            glm::vec3 direct_light = sample_lights(rec, scene, gen_local);
            vertex.direct_light = direct_light;

            path.push_back(vertex);

            // Generate a new direction
            glm::vec3 new_direction = cosine_weighted_sample(rec.m_Normal, gen_local);
            float cos_theta = glm::dot(rec.m_Normal, new_direction);
            glm::vec3 brdf = rec.m_Material->m_Albedo / PI;

            throughput *= brdf * cos_theta;

            current_ray = Ray(rec.m_Point + EPSILON * new_direction, new_direction);

            // Russian roulette
            if (depth >= 5) {
                float max_component = std::max({ throughput.r, throughput.g, throughput.b });
                float termination_probability = std::min(max_component, 0.95f);
                if (gen_local.next_float() >= termination_probability)
                    break;
                throughput /= termination_probability;
            }

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
    std::vector<HittablePtr> emissive_objects = scene.GetEmissiveObjects();

    if (emissive_objects.empty()) {
        return {};
    }

    // Select a random light source
    size_t light_index = gen_local.next() % emissive_objects.size();
    auto selected_light = emissive_objects[light_index];

    glm::vec3 sampled_point;
    glm::vec3 sampled_normal;
    glm::vec3 emission;

    // Handle different types of objects
    if (auto rect = std::dynamic_pointer_cast<Rectangle>(selected_light)) {
        sampled_point = rect->SamplePoint(gen_local);
        sampled_normal = rect->GetNormal();
        emission = rect->GetMaterial()->m_Emission;
    }
    else {
        return {};
    }

    // Generate a new direction
    glm::vec3 direction = cosine_weighted_sample(sampled_normal, gen_local);

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

            // Generate a new direction
            glm::vec3 new_direction = cosine_weighted_sample(rec.m_Normal, gen_local);
            float cos_theta = glm::dot(rec.m_Normal, new_direction);
            glm::vec3 brdf = rec.m_Material->m_Albedo / PI;

            throughput *= brdf * cos_theta;

            current_ray = Ray(rec.m_Point + EPSILON * new_direction, new_direction);

            // Russian roulette
            if (depth >= 5) {
                float max_component = std::max({ throughput.r, throughput.g, throughput.b });
                float termination_probability = std::min(max_component, 0.95f);
                if (gen_local.next_float() >= termination_probability)
                    break;
                throughput /= termination_probability;
            }

            depth++;
        }
        else {
            break;
        }
    }

    return path;
}

// Connect paths between eye and light vertices
glm::vec3 Renderer::connect_paths(const PathVertex& eye_vertex, const PathVertex& light_vertex, const Scene& scene) 
{
    glm::vec3 direction = light_vertex.position - eye_vertex.position;
    float distance_squared = glm::dot(direction, direction);
    float distance = std::sqrt(distance_squared);
    direction = glm::normalize(direction);

    Ray connecting_ray(eye_vertex.position + EPSILON * direction, direction);

    HitRecord rec;
    if (scene.Hit(connecting_ray, EPSILON, distance - EPSILON, rec)) {
        return glm::vec3(0.0f);
    }

    // Calculate the BRDF for the eye vertex
    glm::vec3 f_e = eye_vertex.material->m_Albedo / PI;

    // Calculate the BRDF for the light vertex
    glm::vec3 f_l = light_vertex.material->m_Albedo / PI;

    float cos_theta_eye = glm::dot(eye_vertex.normal, direction);
    float cos_theta_light = glm::dot(light_vertex.normal, -direction);

    // Check if angles are valid
    if (cos_theta_eye <= 0.0f || cos_theta_light <= 0.0f) {
        return glm::vec3(0.0f);
    }

    // Calculate the contribution
    glm::vec3 contribution = f_e * f_l * cos_theta_eye * cos_theta_light / distance_squared;

    // Apply throughputs
    contribution *= eye_vertex.throughput * light_vertex.throughput;

    // Calculate the MIS weight
    float weight = mis_weight(eye_vertex, light_vertex, scene);

    return contribution * weight;
}

// Double trace
glm::vec3 Renderer::double_trace(const Ray& ray, const Scene& scene, XorShift& gen_local) 
{
    glm::vec3 radiance(0.0f);

    // Trace the path from the camera
    std::vector<PathVertex> eye_path = trace_eye_path(ray, scene, gen_local);

    // Trace the path from the light
    std::vector<PathVertex> light_path = trace_light_path(scene, gen_local);

    // Connect paths with MIS
    for (const auto& eye_vertex : eye_path) 
    {
        // NEE contribution
        radiance += eye_vertex.throughput * eye_vertex.direct_light;

        for (const auto& light_vertex : light_path) {
            glm::vec3 contrib = connect_paths(eye_vertex, light_vertex, scene);

            radiance += contrib;
        }

        // Add emission if the vertex is a light source
        if (eye_vertex.is_light) {
            radiance += eye_vertex.throughput * eye_vertex.emission;
        }
    }

    return radiance;
}

// Sample lights
glm::vec3 Renderer::sample_lights(const HitRecord& rec, const Scene& scene, XorShift& gen_local) 
{
    glm::vec3 direct_light(0.0f);

    // Get the light sources in the scene
    std::vector<HittablePtr> emissive_objects = scene.GetEmissiveObjects();

    if (emissive_objects.empty()) {
        return direct_light;
    }

    // Select a random light source
    size_t light_index = gen_local.next() % emissive_objects.size();
    auto light = emissive_objects[light_index];

    // Sample a point on the light source
    glm::vec3 light_point, light_normal;
    float pdf = 1.0f;

    if (auto rect = std::dynamic_pointer_cast<Rectangle>(light)) {
        light_point = rect->SamplePoint(gen_local);
        light_normal = rect->GetNormal();
        pdf = 1.0f / rect->Area();
    } else {
        return direct_light;
    }

    // Direction to the light
    glm::vec3 to_light = light_point - rec.m_Point;
    float distance_squared = glm::dot(to_light, to_light);
    glm::vec3 direction = glm::normalize(to_light);

    // Check for shadows (shadow ray)
    Ray shadow_ray(rec.m_Point + EPSILON * direction, direction);
    HitRecord shadow_rec;
    if (scene.Hit(shadow_ray, EPSILON, std::sqrt(distance_squared) - EPSILON, shadow_rec)) {
        return direct_light;
    }

    // Calculate the contribution
    float cos_theta_surface = glm::dot(rec.m_Normal, direction);
    float cos_theta_light = glm::dot(light_normal, -direction);

    if (cos_theta_surface <= 0.0f || cos_theta_light <= 0.0f) {
        return direct_light;
    }

    MaterialPtr light_material = light->GetMaterial();
    glm::vec3 emission = light_material->m_Emission;
    glm::vec3 brdf = rec.m_Material->m_Albedo / PI;

    direct_light = emission * brdf * cos_theta_surface * cos_theta_light / (distance_squared * pdf * emissive_objects.size());

    return direct_light;
}

// Update mean and variance
void Renderer::update_statistics(glm::vec3& mean, glm::vec3& M2, int n, const glm::vec3& new_sample) 
{
    glm::vec3 delta = new_sample - mean;
    mean += delta / static_cast<float>(n);
    glm::vec3 delta2 = new_sample - mean;
    M2 += delta * delta2;
}

// Generate a unique filename
std::string Renderer::generate_unique_filename(const std::string& base, const std::string& ext) 
{
    std::string filename;
    int file_index = 1;
    do {
        std::ostringstream oss;
        oss << base << "_" << std::setw(3) << std::setfill('0') << file_index++ << "." << ext;
        filename = oss.str();
    } while (std::filesystem::exists(filename));
    return filename;
}

float Renderer::compute_pdf_light(const PathVertex& eye_vertex, const PathVertex& light_vertex, const Scene& scene) const 
{
    std::vector<HittablePtr> emissive_objects = scene.GetEmissiveObjects();
    size_t num_lights = emissive_objects.size();

    float pdf_choose_light = 1.0f / static_cast<float>(num_lights);

    float area = light_vertex.geometry->Area();
    float pdf_sample_light = 1.0f / area;

    float pdf_light = pdf_choose_light * pdf_sample_light;

    return pdf_light;
}

float Renderer::compute_pdf_bsdf(const PathVertex& eye_vertex, const PathVertex& light_vertex) const 
{
    glm::vec3 direction = glm::normalize(light_vertex.position - eye_vertex.position);

    float cos_theta = glm::dot(eye_vertex.normal, direction);
    if (cos_theta <= 0.0f) {
        return 0.0f;
    }

    float pdf_bsdf = cos_theta / PI;

    return pdf_bsdf;
}
