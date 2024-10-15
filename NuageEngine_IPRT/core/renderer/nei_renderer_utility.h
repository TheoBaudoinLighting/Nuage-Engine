// nei_renderer_utility.h
#pragma once

// Class to hold renderer parameters
class RendererParameters
{
public:
    int image_width;          // Width of the rendered image
    int image_height;         // Height of the rendered image
    int samples_per_pixel;    // Number of samples per pixel
    int max_depth;            // Maximum depth for path tracing
    float russian_roulette;   // Probability for Russian roulette termination
};
