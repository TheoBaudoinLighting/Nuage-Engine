// nei_utility.cpp
#include "nei_utility.h"

std::mt19937 gen(std::random_device{}());
std::uniform_real_distribution<float> dis(0.0f, 1.0f);