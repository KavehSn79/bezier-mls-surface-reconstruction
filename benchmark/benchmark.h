#pragma once

#include <array>
#include <vector>

#include "../spatial/spatial_data_structure.h"

void runBenchmark(std::vector<Point> const& points,
                  int bucketSize = 16,
                  int maxDepth = 64);