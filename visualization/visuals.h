#pragma once

#include <cstddef>
#include <memory>
#include <queue>
#include <utility>
#include <array>
#include <vector>

#include "../spatial/spatial_data_structure.h"
#include "../spatial/geometry.h"  // for Point
#include "polyscope/curve_network.h"



// Axis-aligned cutting plane.
// Positioned at `offset` along that axis.
// Sized to span the bounding box.
struct Plane {
    int axis;                            // 0 = x, 1 = y, 2 = z
    float offset;                        // coordinate along `axis`
    std::array<float, 3> boxMin;         // min corner of the point-cloud
    std::array<float, 3> boxMax;         // max corner of the point-cloud

    std::array<Point, 4> corners() const; // Four corners of the rectangle this plane carves out
};

struct DepthPlane {
    Plane plane;
    int depth;
};

void computeAABB(std::vector<Point> const& points,
                 std::array<float, 3>& boxMin,
                 std::array<float, 3>& boxMax);

std::vector<DepthPlane> generateKdTreePlanes(KdTree const& tree,
                                             std::array<float, 3> const& boxMin,
                                             std::array<float, 3> const& boxMax);

void renderPlanes(std::string const& name,
                  std::vector<DepthPlane> const& planes,
                  bool fixedDepth, int depth);

void clearPlanes(std::string const& name);

void renderQueryHighlight(std::string const& pointCloudName,
                          std::size_t numPoints,
                          std::vector<std::size_t> const& highlightedIndices);

void clearQueryHighlight(std::string const& pointCloudName);

polyscope::CurveNetwork *renderGrid(const std::string &name, int m, int n,
                                    const std::vector<Point> &points,
                                    std::vector<Point> &vertices,
                                    std::vector<std::array<size_t, 2>> &edges);