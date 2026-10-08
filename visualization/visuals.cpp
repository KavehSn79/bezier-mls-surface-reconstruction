#include <algorithm>
#include <numeric>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"


#include "visuals.h"

std::array<Point, 4> Plane::corners() const {
    // The two axes that span the plane (everything except `axis`).
    int u = (axis + 1) % 3;
    int v = (axis + 2) % 3;

    Point c0{}, c1{}, c2{}, c3{};

    auto set = [&](Point &p, float uVal, float vVal) {
        p[axis] = offset;
        p[u] = uVal;
        p[v] = vVal;
    };

    set(c0, boxMin[u], boxMin[v]);
    set(c1, boxMax[u], boxMin[v]);
    set(c2, boxMax[u], boxMax[v]);
    set(c3, boxMin[u], boxMax[v]);

    return {c0, c1, c2, c3};
}


void computeAABB(std::vector<Point> const& points,
                 std::array<float, 3>& boxMin,
                 std::array<float, 3>& boxMax) {
    if (points.empty()) {
        boxMin = {0, 0, 0};
        boxMax = {0, 0, 0};
        return;
    }
    boxMin = boxMax = points[0];
    for (auto const& p : points) {
        for (int a = 0; a < 3; ++a) {
            if (p[a] < boxMin[a]) boxMin[a] = p[a];
            if (p[a] > boxMax[a]) boxMax[a] = p[a];
        }
    }
}

polyscope::CurveNetwork *renderGrid(const std::string &name, int m, int n,
                                    const std::vector<Point> &points,
                                    std::vector<Point> &vertices,
                                    std::vector<std::array<size_t, 2>> &edges) {
    vertices.clear();
    edges.clear();

    std::array<float, 3> boxMin;
    std::array<float, 3> boxMax;

    computeAABB(points, boxMin, boxMax);

    float x_min = boxMin[0];
    float y_min = boxMin[1];
    float x_max = boxMax[0];
    float y_max = boxMax[1];

    float diff_x = (x_max - x_min) / m;
    float diff_y = (y_max - y_min) / n;

    for (int y = 0; y <= n; y++) {
        for (int x = 0; x <= m; x++) {
            vertices.push_back({x_min + x * diff_x, y_min + y * diff_y, 0.0f});
        }
    }

    auto get_vertex_idx = [&](int x, int y) -> size_t {
        return y * (m + 1) + x;
    };

    for (int y = 0; y <= n; y++) {
        for (int x = 0; x < m; x++) {
            edges.push_back({get_vertex_idx(x, y), get_vertex_idx(x + 1, y)});
        }
    }

    for (int y = 0; y < n; y++) {
        for (int x = 0; x <= m; x++) {
            edges.push_back({get_vertex_idx(x, y), get_vertex_idx(x, y + 1)});
        }
    }

    return polyscope::registerCurveNetwork(name, vertices, edges);
}

std::vector<DepthPlane> generateKdTreePlanes(KdTree const& tree,
                                             std::array<float, 3> const& boxMin,
                                             std::array<float, 3> const& boxMax) {
    std::vector<DepthPlane> out;
    tree.visitInternal(boxMin, boxMax,
        [&](int axis, float splitvalue, int depth,
            std::array<float, 3> nodeMin, std::array<float, 3> nodeMax) {
            Plane p{axis, splitvalue, nodeMin, nodeMax};
            out.push_back({p, depth});
        });
    return out;
}

void renderPlanes(std::string const& name,
                  std::vector<DepthPlane> const& planes,
                  bool fixedDepth, int depth) {
    std::vector<std::array<float, 3>> vertices;
    std::vector<std::array<size_t, 4>> faces;

    vertices.reserve(planes.size() * 4);
    faces.reserve(planes.size());

    for (auto const& dp : planes) {
        if (fixedDepth && dp.depth != depth) continue;
        auto c = dp.plane.corners();
        size_t base = vertices.size();
        vertices.push_back(c[0]);
        vertices.push_back(c[1]);
        vertices.push_back(c[2]);
        vertices.push_back(c[3]);
        faces.push_back({base, base + 1, base + 2, base + 3});
    }

    polyscope::registerSurfaceMesh(name, vertices, faces);
}

void clearPlanes(std::string const& name) {
    if (polyscope::hasSurfaceMesh(name)) {
        polyscope::removeSurfaceMesh(name);
    }
}


void renderQueryHighlight(std::string const& pointCloudName,
                          std::size_t numPoints,
                          std::vector<std::size_t> const& highlightedIndices) {
    auto* cloud = polyscope::getPointCloud(pointCloudName);
    if (cloud == nullptr) return;

    std::vector<std::array<float, 3>> colors(numPoints, {0.7f, 0.7f, 0.7f}); // base grey
    for (std::size_t idx : highlightedIndices) {
        if (idx < numPoints) colors[idx] = {1.0f, 0.2f, 0.2f}; // red
    }

    cloud->addColorQuantity("query_highlight", colors)->setEnabled(true);
}

void clearQueryHighlight(std::string const& pointCloudName) {
    auto* cloud = polyscope::getPointCloud(pointCloudName);
    if (cloud == nullptr) return;
    cloud->removeQuantity("query_highlight");
}