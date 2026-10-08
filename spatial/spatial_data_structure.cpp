#include "spatial_data_structure.h"

#include <algorithm>
#include <numeric>

// SpatialDataStructure

SpatialDataStructure::SpatialDataStructure(std::vector<Point> const &points)
    : m_points(points) {}

std::vector<std::size_t>
SpatialDataStructure::collectInRadius(Point const &p, float radius) const {
    std::vector<std::size_t> result;
    for (std::size_t i = 0; i < m_points.size(); ++i) {
        float distance = EuclideanDistance::measure(p, m_points[i]);
        if (distance <= radius)
            result.push_back(i);
    }
    return result;
}

std::vector<std::size_t>
SpatialDataStructure::collectKNearest(Point const &p, unsigned int k) const {
    std::vector<std::pair<float, std::size_t>> dists;
    dists.reserve(m_points.size());
    for (std::size_t i = 0; i < m_points.size(); ++i) {
        dists.emplace_back(EuclideanDistance::measure(p, m_points[i]), i);
    }

    unsigned int kk = std::min<unsigned int>(k, dists.size());
    std::partial_sort(dists.begin(), dists.begin() + kk, dists.end());

    std::vector<std::size_t> result;
    result.reserve(kk);
    for (unsigned int i = 0; i < kk; ++i)
        result.push_back(dists[i].second);
    return result;
}

// KdTree

KdTree::KdTree(std::vector<Point> const &points, int bucketSize, int maxDepth)
    : SpatialDataStructure(points),
        m_bucketSize(bucketSize),
        m_maxDepth(maxDepth) {
    if (points.empty())
        return;

    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);

    m_root = buildRecursive(indices, 0, indices.size(), 0);
}

std::vector<std::size_t> KdTree::collectInRadius(Point const &p,
                                                 float radius) const {
    std::vector<std::size_t> result;
    collectInRadiusRecursive(m_root.get(), p, radius, result);
    return result;
}

std::vector<std::size_t> KdTree::collectKNearest(Point const &p,
                                                 unsigned int k) const {
    std::vector<std::size_t> result;
    std::priority_queue<std::pair<float, std::size_t>> nearest;

    collectKNearestRecursive(m_root.get(), p, k, nearest);

    result.reserve(nearest.size());
    while (!nearest.empty()) {
        result.push_back(nearest.top().second);
        nearest.pop();
    }
    std::reverse(result.begin(), result.end());
    return result;
}

std::unique_ptr<KdTree::Node>
KdTree::buildRecursive(std::vector<std::size_t> &indices, std::size_t lo,
                       std::size_t hi, int depth) {
    if (hi - lo <= static_cast<std::size_t>(m_bucketSize) || depth >= m_maxDepth) {
        auto leaf = std::make_unique<Node>();
        leaf->axis = -1;
        leaf->indices = std::vector<std::size_t>(indices.begin() + lo,
                                                 indices.begin() + hi);
        return leaf;
    }

    int axis = depth % 3;
    std::size_t mid = lo + (hi - lo) / 2;

    std::nth_element(indices.begin() + lo, indices.begin() + mid,
                     indices.begin() + hi,
                     [&](std::size_t a, std::size_t b) {
                         return getPoints()[a][axis] < getPoints()[b][axis];
                     });

    auto inner_node = std::make_unique<Node>();
    inner_node->axis = axis;
    inner_node->splitvalue = getPoints()[indices[mid]][axis];
    inner_node->left = buildRecursive(indices, lo, mid, depth + 1);
    inner_node->right = buildRecursive(indices, mid, hi, depth + 1);
    return inner_node;
}

void KdTree::collectInRadiusRecursive(
    Node const *node, Point const &p, float radius,
    std::vector<std::size_t> &result) const {
    if (node == nullptr)
        return;

    if (node->axis == -1) {
        for (std::size_t index : node->indices) {
            float distance = EuclideanDistance::measure(p, getPoints()[index]);
            if (distance <= radius)
                result.push_back(index);
        }
        return;
    }

    int axis = node->axis;
    float diff = p[axis] - node->splitvalue;

    Node const *nearChild;
    Node const *farChild;
    if (diff <= 0) {
        nearChild = node->left.get();
        farChild = node->right.get();
    } else {
        nearChild = node->right.get();
        farChild = node->left.get();
    }

    collectInRadiusRecursive(nearChild, p, radius, result);

    if (radius >= std::abs(diff))
        collectInRadiusRecursive(farChild, p, radius, result);
}

void KdTree::collectKNearestRecursive(
    Node const *node, Point const &p, unsigned int k,
    std::priority_queue<std::pair<float, std::size_t>> &nearest) const {
    if (node == nullptr || k == 0)
        return;

    if (node->axis == -1) {
        for (std::size_t index : node->indices) {
            float distance = EuclideanDistance::measure(p, getPoints()[index]);
            if (nearest.size() < k)
                nearest.push({distance, index});
            else if (distance < nearest.top().first) {
                nearest.pop();
                nearest.push({distance, index});
            }
        }
        return;
    }

    int axis = node->axis;
    float diff = p[axis] - node->splitvalue;

    Node const *nearChild;
    Node const *farChild;
    if (diff <= 0) {
        nearChild = node->left.get();
        farChild = node->right.get();
    } else {
        nearChild = node->right.get();
        farChild = node->left.get();
    }

    collectKNearestRecursive(nearChild, p, k, nearest);

    if (nearest.size() < k || std::abs(diff) < nearest.top().first) {
        collectKNearestRecursive(farChild, p, k, nearest);
    }
}

std::vector<Point>
KdTree::GetPointsFromIndices(const std::vector<std::size_t> &indices) const {

    std::vector<Point> resultPoints;
    resultPoints.reserve(indices.size());

    for (std::size_t i = 0; i < indices.size(); ++i) {
        resultPoints.push_back(getPoints()[indices[i]]);
    }

    return resultPoints;
}

std::vector<Point> KdTree::collectInRadiusXY(Point const &p,
                                             float radius) const {
    std::vector<Point> result;
    collectInRadiusXYRecursive(m_root.get(), p, radius, result);
    return result;
}

void KdTree::collectInRadiusXYRecursive(Node const *node, Point const &p,
                                        float radius,
                                        std::vector<Point> &result) const {

    if (node == nullptr)
        return;

    if (node->axis == -1) {
        for (std::size_t index : node->indices) {
            float distance =
                EuclideanDistance::x_y_domain_measure(p, getPoints()[index]);

            if (distance <= radius)
                result.push_back(getPoints()[index]);
        }
        return;
    }

    int axis = node->axis;
    // If the split axis is z, we cannot prune using an x-y radius query.
    if (axis == 2) {
        collectInRadiusXYRecursive(node->left.get(), p, radius, result);
        collectInRadiusXYRecursive(node->right.get(), p, radius, result);
        return;
    }

    float diff = p[axis] - node->splitvalue;

    Node const *nearChild;
    Node const *farChild;

    if (diff <= 0) {
        nearChild = node->left.get();
        farChild = node->right.get();
    } else {
        nearChild = node->right.get();
        farChild = node->left.get();
    }

    collectInRadiusXYRecursive(nearChild, p, radius, result);

    if (std::abs(diff) <= radius) {
        collectInRadiusXYRecursive(farChild, p, radius, result);
    }
}