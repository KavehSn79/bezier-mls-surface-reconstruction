#pragma once

#include <cstddef>
#include <memory>
#include <queue>
#include <utility>
#include <vector>
#include <array>
#include "math.h"

#include "geometry.h"

class SpatialDataStructure {
  public:
    SpatialDataStructure(std::vector<Point> const &points);
    virtual ~SpatialDataStructure() = default;

    std::vector<Point> const &getPoints() const { return m_points; }

    virtual std::vector<std::size_t> collectInRadius(Point const &p,
                                                     float radius) const;
    virtual std::vector<std::size_t> collectKNearest(Point const &p,
                                                     unsigned int k) const;

  private:
    std::vector<Point> m_points;
};

class KdTree : public SpatialDataStructure {
  public:
    KdTree(std::vector<Point> const &points, int bucketSize, int maxDepth);

    std::vector<std::size_t> collectInRadius(Point const &p,
                                             float radius) const override;

    std::vector<std::size_t> collectKNearest(Point const &p,
                                             unsigned int k) const override;

    std::vector<Point>
    GetPointsFromIndices(const std::vector<std::size_t> &indices) const;

    template <typename Visitor>
    void visitInternal(std::array<float, 3> boxMin, std::array<float, 3> boxMax,
                       Visitor visitor) const {
        visitInternalRecursive(m_root.get(), 0, boxMin, boxMax, visitor);
    }

    std::vector<Point> collectInRadiusXY(Point const &p, float radius) const;

  private:
    struct Node {
        int axis = -1;
        float splitvalue = 0;
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
        std::vector<std::size_t> indices;
    };

    int m_bucketSize;
    int m_maxDepth;
    std::unique_ptr<Node> m_root;

    std::unique_ptr<Node> buildRecursive(std::vector<std::size_t> &indices,
                                         std::size_t lo, std::size_t hi,
                                         int depth);

    void collectInRadiusRecursive(Node const *node, Point const &p,
                                  float radius,
                                  std::vector<std::size_t> &result) const;

    void collectKNearestRecursive(
        Node const *node, Point const &p, unsigned int k,
        std::priority_queue<std::pair<float, std::size_t>> &nearest) const;

    template <typename Visitor>
    void visitInternalRecursive(Node const *node, int depth,
                                std::array<float, 3> boxMin,
                                std::array<float, 3> boxMax,
                                Visitor &visitor) const {
        if (!node || node->axis == -1)
            return;

        visitor(node->axis, node->splitvalue, depth, boxMin, boxMax);

        auto leftMax = boxMax;
        leftMax[node->axis] = node->splitvalue;

        auto rightMin = boxMin;
        rightMin[node->axis] = node->splitvalue;

        visitInternalRecursive(node->left.get(), depth + 1, boxMin, leftMax,
                               visitor);
        visitInternalRecursive(node->right.get(), depth + 1, rightMin, boxMax,
                               visitor);
    }

    void collectInRadiusXYRecursive(Node const *node, Point const &p,
                                    float radius,
                                    std::vector<Point> &result) const;
};