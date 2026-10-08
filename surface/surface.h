#pragma once

#include <utility>
#include <array>
#include <vector>
#include <Eigen/Dense>
#include "../spatial/geometry.h"
#include "polyscope/curve_network.h"
#include "../spatial/spatial_data_structure.h"


struct BezierSurfaceData {
    std::vector<Point> vertices;
    std::vector<Point> tangentsU;
    std::vector<Point> tangentsV;
    std::vector<Normal> normals;
};

struct MLSSurfaceData {
    std::vector<Point> vertices;
    std::vector<Normal> normals;
};

void build_A_b(Eigen::MatrixXd &oA, Eigen::VectorXd &ob, const std::vector<Point>& iPoints);

std::vector<Point> OLS(const std::vector<Point> &points,
         std::vector<Point> vertices);

Point WLS(const KdTree &sds, Point vertex, float radius);

std::vector<Point> MLS(const KdTree &sds, std::vector<Point> vertices,
                       float radius);

Point deCasteljau1D(const std::vector<Point>& controlPoints,
    float t, Point& tangent);

Point deCasteljau2D(const std::vector<Point>& controlPoints,
    int width,
    int height,
    float u,
    float v, Point& oTangentU, Point& oTangentV);



BezierSurfaceData generateBezierSurfaceData(
    const std::vector<Point>& controlPoints,
    int width,
    int height,
    int k
);




MLSSurfaceData generateMLSSurfaceData(
    const KdTree& sds,
    const std::vector<Point>& vertices,
    float radius
);