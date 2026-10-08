#include "surface.h"



void build_A_b(Eigen::MatrixXd& oA, Eigen::VectorXd& ob,
    const std::vector<Point>& iPoints) {
    int numRows = iPoints.size();

    oA.resize(numRows, 6);
    ob.resize(numRows);

    for (int i = 0; i < numRows; ++i) {
        float x_i = iPoints[i][0];
        float y_i = iPoints[i][1];
        float z_i = iPoints[i][2];

        oA(i, 0) = 1.0f;
        oA(i, 1) = x_i;
        oA(i, 2) = y_i;
        oA(i, 3) = x_i * x_i;
        oA(i, 4) = x_i * y_i;
        oA(i, 5) = y_i * y_i;

        ob(i) = z_i;
    }
}





std::vector<Point> OLS(const std::vector<Point> &points,
                       std::vector<Point> vertices) {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    build_A_b(A, b, points);

    Eigen::MatrixXd M = A.transpose() * A;
    Eigen::VectorXd r = A.transpose() * b;

    Eigen::VectorXd coeffs = M.ldlt().solve(r);

    for (size_t i = 0; i < vertices.size(); ++i) {
        float x_i = vertices[i][0];
        float y_i = vertices[i][1];
        float z_i = coeffs[0] + coeffs[1] * x_i + coeffs[2] * y_i +
                    coeffs[3] * x_i * x_i + coeffs[4] * x_i * y_i +
                    coeffs[5] * y_i * y_i;
        vertices[i][2] = z_i;
    }
    return vertices;
}

Point WLS(const KdTree &sds, Point vertex, float radius) {
    std::vector<Point> nearPointsToGrid = sds.collectInRadiusXY(vertex, radius);
    
    size_t N = nearPointsToGrid.size();
    if (N < 6)
        return vertex;

    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    build_A_b(A, b, nearPointsToGrid);
    
    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(N, N);
    for (size_t i = 0; i < N; ++i) {
        float distance =
            EuclideanDistance::x_y_domain_measure(vertex, nearPointsToGrid[i]);
        float q = distance / radius;
        if (q > 1)
            W(i, i) = 0;
        else {
            W(i, i) = pow((1 - q), 4) * (4 * q + 1);
        }
    }
    Eigen::MatrixXd M = A.transpose() * W * A;
    Eigen::VectorXd r = A.transpose() * W * b;

    Eigen::VectorXd coeffs = M.ldlt().solve(r);
    float x_i = vertex[0];
    float y_i = vertex[1];
    float z_i = coeffs[0] + coeffs[1] * x_i + coeffs[2] * y_i +
                coeffs[3] * x_i * x_i + coeffs[4] * x_i * y_i +
                coeffs[5] * y_i * y_i;
    vertex[2] = z_i;
    return vertex;
}

std::vector<Point> MLS(const KdTree& sds, std::vector<Point> vertices,
    float radius) {
    std::vector<Point> surfaceVertices;
    surfaceVertices.reserve(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        Point vertex = WLS(sds, vertices[i], radius);
        surfaceVertices.push_back(vertex);
    }
    return surfaceVertices;
}

Point deCasteljau1D(const std::vector<Point>& controlPoints, float t, Point& oTangent)
{
    if (controlPoints.empty()) {
        return {0.0f, 0.0f, 0.0f};
    }
    
    std::vector<Point> points = controlPoints;

    for (size_t level = 1; level < controlPoints.size(); ++level) {
        if (controlPoints.size() - level == 1) {
            Point p0 = points[0];
            Point p1 = points[1];

            oTangent = {
                p1[0] - p0[0],
                p1[1] - p0[1],
                p1[2] - p0[2]
            };
        }
        for (size_t i = 0; i < controlPoints.size() - level; ++i) {
            Point p0 = points[i];
            Point p1 = points[i + 1];
            points[i] = {
                (1.0f - t)* p0[0] + t * p1[0],
                (1.0f - t)* p0[1] + t * p1[1],
                (1.0f - t)* p0[2] + t * p1[2]
            };
        }
    }

    return points[0];
}

Point deCasteljau2D(const std::vector<Point>& controlPoints,
    int width,
    int height,
    float u,
    float v,
    Point& oTangentU,
    Point& oTangentV) {
    std::vector<Point> intermediatePoints;
    std::vector<Point> intermediateTangentsU;

    intermediatePoints.reserve(height);
    intermediateTangentsU.reserve(height);

    for (int row = 0; row < height; ++row) {
        std::vector<Point> rowPoints;
        rowPoints.reserve(width);

        for (int col = 0; col < width; ++col) {
            int idx = row * width + col;
            rowPoints.push_back(controlPoints[idx]);
        }

        Point rowTangentU;
        Point q = deCasteljau1D(rowPoints, u, rowTangentU);

        intermediatePoints.push_back(q);
        intermediateTangentsU.push_back(rowTangentU);
    }

    Point dummy;
    oTangentU = deCasteljau1D(intermediateTangentsU, v, dummy);

    return deCasteljau1D(intermediatePoints, v, oTangentV);
}

BezierSurfaceData generateBezierSurfaceData(
    const std::vector<Point>& controlPoints,
    int width,
    int height,
    int k
) {
    int controlWidth = width + 1;
    int controlHeight = height + 1;

    int sampleWidth = k * width;
    int sampleHeight = k * height;

    BezierSurfaceData data;

    size_t numVertices =
        (sampleWidth + 1) * (sampleHeight + 1);

    data.vertices.reserve(numVertices);
    data.tangentsU.reserve(numVertices);
    data.tangentsV.reserve(numVertices);
    data.normals.reserve(numVertices);

    for (int y = 0; y <= sampleHeight; ++y) {

        float v =
            static_cast<float>(y) /
            static_cast<float>(sampleHeight);

        for (int x = 0; x <= sampleWidth; ++x) {

            float u =
                static_cast<float>(x) /
                static_cast<float>(sampleWidth);

            Point tangentU;
            Point tangentV;

            Point p = deCasteljau2D(
                controlPoints,
                controlWidth,
                controlHeight,
                u,
                v,
                tangentU,
                tangentV
            );

            Normal n = {
                tangentU[1] * tangentV[2] -
                    tangentU[2] * tangentV[1],

                tangentU[2] * tangentV[0] -
                    tangentU[0] * tangentV[2],

                tangentU[0] * tangentV[1] -
                    tangentU[1] * tangentV[0]
            };

            float len = std::sqrt(
                n[0] * n[0] +
                n[1] * n[1] +
                n[2] * n[2]
            );

            if (len > 1e-8f) {
                n[0] /= len;
                n[1] /= len;
                n[2] /= len;
            }

            data.vertices.push_back(p);
            data.tangentsU.push_back(tangentU);
            data.tangentsV.push_back(tangentV);
            data.normals.push_back(n);
        }
    }

    return data;
}


MLSSurfaceData generateMLSSurfaceData(
    const KdTree& sds,
    const std::vector<Point>& vertices,
    float radius
) {
    MLSSurfaceData data;

    data.vertices.reserve(vertices.size());
    data.normals.reserve(vertices.size());

    for (size_t i = 0; i < vertices.size(); ++i) {
        Point vertex = vertices[i];

        std::vector<Point> nearPointsToGrid =
            sds.collectInRadiusXY(vertex, radius);

        size_t N = nearPointsToGrid.size();

        if (N < 6) {
            data.vertices.push_back(vertex);
            data.normals.push_back({ 0.0f, 0.0f, 1.0f });
            continue;
        }

        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        build_A_b(A, b, nearPointsToGrid);

        Eigen::MatrixXd W = Eigen::MatrixXd::Zero(N, N);

        for (size_t j = 0; j < N; ++j) {
            float distance =
                EuclideanDistance::x_y_domain_measure(vertex, nearPointsToGrid[j]);

            float q = distance / radius;

            if (q > 1.0f) {
                W(j, j) = 0.0f;
            }
            else {
                W(j, j) = std::pow(1.0f - q, 4) * (4.0f * q + 1.0f);
            }
        }

        Eigen::MatrixXd M = A.transpose() * W * A;
        Eigen::VectorXd r = A.transpose() * W * b;

        Eigen::VectorXd coeffs = M.ldlt().solve(r);

        float x = vertex[0];
        float y = vertex[1];

        float z =
            coeffs[0] +
            coeffs[1] * x +
            coeffs[2] * y +
            coeffs[3] * x * x +
            coeffs[4] * x * y +
            coeffs[5] * y * y;

        vertex[2] = z;

        float dzdx =
            coeffs[1] +
            2.0f * coeffs[3] * x +
            coeffs[4] * y;

        float dzdy =
            coeffs[2] +
            coeffs[4] * x +
            2.0f * coeffs[5] * y;

        Point tangentX = { 1.0f, 0.0f, dzdx };
        Point tangentY = { 0.0f, 1.0f, dzdy };

        Normal n = {
            tangentX[1] * tangentY[2] - tangentX[2] * tangentY[1],
            tangentX[2] * tangentY[0] - tangentX[0] * tangentY[2],
            tangentX[0] * tangentY[1] - tangentX[1] * tangentY[0]
        };

        float len = std::sqrt(
            n[0] * n[0] +
            n[1] * n[1] +
            n[2] * n[2]
        );

        if (len > 1e-8f) {
            n[0] /= len;
            n[1] /= len;
            n[2] /= len;
        }
        else {
            n = { 0.0f, 0.0f, 1.0f };
        }

        data.vertices.push_back(vertex);
        data.normals.push_back(n);
    }

    return data;
}