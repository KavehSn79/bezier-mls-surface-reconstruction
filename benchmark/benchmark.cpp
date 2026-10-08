#include "benchmark.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <queue>
#include <random>

void runBenchmark(std::vector<Point> const& points,
                  int bucketSize, int maxDepth) {
    if (points.empty()) return;

    SpatialDataStructure brute(points);
    KdTree kd(points, bucketSize, maxDepth);

    std::cout << "Bucket size:      " << bucketSize << "\n";
    std::cout << "Maximum depth:    " << maxDepth << "\n";

    constexpr int NUM_QUERIES = 1000;
    constexpr float RADIUS = 0.05f;
    constexpr unsigned int K = 10;

    std::mt19937 rng(42);
    std::uniform_int_distribution<std::size_t> dist(0, points.size() - 1);

    std::vector<Point> queries;
    queries.reserve(NUM_QUERIES);

    for (int i = 0; i < NUM_QUERIES; ++i) {
        queries.push_back(points[dist(rng)]);
    }

    std::size_t bruteRadiusCount = 0;
    std::size_t kdRadiusCount = 0;
    bool radiusCorrect = true;

    auto startRadius = std::chrono::high_resolution_clock::now();

    for (Point const &q : queries) {
        auto bruteResult = brute.collectInRadius(q, RADIUS);
        auto kdResult = kd.collectInRadius(q, RADIUS);

        std::sort(bruteResult.begin(), bruteResult.end());
        std::sort(kdResult.begin(), kdResult.end());

        if (bruteResult != kdResult) {
            radiusCorrect = false;
            std::cout << "Radius search mismatch!" << std::endl;
        }

        bruteRadiusCount += bruteResult.size();
        kdRadiusCount += kdResult.size();
    }

    auto endRadius = std::chrono::high_resolution_clock::now();

    std::size_t bruteKnnCount = 0;
    std::size_t kdKnnCount = 0;
    bool knnCorrect = true;

    auto startKnn = std::chrono::high_resolution_clock::now();

    for (Point const &q : queries) {
        auto bruteResult = brute.collectKNearest(q, K);
        auto kdResult = kd.collectKNearest(q, K);

        auto getSortedDistances = [&](std::vector<std::size_t> const &result) {
            std::vector<float> distances;
            distances.reserve(result.size());

            for (std::size_t index : result) {
                distances.push_back(
                    EuclideanDistance::measure(q, points[index]));
            }

            std::sort(distances.begin(), distances.end());
            return distances;
        };

        auto bruteDistances = getSortedDistances(bruteResult);
        auto kdDistances = getSortedDistances(kdResult);

        if (bruteDistances.size() != kdDistances.size()) {
            knnCorrect = false;
            std::cout << "kNN size mismatch!" << std::endl;
        } else {
            for (std::size_t i = 0; i < bruteDistances.size(); ++i) {
                if (std::abs(bruteDistances[i] - kdDistances[i]) > 1e-5f) {
                    knnCorrect = false;
                    std::cout << "kNN distance mismatch!" << std::endl;
                    break;
                }
            }
        }

        bruteKnnCount += bruteResult.size();
        kdKnnCount += kdResult.size();
    }

    auto endKnn = std::chrono::high_resolution_clock::now();

    auto correctnessMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                             endKnn - startRadius)
                             .count();

    // brute radius Runtime
    auto startBruteRadius = std::chrono::high_resolution_clock::now();
    for (Point const &q : queries) {
        brute.collectInRadius(q, RADIUS);
    }
    auto endBruteRadius = std::chrono::high_resolution_clock::now();

    // kd radius Runtime
    auto startKdRadius = std::chrono::high_resolution_clock::now();
    for (Point const &q : queries) {
        kd.collectInRadius(q, RADIUS);
    }
    auto endKdRadius = std::chrono::high_resolution_clock::now();

    // brute kNN Runtime
    auto startBruteKnn = std::chrono::high_resolution_clock::now();
    for (Point const &q : queries) {
        brute.collectKNearest(q, K);
    }
    auto endBruteKnn = std::chrono::high_resolution_clock::now();

    // kd kNN Runtime
    auto startKdKnn = std::chrono::high_resolution_clock::now();
    for (Point const &q : queries) {
        kd.collectKNearest(q, K);
    }
    auto endKdKnn = std::chrono::high_resolution_clock::now();

    auto bruteRadiusMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                             endBruteRadius - startBruteRadius)
                             .count();

    auto kdRadiusMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                          endKdRadius - startKdRadius)
                          .count();

    auto bruteKnnMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                          endBruteKnn - startBruteKnn)
                          .count();

    auto kdKnnMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                       endKdKnn - startKdKnn)
                       .count();

    std::cout << "\n===== Runtime Benchmark =====\n";
    std::cout << "Number of points: " << points.size() << "\n";
    std::cout << "Number of queries: " << NUM_QUERIES << "\n";
    std::cout << "Radius: " << RADIUS << "\n";
    std::cout << "k: " << K << "\n\n";

    std::cout << "Correctness check:\n";
    std::cout << "Radius search: " << (radiusCorrect ? "PASSED" : "FAILED")
              << "\n";
    std::cout << "kNN search:    " << (knnCorrect ? "PASSED" : "FAILED")
              << "\n";
    std::cout << "Correctness check time: " << correctnessMs << " ms\n\n";

    std::cout << "Radius search:\n";
    std::cout << "Brute force: " << bruteRadiusMs
              << " ms, results: " << bruteRadiusCount << "\n";
    std::cout << "KD-tree:     " << kdRadiusMs
              << " ms, results: " << kdRadiusCount << "\n\n";

    std::cout << "k-nearest search:\n";
    std::cout << "Brute force: " << bruteKnnMs
              << " ms, results: " << bruteKnnCount << "\n";
    std::cout << "KD-tree:     " << kdKnnMs << " ms, results: " << kdKnnCount
              << "\n";
    std::cout << "=============================\n\n";
}