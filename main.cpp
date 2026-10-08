#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "args/args.hxx"

#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/pick.h"
#include "imgui.h"
#include "polyscope/surface_mesh.h"

#include "portable-file-dialogs.h"

#include "benchmark/benchmark.h"
#include "spatial/geometry.h"
#include "spatial/spatial_data_structure.h"

#include "visualization/visuals.h"
#include "surface/surface.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/vector_quantity.h"

std::vector<DepthPlane> gPlanes;


polyscope::PointCloud *pc = nullptr;
std::unique_ptr<KdTree> sds;

bool gRunBenchmark = false;

static bool nextDataLine(std::ifstream &in, std::string &line) {
    while (std::getline(in, line)) {
        std::size_t start = line.find_first_not_of(" \t\r");

        if (start == std::string::npos)
            continue;

        if (line[start] == '#')
            continue;

        line = line.substr(start);
        return true;
    }

    return false;
}

void readOff(const std::string &filename, std::vector<Point> &points,
             std::vector<Normal> &normals) {
    points.clear();
    normals.clear();

    std::ifstream in(filename);

    if (!in) {
        polyscope::warning("Failed to open file: " + filename);
        return;
    }

    std::string line;

    if (!nextDataLine(in, line)) {
        polyscope::warning("Empty OFF file: " + filename);
        return;
    }

    bool hasNormals = false;

    if (!std::isdigit(static_cast<unsigned char>(line[0]))) {
        std::istringstream iss(line);

        std::string keyword;
        iss >> keyword;

        if (keyword == "OFF") {
            hasNormals = false;
        } else if (keyword == "NOFF") {
            hasNormals = true;
        } else {
            polyscope::warning("Unsupported OFF variant: " + keyword);
            return;
        }

        if (!nextDataLine(in, line)) {
            polyscope::warning("Missing counts line in " + filename);
            return;
        }
    }

    int numVertices = 0;
    int numFaces = 0;
    int numEdges = 0;

    {
        std::istringstream iss(line);

        if (!(iss >> numVertices >> numFaces >> numEdges)) {
            polyscope::warning("Malformed counts line in: " + filename);
            return;
        }
    }

    points.reserve(numVertices);

    if (hasNormals)
        normals.reserve(numVertices);

    for (int i = 0; i < numVertices; ++i) {
        float x, y, z;

        if (!(in >> x >> y >> z)) {
            polyscope::warning("Unexpected EOF reading vertex " +
                               std::to_string(i));
            return;
        }

        points.push_back({x, y, z});

        if (hasNormals) {
            float nx, ny, nz;

            if (!(in >> nx >> ny >> nz)) {
                polyscope::warning("Unexpected EOF reading normal " +
                                   std::to_string(i));
                return;
            }

            normals.push_back({nx, ny, nz});
        }
    }
}

void readOff(const std::string &filename, std::vector<Point> &points) {
    std::vector<Normal> dummyNormals;
    readOff(filename, points, dummyNormals);
}

void loadMesh() {
    auto paths =
        pfd::open_file("Load Off", "",
                       std::vector<std::string>{"point data (*.off)", "*.off"},
                       pfd::opt::none)
            .result();

    if (paths.empty())
        return;

    std::filesystem::path path(paths[0]);

    if (path.extension() != ".off")
        return;

    std::cout << "Loading: " << path.string() << std::endl;

    std::vector<Point> points;
    readOff(path.string(), points);

    std::cout << "Parsed " << points.size() << " points" << std::endl;

    if (pc != nullptr) {
        polyscope::removeStructure(pc);
        pc = nullptr;
    }

    pc = polyscope::registerPointCloud("Points", points);
    clearQueryHighlight("Points");

    sds = std::make_unique<KdTree>(points, 16, 20);

    std::array<float, 3> bmin, bmax;
    computeAABB(points, bmin, bmax);
    gPlanes = generateKdTreePlanes(*sds, bmin, bmax);

    std::cout << "kd-tree built" << std::endl;

    if (gRunBenchmark) {
        runBenchmark(points, 16, 20);
    }
}

void callback() {
    static bool showPoints = true;

    static int width = 10;
    static int height = 10;
    static float radius = 0.100f;
    static int subdivision = 4;

    static bool showGrid = false;
    static bool showControlMesh = false;

    static bool showSurfaceMesh = false;
    static bool showNormals = false;
    static bool showTangents = false;

    static bool showMLSSurfaceMesh = false;
    static bool showMLSNormals = false;
    static int MLSSubdivision = 2;
    static float MLSradius = 0.100f;

    static polyscope::CurveNetwork* grid = nullptr;
    static polyscope::CurveNetwork* controlMesh = nullptr;
    static polyscope::PointCloud* surfacePoints = nullptr;
    static polyscope::SurfaceMesh* bezierMesh = nullptr;

    static std::vector<Point> gridVertices;
    static std::vector<std::array<size_t, 2>> gridEdges;
    static std::vector<Point> controlVertices;

    static std::vector<Point> highResGridVertices;
    static polyscope::SurfaceMesh* mlsMesh = nullptr;

    static int lastGridWidth = -1;
    static int lastGridHeight = -1;

    static int lastControlWidth = -1;
    static int lastControlHeight = -1;
    static float lastControlRadius = -1.0f;

    static int lastBezierWidth = -1;
    static int lastBezierHeight = -1;
    static float lastBezierRadius = -1.0f;
    static int lastBezierSubdivision = -1;

    static bool lastShowNormals = false;
    static bool lastShowTangents = false;

    static int lastMLSWidth = -1;
    static int lastMLSHeight = -1;
    static float lastMLSRadius = 1.0f;
    static int lastMLSSubdivision = -1;

    static bool lastShowMLSNormals = false;
    static bool lastShowMLSTangents = false;


    ImGui::Text("Input Data");
    ImGui::Separator();

    if (ImGui::Button("Load Points")) {
        loadMesh();

        if (grid != nullptr) {
            polyscope::removeStructure(grid);
            grid = nullptr;
        }

        if (controlMesh != nullptr) {
            polyscope::removeStructure(controlMesh);
            controlMesh = nullptr;
        }

        if (surfacePoints != nullptr) {
            polyscope::removeStructure(surfacePoints);
            surfacePoints = nullptr;
        }

        if (bezierMesh != nullptr) {
            polyscope::removeStructure(bezierMesh);
            bezierMesh = nullptr;
        }

        if (mlsMesh != nullptr) {
            polyscope::removeStructure(mlsMesh);
            mlsMesh = nullptr;
        }

        gridVertices.clear();
        gridEdges.clear();
        controlVertices.clear();

        lastGridWidth = -1;
        lastGridHeight = -1;

        lastControlWidth = -1;
        lastControlHeight = -1;
        lastControlRadius = -1.0f;

        lastBezierWidth = -1;
        lastBezierHeight = -1;
        lastBezierRadius = -1.0f;
        lastBezierSubdivision = -1;

        lastShowNormals = false;
        lastShowTangents = false;

        lastMLSWidth = -1;
        lastMLSHeight = -1;
        lastMLSRadius = -1.0f;
        lastMLSSubdivision = -1;
        lastShowMLSNormals = false;
    }

    ImGui::Checkbox("Show Points", &showPoints);

    if (pc != nullptr) {
        pc->setEnabled(showPoints);
    }

    ImGui::Text("Control Points");
    ImGui::Separator();

    ImGui::Text("Grid");
    ImGui::SameLine();

    ImGui::SetNextItemWidth(150);
    ImGui::SliderInt("##width", &width, 1, 100);

    ImGui::SameLine();

    ImGui::SetNextItemWidth(150);
    ImGui::SliderInt("##height", &height, 1, 100);

    ImGui::SameLine();
    ImGui::Text("Width, Height");

    ImGui::SetNextItemWidth(225);
    ImGui::SliderFloat("##radius", &radius, 0.001f, 1.0f, "%.3f");

    ImGui::SameLine();
    ImGui::Text("Radius");

    ImGui::Checkbox("Show Grid", &showGrid);

    auto ensureFlatGrid = [&]() {
        if (sds == nullptr) {
            return;
        }

        if (grid == nullptr ||
            lastGridWidth != width ||
            lastGridHeight != height) {

            std::vector<Point> points = sds->getPoints();

            if (grid != nullptr) {
                polyscope::removeStructure(grid);
                grid = nullptr;
            }

            grid = renderGrid(
                "Flat Grid",
                width,
                height,
                points,
                gridVertices,
                gridEdges
            );

            grid->setRadius(0.001);
            grid->setColor(glm::vec3(1.f, 0.f, 0.f));
            grid->setEnabled(showGrid);

            lastGridWidth = width;
            lastGridHeight = height;

            if (controlMesh != nullptr) {
                polyscope::removeStructure(controlMesh);
                controlMesh = nullptr;
            }

            if (surfacePoints != nullptr) {
                polyscope::removeStructure(surfacePoints);
                surfacePoints = nullptr;
            }

            if (bezierMesh != nullptr) {
                polyscope::removeStructure(bezierMesh);
                bezierMesh = nullptr;
            }

            controlVertices.clear();

            lastControlWidth = -1;
            lastControlHeight = -1;
            lastControlRadius = -1.0f;

            lastBezierWidth = -1;
            lastBezierHeight = -1;
            lastBezierRadius = -1.0f;
            lastBezierSubdivision = -1;


        }
        };

    ensureFlatGrid();

    if (grid != nullptr) {
        grid->setEnabled(showGrid);
    }

    ImGui::Checkbox("Show Control Mesh", &showControlMesh);

    auto ensureControlMesh = [&]() {
        if (sds == nullptr) {
            return;
        }

        ensureFlatGrid();

        if (controlMesh == nullptr ||
            surfacePoints == nullptr ||
            lastControlWidth != width ||
            lastControlHeight != height ||
            lastControlRadius != radius) {

            if (controlMesh != nullptr) {
                polyscope::removeStructure(controlMesh);
                controlMesh = nullptr;
            }

            if (surfacePoints != nullptr) {
                polyscope::removeStructure(surfacePoints);
                surfacePoints = nullptr;
            }

            controlVertices = MLS(*sds, gridVertices, radius);

            controlMesh = polyscope::registerCurveNetwork(
                "MLS Control Mesh",
                controlVertices,
                gridEdges
            );

            controlMesh->setRadius(0.001);
            controlMesh->setColor(glm::vec3(0.f, 1.f, 0.f));
            controlMesh->setEnabled(showControlMesh);

            surfacePoints = polyscope::registerPointCloud(
                "MLS Vertices",
                controlVertices
            );

            surfacePoints->setPointRadius(0.005);
            surfacePoints->setPointColor(glm::vec3(0.f, 1.f, 0.f));
            surfacePoints->setEnabled(showControlMesh);

            lastControlWidth = width;
            lastControlHeight = height;
            lastControlRadius = radius;

            if (bezierMesh != nullptr) {
                polyscope::removeStructure(bezierMesh);
                bezierMesh = nullptr;
            }

            lastBezierWidth = -1;
            lastBezierHeight = -1;
            lastBezierRadius = -1.0f;
            lastBezierSubdivision = -1;
        }
        };

    ensureControlMesh();

    if (controlMesh != nullptr) {
        controlMesh->setEnabled(showControlMesh);
    }

    if (surfacePoints != nullptr) {
        surfacePoints->setEnabled(showControlMesh);
    }

    ImGui::Text("Bezier Surface");
    ImGui::Separator();

    ImGui::Checkbox("Show Surface Mesh", &showSurfaceMesh);
    ImGui::Checkbox("Show Normals", &showNormals);
    ImGui::Checkbox("Show Tangents", &showTangents);

    ImGui::SetNextItemWidth(225);
    ImGui::SliderInt("##bezierSubdivision", &subdivision, 1, 20);

    ImGui::SameLine();
    ImGui::Text("Subdivision");

    if (showSurfaceMesh && sds != nullptr) {
        ensureControlMesh();

        if (bezierMesh == nullptr ||
            lastBezierWidth != width ||
            lastBezierHeight != height ||
            lastBezierRadius != radius ||
            lastBezierSubdivision != subdivision ||
            lastShowNormals != showNormals ||
            lastShowTangents != showTangents) {

            if (bezierMesh != nullptr) {
                polyscope::removeStructure(bezierMesh);
                bezierMesh = nullptr;
            }

            BezierSurfaceData bezierData =
                generateBezierSurfaceData(
                    controlVertices,
                    width,
                    height,
                    subdivision
                );

            int sampleWidth = subdivision * width;
            int sampleHeight = subdivision * height;

            std::vector<std::array<size_t, 4>> faces;
            faces.reserve(sampleWidth * sampleHeight);

            auto idx = [&](int x, int y) -> size_t {
                return y * (sampleWidth + 1) + x;
                };

            for (int y = 0; y < sampleHeight; ++y) {
                for (int x = 0; x < sampleWidth; ++x) {
                    faces.push_back({
                        idx(x, y),
                        idx(x + 1, y),
                        idx(x + 1, y + 1),
                        idx(x, y + 1)
                        });
                }
            }

            bezierMesh = polyscope::registerSurfaceMesh(
                "Bezier Surface",
                bezierData.vertices,
                faces
            );
            bezierMesh->setSurfaceColor(glm::vec3{ 0.0f, 1.0f, 0.0f });
            bezierMesh->setSmoothShade(false);
            bezierMesh->setEnabled(true);

            bezierMesh->addVertexVectorQuantity(
                "Normals",
                bezierData.normals
            )->setEnabled(showNormals);

            bezierMesh->addVertexVectorQuantity(
                "Tangents U",
                bezierData.tangentsU
            )->setEnabled(showTangents);

            bezierMesh->addVertexVectorQuantity(
                "Tangents V",
                bezierData.tangentsV
            )->setEnabled(showTangents);

            lastBezierWidth = width;
            lastBezierHeight = height;
            lastBezierRadius = radius;
            lastBezierSubdivision = subdivision;

            lastShowNormals = showNormals;
            lastShowTangents = showTangents;
        }

        bezierMesh->setEnabled(true);
    }
    else {
        if (bezierMesh != nullptr) {
            bezierMesh->setEnabled(false);
        }
    }

    ImGui::Text("MLS Surface");
    ImGui::Separator();

    ImGui::Checkbox("Show Surface Mesh##MLS", &showMLSSurfaceMesh);
    ImGui::Checkbox("Show Normals##MLS", &showMLSNormals);

    ImGui::SetNextItemWidth(225);
    ImGui::SliderInt("##MLSSubdivision", &MLSSubdivision, 1, 20);
    ImGui::SameLine();
    ImGui::Text("Subdivision");

    ImGui::SetNextItemWidth(225);
    ImGui::SliderFloat("##MLSRadiusSlider", &MLSradius, 0.001, 1);
    ImGui::SameLine();
    ImGui::Text("Radius");

    if (showMLSSurfaceMesh && sds != nullptr) {

        if (mlsMesh == nullptr ||
            lastMLSWidth != width ||
            lastMLSHeight != height ||
            lastMLSRadius != MLSradius ||
            lastMLSSubdivision != MLSSubdivision ||
            lastShowMLSNormals != showMLSNormals) {

            if (mlsMesh != nullptr) {
                polyscope::removeStructure(mlsMesh);
                mlsMesh = nullptr;
            }

            std::vector<Point> points = sds->getPoints();

            int highResWidth = MLSSubdivision * width;
            int highResHeight = MLSSubdivision * height;

            std::vector<Point> highResGridVertices;
            std::vector<std::array<size_t, 2>> highResGridEdges;

            polyscope::CurveNetwork* tempGrid = renderGrid(
                "Temporary MLS Grid",
                highResWidth,
                highResHeight,
                points,
                highResGridVertices,
                highResGridEdges
            );

            if (tempGrid != nullptr) {
                polyscope::removeStructure(tempGrid);
                tempGrid = nullptr;
            }

            MLSSurfaceData mlsData =
                generateMLSSurfaceData(
                    *sds,
                    highResGridVertices,
                    MLSradius
                );

            std::vector<std::array<size_t, 4>> faces;
            faces.reserve(highResWidth* highResHeight);

            auto idx = [&](int x, int y) -> size_t {
                return y * (highResWidth + 1) + x;
                };

            for (int y = 0; y < highResHeight; ++y) {
                for (int x = 0; x < highResWidth; ++x) {
                    faces.push_back({
                        idx(x, y),
                        idx(x + 1, y),
                        idx(x + 1, y + 1),
                        idx(x, y + 1)
                        });
                }
            }

            mlsMesh = polyscope::registerSurfaceMesh(
                "MLS Surface",
                mlsData.vertices,
                faces
            );

            mlsMesh->addVertexVectorQuantity(
                "MLS Normals",
                mlsData.normals
            )->setEnabled(showMLSNormals);

            std::cout << "highResGridVertices: " << highResGridVertices.size() << std::endl;
            std::cout << "mlsVertices: " << mlsData.vertices.size() << std::endl;
            std::cout << "mlsNormals: " << mlsData.normals.size() << std::endl;
            std::cout << "faces: " << faces.size() << std::endl;

            mlsMesh->setSurfaceColor(glm::vec3{ 0.2f, 0.6f, 1.0f });
            mlsMesh->setSmoothShade(false);
            mlsMesh->setEnabled(true);

            lastMLSWidth = width;
            lastMLSHeight = height;
            lastMLSRadius = MLSradius;
            lastMLSSubdivision = MLSSubdivision;
            lastShowMLSNormals = showMLSNormals;
        }

        mlsMesh->setEnabled(true);
    }
    else {
        if (mlsMesh != nullptr) {
            mlsMesh->setEnabled(false);
        }
    }


}



int main(int argc, char **argv) {
    args::ArgumentParser parser("Computer Graphics 2 Sample Code.");

    args::Flag benchmarkFlag(parser, "benchmark",
                             "Run benchmark after loading an OFF file",
                             {'b', "benchmark"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help &) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    gRunBenchmark = benchmarkFlag;

    polyscope::options::groundPlaneMode =
        polyscope::GroundPlaneMode::ShadowOnly;

    polyscope::options::shadowBlurIters = 6;

    polyscope::init();

    polyscope::view::lookAt(glm::vec3(0.f, 0.f, 2.f), // camera position
                            glm::vec3(0.f, 0.f, 0.f)   // look at
    );

    polyscope::state::userCallback = callback;

    polyscope::show();

    return 0;
}