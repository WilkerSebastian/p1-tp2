#include "p1/ParticleKdTree.h"
#include "p1/Utils.h"
#include "p1/Vec2.h"
#include "p1/Bounds2.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace tcii::p1;
using namespace tcii::cg;
using namespace tcii::physx;

template <typename PointType, typename R_Type, typename PointArrayT>
void printKnnSearchResults(
    const PointType& queryPoint,
    const std::vector<std::pair<R_Type, unsigned int>>& results,
    const PointArrayT& point_array)
{

  std::cout << "Query Point: " << queryPoint << "\n";

  if (results.empty()) {
    std::cout << "No neighbors found!\n";
    return;
  }
  
  std::cout << "Found " << results.size() << " Neighbors (sorted by distance):\n";

  for (const auto& neighbor_pair : results) 
    std::cout << "Dist: " << neighbor_pair.first
              << ", Index: " << neighbor_pair.second
              << ", Coords: " << point_array[neighbor_pair.second] << "\n";
    
}

template <size_t D, typename R, typename TreeType>
void testKdTreeFeatures(
  TreeType& tree,
  const std::string& testName,
  typename TreeType::PointFunc filter = nullptr)
{
  std::cout << "\n\n" << std::string(60, '=') << "\n";
  std::cout << "Starting Test: " << testName << "\n";
  std::cout << std::string(60, '-') << "\n";

  std::cout << std::fixed << std::setprecision(4);

  const auto& points = tree.points();

  if (points.size() == 0) {
    std::cout << "Tree is empty. Skipping search tests!\n";
    std::cout << std::string(60, '=') << "\n";
    return;
  }

  std::cout << "Tree Bounds: " << tree.bounds() << "\n";
  std::cout << "Node Count: " << tree.nodeCount() << ", Leaf Count: " << tree.leafCount() << "\n";

  unsigned queryIndex = points.size() / 3;
  typename TreeType::Point queryPoint = points[queryIndex];

  std::cout << "\nTesting findNeighbors for point at index " << queryIndex << "\n";
  
  unsigned k_values[] = {1, 5, 10};

  for (unsigned k : k_values) {

    if (k > points.size()) 
      k = points.size();

    if (k == 0) 
      continue;
    
    std::cout << "\nSearching for k = " << k << " neighbors...\n";

    if (filter)  
      std::cout << "(Using filter)\n";
    
    typename TreeType::KNN knn_results = tree.findNeighbors(queryPoint, k, filter);
    printKnnSearchResults(queryPoint, knn_results.getSortedResults(), points);

  }

  std::cout << "\nTesting forEachNeighbor for point at index " << queryIndex << "\n";

  R radius_values[] = {static_cast<R>(0.1), static_cast<R>(0.25)};

  for (R radius : radius_values) {

    std::cout << "\nSearching with radius = " << radius << "...\n";

    if (filter)
      std::cout << "(Using filter)\n";
    
    unsigned neighbors_found_count = 0;

    tree.forEachNeighbor(queryPoint, radius,
        [&](const auto& point_array_ref, unsigned int found_idx) {
            if (neighbors_found_count == 0) {
                  std::cout << "Found Neighbors:\n";
            }
            std::cout << "Index: " << found_idx << ", Coords: " << point_array_ref[found_idx] << "\n";
            neighbors_found_count++;
            return true; // Continue searching
        },
        filter
    );

    if (neighbors_found_count == 0) 
      std::cout << "No neighbors found in radius!\n";
    
  }

  std::cout << std::string(60, '=') << "\n";

}

template <typename Tree>
inline void
testTree(Tree&& tree)
{
  auto& points = tree.points();
  auto& bounds = tree.bounds();

  std::cout.precision(3);
  std::cout << "Tree bounds: " << bounds << '\n';
  std::cout << "Points: " << points.size() << '\n';

  typename Tree::Point::value_type d{};

  for (unsigned i = 0; i < points.size(); ++i)
    d += distance(points[i], bounds);
  std::cout << "Distance sum: " << d << '\n';
  (void)tree.findNeighbors(0, 10);
  tree.forEachNeighbor(0,
    0.5f,
    [](const auto& points, unsigned index)
    {
      std::cout << points[index] << '\n';
      return true;
    });
}

inline void testKdTreePointVectorFloat() {

  constexpr size_t D = 3;
  using R = float;
  using PointArrayT = PointVector<D, R>;
  
  unsigned num_points_small = 500;
  auto points_small = prand<D, R>(num_points_small, R{0.0}, R{1.0});
  KdTree<D, R, PointArrayT> tree_small(std::move(points_small));
  testKdTreeFeatures<D,R>(tree_small, "3D KdTree<float> with PointVector (500 points)");

  unsigned num_points_large = 200000; 
  auto points_large = prand<D, R>(num_points_large, R{0.0}, R{1.0});
  KdTree<D, R, PointArrayT> tree_large(std::move(points_large));
  testKdTreeFeatures<D,R>(tree_large, "3D KdTree<float> with PointVector (200k points)");

}

inline void testKdTreePointVectorDouble() {

  constexpr size_t D = 3;
  using R = double; 
  using PointArrayT = PointVector<D, R>;
  
  unsigned num_points = 500;
  auto points = prand<D, R>(num_points, R{0.0}, R{1.0});
  KdTree<D, R, PointArrayT> tree_double(std::move(points));
  testKdTreeFeatures<D,R>(tree_double, "3D KdTree<double> with PointVector (500 points)");

}

inline void test2DKdTreePointVectorFloat() {

  constexpr size_t D = 2;
  using R = float;
  using PointArrayT = PointVector<D, R>;

  unsigned num_points = 500;
  auto points_2d = prand<D, R>(num_points, R{0.0}, R{1.0});
  KdTree<D, R, PointArrayT> tree_2d(std::move(points_2d));
  testKdTreeFeatures<D,R>(tree_2d, "2D KdTree<float> with PointVector (500 points)");

}

inline void testParticleKdTreeWithColorFilter() {

  auto ps = ParticleSystem<Vec3f>::New();
  ps->setParticleBuffer(500);
  Vec3f red{1,0,0}, green{0,1,0}, blue{0,0,1};
  Vec3f colors[] = {red, green, blue};

  for (unsigned i = 0; i < 500; ++i) 
    ps->particles()->add(prand<3, float>(0.f, 1.f), colors[i % 3]); 

  auto particle_tree = ParticleKdTree<Vec3f>{*ps};

  Vec3f targetColor = red;
  typename decltype(particle_tree)::PointFunc colorFilter = 
      [&](const auto& point_array, unsigned int index) -> bool {
          return point_array.getColor(index) == targetColor;
      };

  testKdTreeFeatures<3, float>(particle_tree, "ParticleKdTree WITH Color Filter (searching for RED)", colorFilter);
  testKdTreeFeatures<3, float>(particle_tree, "ParticleKdTree WITHOUT Filter", nullptr); 

  destroy(ps);

}


int
main()
{

  srand(time(0));

  constexpr size_t D{3};
  using R = float;
  using A = PointVector<D, R>;
  constexpr unsigned np{100};

  testTree(KdTree<D, R, A>{prand<D, R>(np)});

  auto ps = ParticleSystem<Vec3f, Vec3f>::New();

  ps->setParticleBuffer(np);
  for (unsigned i = 0; i < np; ++i)
    ps->particles()->add(prand<3, float>(0.f, 1.f), {}, {});

  auto tree = ParticleKdTree<Vec3f, Vec3f>{*ps};

  destroy(ps);
  testTree(std::move(tree));

  testKdTreePointVectorFloat();

  testKdTreePointVectorDouble();

  test2DKdTreePointVectorFloat();

  testParticleKdTreeWithColorFilter();

  puts("Press any key to exit...");
  (void)getchar();
  return 0;
}
