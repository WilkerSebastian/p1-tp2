#include "p1/ParticleKdTree.h"
#include "p1/Utils.h"
#include <iostream>

using namespace tcii::p1;

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

int
main()
{
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
  puts("Press any key to exit...");
  (void)getchar();
  return 0;
}
