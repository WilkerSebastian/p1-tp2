#ifndef __Utils_h
#define __Utils_h

#include "graphics/Bounds3.h"
#include "Concepts.h"
#include <vector>

namespace tcii::p1
{ // begin namespace tcii::p1

using namespace cg;

template <size_t D, IsReal R, IsPointArray<Vec, R, D> A>
auto
computeBounds(const A& points)
{
  Bounds<D, R> b;

  for (size_t n = points.size(), i = 0; i < n; ++i)
    b.inflate(points[i]);
  return b;
}

template <size_t D, IsReal R>
auto
distance(const Vec<D, R>& p, const Bounds<D, R>& bounds)
{
  R d{};

  for (size_t i = 0; i < D; ++i)
    if (auto t = p[i] - bounds[0][i]; t < 0)
      d += t * t;
    else if (auto t = p[i] - bounds[1][i]; t > 0)
      d += t * t;
  return (R)std::sqrt(d);
}

template <IsReal R>
inline auto
srand()
{
  return (R)((double)rand() / RAND_MAX);
}

template <IsReal R>
inline auto
srand(R min, R max)
{
  return min + srand<R>() * (max - min);
}

template <size_t D, IsReal R>
inline auto
prand(R min, R max)
{
  Vec<D, R> p;

  for (size_t i = 0; i < D; ++i)
    p[i] = srand(min, max);
  return p;
}

template <size_t D, IsReal R>
using PointVector = std::vector<Vec<D, R>>;

template <size_t D, IsReal R>
inline auto
prand(size_t n, R min = 0, R max = 1)
{
  assert(n > 0);

  PointVector<D, R> v(n);

  for (size_t i = 0; i < n; ++i)
    v[i] = prand<D>(min, max);
  return v;
}

template <typename R> using Points3 = PointVector<3, R>;

} // end namespace tcii::p1

#endif // __Utils_h
