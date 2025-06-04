#ifndef __KdTree_h
#define __KdTree_h

#include "KNN.h"
#include "Utils.h"
#include <functional>
#include <utility>
#include <algorithm>
#include <numeric>

namespace tcii::p1
{ // begin namespace tcii::p1

using namespace cg;

template <size_t D, IsReal R, IsPointArray<Vec, R, D> A>
class PointHolder
{
public:
  using Point = Vec<D, R>;

  PointHolder(A&& points):
    _points{std::move(points)}
  {
    // do nothing
  }

  auto& points() const
  {
    return _points;
  }

  auto& points() 
  { 
    return _points; 
  }

private:
  A _points;

}; // PointHolder

template <size_t D, typename R, typename A>
class KdTree: public PointHolder<D, R, A>
{
public:
  using Base = PointHolder<D, R, A>;
  using Bounds = cg::Bounds<D, R>;
  using Point = typename Base::Point;
  using PointFunc = std::function<bool(const A&, unsigned)>;
  using KNN = p1::KNN<R>;

  const struct Params
  {
    unsigned maxPointsPerNode;
    unsigned minPointsPerNode;
    unsigned maxDepth;

  } params;

  static constexpr Params dflParams()
  {
    return {20, 5, 8};
  }

  KdTree(A&& points, const Params& params = dflParams());

  KdTree(const KdTree&) = delete;
  KdTree& operator =(const KdTree&) = delete;

  KdTree(KdTree&& other) noexcept:
    Base{std::move(other)},
    _root{std::exchange(other._root, nullptr)},
    _nodeCount{other._nodeCount},
    _leafCount{other._leafCount},
    params{other.params}
  {
    // do nothing
  }

  ~KdTree()
  {
    delete _root;
  }

  auto& bounds() const
  {
    return _root->bounds;
  }

  auto nodeCount() const
  {
    return _nodeCount;
  }

  auto leafCount() const
  {
    return _leafCount;
  }

  KNN findNeighbors(const Point& point,
    unsigned k,
    PointFunc filter = {}) const;

  auto findNeighbors(unsigned index, unsigned k, PointFunc filter = {}) const
  {
    assert(index < this->points().size());
    return findNeighbors(this->points()[index], k, filter);
  }

  void forEachNeighbor(const Point& point,
    R radius,
    PointFunc f,
    PointFunc filter = {}) const;

  void forEachNeighbor(unsigned index,
    R radius,
    PointFunc f,
    PointFunc filter = {}) const
  {
    assert(index < this->points().size());
    forEachNeighbor(this->points()[index], radius, f, filter);
  }

private:
  struct Node
  {
    Bounds bounds;
    unsigned depth;
    Node* children[2]; 
    int splitAxis;
    R splitValue;
    unsigned firstPointIndex;
    unsigned numPoints;

    Node(const Bounds& bounds, unsigned depth): 
      bounds(bounds), depth(depth), splitAxis(-1), firstPointIndex(0), numPoints(0)
    {
      children[0] = nullptr;
      children[1] = nullptr;
    }

    ~Node()
    {
      delete children[0]; 
      delete children[1]; 
    }

    bool isLeaf() const
    {
      return children[0] == nullptr;
    }

  }; // Node

  std::vector<unsigned int> _point_indices; 

  Node* buildRecursive(Bounds currentBounds, unsigned depth,
  unsigned firstIdx, unsigned numPtsInNode,
  int initial_axis);

  Node* _root{};
  unsigned _nodeCount{};
  unsigned _leafCount{};

}; // KdTree

template <size_t D, typename R, typename A>
KdTree<D, R, A>::KdTree(A&& points, const Params& params):
  Base{std::move(points)},
  params{params},
  _root{nullptr},
  _nodeCount{0},
  _leafCount{0}
{

  if (this->points().size() == 0)
    return;

  Bounds rootBounds = computeBounds<D, R>(this->points());

  int initialAxis = 0; 

  if constexpr (D > 0) { 

    R maxSpan = R{-1}; 

    for (size_t i = 0; i < D; ++i) {

      R span = rootBounds.max()[i] - rootBounds.min()[i];

      if (span > maxSpan) {

        maxSpan = span;
        initialAxis = static_cast<int>(i);

      }

    }

  }

  if (this->points().size() != 0) { 
    _point_indices.resize(this->points().size());
    std::iota(_point_indices.begin(), _point_indices.end(), 0);
  }

  if (!_point_indices.empty()) 
    _root = buildRecursive(rootBounds, 0, 0, _point_indices.size(), initialAxis);

}

template <size_t D, typename R, typename A>
typename KdTree<D, R, A>::Node*
KdTree<D, R, A>::buildRecursive(Bounds currentBounds, unsigned depth,
  unsigned firstIdx, unsigned numPtsInNode,
  int initial_axis)
{
  
  Node* node = new Node(currentBounds, depth);
  _nodeCount++;

  if (numPtsInNode == 0) { 
    node->firstPointIndex = firstIdx; 
    node->numPoints = 0;
    _leafCount++;
    return node; 
  }

  if (numPtsInNode <= params.maxPointsPerNode || depth >= params.maxDepth) {
    node->firstPointIndex = firstIdx;
    node->numPoints = numPtsInNode;
    _leafCount++;
    return node;
  }

  node->splitAxis = (initial_axis + depth) % D;

  unsigned m_count = (numPtsInNode + 1) / 2;
  unsigned median_offset = m_count - 1;      
  unsigned medianAbsoluteIndex = firstIdx + median_offset;

  std::nth_element(
    _point_indices.begin() + firstIdx,
    _point_indices.begin() + medianAbsoluteIndex, 
    _point_indices.begin() + firstIdx + numPtsInNode,
    [&](unsigned int index1, unsigned int index2) { 
      return this->points()[index1][node->splitAxis] < this->points()[index2][node->splitAxis];
    }
  );

  unsigned int actualMedianPointOriginalIndex = _point_indices[medianAbsoluteIndex];
  node->splitValue = this->points()[actualMedianPointOriginalIndex][node->splitAxis];

  Bounds leftChildBounds = currentBounds;
  leftChildBounds.max()[node->splitAxis] = node->splitValue;

  Bounds rightChildBounds = currentBounds;
  rightChildBounds.min()[node->splitAxis] = node->splitValue;

  node->children[0] = buildRecursive(leftChildBounds, 
                                      depth + 1,
                                      firstIdx, m_count,
                                      initial_axis);

  node->children[1] = buildRecursive(rightChildBounds, 
                                      depth + 1,
                                      firstIdx + m_count, 
                                      numPtsInNode - m_count,
                                      initial_axis);

  return node;

}

template <size_t D, typename R, typename A>
auto
KdTree<D, R, A>::findNeighbors(const Point& point,
  unsigned k,
  PointFunc filter) const -> KNN
{
  KNN knn{k};

  // TODO
  return knn;
}

template <size_t D, typename R, typename A>
void
KdTree<D, R, A>::forEachNeighbor(const Point& point,
  R radius,
  PointFunc f,
  PointFunc filter) const
{
  f(this->points(), 0);
  // TODO
}

template <typename R, typename A> using KdTree3 = KdTree<3, R, A>;

} // end namespace tcii::p1

#endif // __KdTree_h
