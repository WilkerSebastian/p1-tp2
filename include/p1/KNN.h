#ifndef __KNN_h
#define __KNN_h

#include "Concepts.h"
#include <vector>
#include <queue>     
#include <algorithm> 
#include <limits>   

namespace tcii::p1
{ // begin namespace tcii::p1

template <IsReal R>
class KNN
{
public:
  struct Neighbor
  {
    R distance;
    unsigned int index;

    bool operator<(const Neighbor& other) const
    {
      return distance < other.distance;
    }
  };

  explicit KNN(unsigned int k) : _k{k}
  {
    if (_k == 0) 
      _k = 1;
  }

  void add(R dist, unsigned int point_idx)
  {
    if (_found_neighbors.size() < _k)
      _found_neighbors.push({dist, point_idx});
    
    else if (dist < _found_neighbors.top().distance) 
    {
      _found_neighbors.pop();
      _found_neighbors.push({dist, point_idx});
    }
  }

  R getWorstDistance() const
  {
    if (_found_neighbors.size() < _k)
      return std::numeric_limits<R>::max(); 
    
    return _found_neighbors.top().distance; 
  }

  bool isFull() const
  {
    return _found_neighbors.size() == _k;
  }

  size_t count() const
  {
    return _found_neighbors.size();
  }

  std::vector<Neighbor> getResults() const
  {
    std::vector<Neighbor> results;
    results.reserve(_found_neighbors.size());
    
    std::priority_queue<Neighbor> temp_pq = _found_neighbors;
    while (!temp_pq.empty())
    {
      results.push_back(temp_pq.top());
      temp_pq.pop();
    }

    std::reverse(results.begin(), results.end());
    return results;
  }

  std::vector<std::pair<R, unsigned int>> getSortedResults() const 
  {
    std::vector<Neighbor> temp_results = getResults(); 
    std::vector<std::pair<R, unsigned int>> final_results;
    final_results.reserve(temp_results.size());

    for(const auto& n : temp_results) 
      final_results.push_back({n.distance, n.index});
    
    return final_results;
  }

private:
  unsigned int _k;
  std::priority_queue<Neighbor> _found_neighbors;

}; // KNN

} // end namespace tcii::p1

#endif // __KNN_h
