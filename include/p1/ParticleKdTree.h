#ifndef __ParticleKdTree_h
#define __ParticleKdTree_h

#include "KdTree.h"
#include "ex3/ParticleSystem.h"

namespace tcii::p1
{ // begin namespace tcii::p1

using namespace cg;
using namespace physx;

template <typename... Fields>
class ParticleArray
{
public:
  ParticleArray(const ParticleSystem<Fields...>& ps):
    _particles{ps.particles()}
  {
    assert(_particles != nullptr);
  }

  ParticleArray(ParticleArray&&) noexcept = default;

  size_t size() const
  {
    return _particles->particleCount();
  }

  auto& operator [](size_t i) const
  {
    return _particles->position(unsigned(i));
  }

private:
  using PA = ParticleBuffer<Fields...>;

  ObjectPtr<PA> _particles;

}; // ParticleArray

template <typename... Fields>
using ParticleKdTree = KdTree3<float, ParticleArray<Fields...>>;

} // end namespace tcii::p1

#endif // __ParticleKdTree_h
