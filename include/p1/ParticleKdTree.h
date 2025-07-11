/*
* Autor: Wilker Sebastian Afonso Pereira
* GitHub: https://github.com/WilkerSebastian/p1-tp2
*
* Prova 1 de Tópicos em Computação 2
*/
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

  const auto& getColor(size_t i) const
  {
    return _particles->getColor(unsigned(i));
  }

private:
  using PA = ParticleBuffer<Fields...>;

  ObjectPtr<PA> _particles;

}; // ParticleArray

template <typename... Fields>
using ParticleKdTree = KdTree3<float, ParticleArray<Fields...>>;

} // end namespace tcii::p1

#endif // __ParticleKdTree_h
