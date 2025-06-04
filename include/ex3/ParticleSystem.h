#ifndef __ParticleSystem_h
#define __ParticleSystem_h

#include "physx/Actor.h"
#include "ParticleBuffer.h"

namespace tcii::physx
{ // begin namespace tcii::physx

using namespace cg;

class ParticleSystemBase: public Actor
{
protected:
  ParticleSystemBase()
  {
    setName("PS %d", ++_nextId);
  }

private:
  inline static int _nextId;

}; // ParticleSystemBase
  
template <typename... Fields>
class ParticleSystem: public ParticleSystemBase
{
public:
  static ParticleSystem* New()
  {
    return new ParticleSystem;
  }

  void setParticleBuffer(unsigned capacity);

  auto* particles() const
  {
    return (PA*)_particles;
  }

  Bounds3f bounds() const override;
  void render(Renderer&) const override;

private:
  using PA = ParticleBuffer<Fields...>;

  ObjectPtr<PA> _particles;

}; // ParticleSystem

template <typename... Fields>
void
ParticleSystem<Fields...>::setParticleBuffer(unsigned capacity)
{
  if (_particles)
  {
    if (_particles->capacity() == capacity)
      return;
    _particles = nullptr;
  }
  else if (!capacity)
    return;
  _particles = new PA{capacity};
}

template <typename... Fields>
Bounds3f
ParticleSystem<Fields...>::bounds() const
{
  return _particles ? _particles->bounds() : Bounds3f{};
}

template <typename... Fields>
void
ParticleSystem<Fields...>::render(Renderer& renderer) const
{
  // TODO
}

} // end namespace tcii::physx

#endif // __ParticleSystem_h
