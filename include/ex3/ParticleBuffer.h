#ifndef __ParticleBuffer_h
#define __ParticleBuffer_h

#include "graphics/Bounds3.h"
#include "util/SharedObject.h"
#include "util/SoA.h"
#include "physx/Renderer.h"
#include "p1/DefaultSoAAllocator.h"

namespace tcii::physx
{ // begin namespace tcii::physx

using namespace cg;

template <typename... Fields> class ParticleSystem;

template <typename... Fields>
class ParticleBuffer: public SharedObject
{
public:

  using ParticleSoA = cg::SoA<DefaultSoAAllocator, unsigned int, Vec3f, Fields...>;

  auto particleCount() const
  {
    return _particleCount;
  }

  auto capacity() const
  {
    return _capacity;
  }

  void add(const Vec3f& p, const Fields&... fields)
  {
    
    if (_particleCount < _capacity)
    {
      _data.set(_particleCount, p, fields...);
      _particleCount++;
    }

  }

  void clear()
  {
    _particleCount = 0;
  }

  auto& position(unsigned i) const
  {
    return _data.template get<0>(i);
  }

  Bounds3f bounds() const;
  void render(Renderer&) const;

private:

  unsigned int _capacity{};
  unsigned _particleCount{};

  ParticleSoA _data;

  ParticleBuffer(unsigned capacity): 
  _capacity(capacity),
  _particleCount(0),
  _data(capacity)
  {
    // TODO
  }

  friend ParticleSystem<Fields...>;

}; // ParticleBuffer

template <typename... Fields>
Bounds3f
ParticleBuffer<Fields...>::bounds() const
{
  Bounds3f b;
  
  if (_particleCount == 0)
    return b; 

  for (unsigned i = 0; i < _particleCount; ++i)
    b.inflate(position(i)); 
  
  return b;
}

template <typename... Fields>
void
ParticleBuffer<Fields...>::render(Renderer& renderer) const
{
  // TODO
}

} // end namespace tcii::physx

#endif // __ParticleBuffer_h
