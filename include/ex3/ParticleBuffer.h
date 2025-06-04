#ifndef __ParticleBuffer_h
#define __ParticleBuffer_h

#include "graphics/Bounds3.h"
#include "util/SharedObject.h"
#include "util/SoA.h"
#include "physx/Renderer.h"

namespace tcii::physx
{ // begin namespace tcii::physx

using namespace cg;

template <typename... Fields> class ParticleSystem;

template <typename... Fields>
class ParticleBuffer: public SharedObject
{
public:
  auto particleCount() const
  {
    return _particleCount;
  }

  auto capacity() const
  {
    // TODO
    return 0u;
  }

  void add(const Vec3f& p, const Fields&... fields)
  {
    // TODO
    ++_particleCount;
  }

  void clear()
  {
    // TODO
  }

  auto& position(unsigned i) const
  {
    // TODO
    static Vec3f _dummy;
    return _dummy;
  }

  Bounds3f bounds() const;
  void render(Renderer&) const;

private:
  unsigned _particleCount{};

  ParticleBuffer(unsigned capacity)
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

  // TODO
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
