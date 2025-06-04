#ifndef __Actor_h
#define __Actor_h

#include "graphics/Bounds3.h"
#include "util/NameableObject.h"
#include "ObjectList.h"

namespace tcii::physx
{ // begin namespace tcii::physx

using namespace cg;

class Renderer;
class Scene;

class Actor: ObjectListNode<Actor>, public NameableObject
{
public:
  virtual ~Actor() = default;

  auto scene() const
  {
    return _scene;
  }

  template <typename T> bool is() const
  {
    return dynamic_cast<const T*>(this);
  }

  virtual Bounds3f bounds() const = 0;
  virtual void render(Renderer&) const = 0;

private:
  Scene* _scene{};

  friend Scene;

}; // Actor

} // end namespace tcii::physx

#endif // __Actor_h
