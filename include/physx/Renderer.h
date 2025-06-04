#ifndef __Renderer_h
#define __Renderer_h

#include "graphics/Vec3.h"

namespace tcii::physx
{ // begin namespace tcii::physx

class Renderer
{
public:
  auto& draw(const char* label, const cg::Vec3f& p)
  {
    std::cout << label << p << '\n';
    return *this;
  }

  auto& draw(const char* label, float s)
  {
    std::cout << label << ':' << s << '\n';
    return *this;
  }

}; // Renderer

} // end namespace tcii::physx

#endif // __Renderer_h
