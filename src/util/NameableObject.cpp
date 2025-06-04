#include "util/NameableObject.h"
#include <cstdarg>
#include <cstdio>

namespace tcii::cg
{ // begin namespace tcii::cg

void
NameableObject::setName(const char* format, ...)
{
  constexpr auto bufferSize = 128;
  char buffer[bufferSize];
  va_list args;

  va_start(args, format);
  vsnprintf(buffer, bufferSize, format, args);
  _name = buffer;
}

} // end namespace tcii::cg
