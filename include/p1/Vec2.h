#ifndef __Vec2_h
#define __Vec2_h

#include <cassert>
#include <cmath>
#include <iostream>
#include <concepts> 

namespace tcii::cg
{ 
template <size_t D, typename real> struct Vec; 

#ifndef ASSERT_REAL
#define ASSERT_REAL(T, msg) static_assert(std::floating_point<T>, msg)
#endif

template <typename real>
struct Vec<2, real>
{
    ASSERT_REAL(real, "Vec2: floating-point type expected");

    using value_type = real;

    real x;
    real y;

    Vec() : x(real{0}), y(real{0}) {}
    Vec(real x_val, real y_val) : x(x_val), y(y_val) {}

    auto& operator[](size_t i)
    {
        assert(i < 2);
        return (&x)[i]; 
    }

    auto operator[](size_t i) const
    {
        return const_cast<Vec*>(this)->operator[](i); 
    }

    real length() const
    {
        return std::sqrt(x * x + y * y);
    }

    Vec versor() const
    {
        real len = length();

        if (len == real{0}) 
            return Vec(real{0}, real{0}); 

        return (real{1} / len) * (*this);
    }

}; 

template <typename real> using Vec2 = Vec<2, real>;

template <typename real>
inline Vec2<real>
operator+(const Vec2<real>& u, const Vec2<real>& v)
{
  return {u.x + v.x, u.y + v.y};
}

template <typename real>
inline Vec2<real>
operator-(const Vec2<real>& u, const Vec2<real>& v)
{
  return {u.x - v.x, u.y - v.y};
}

template <typename real>
inline Vec2<real>
operator*(real s, const Vec2<real>& v)
{
  return {s * v.x, s * v.y};
}

template <typename real>
inline Vec2<real>
operator*(const Vec2<real>& v, real s)
{
  return s * v;
}

template <typename T>
std::ostream&
operator<<(std::ostream& os, const Vec2<T>& v)
{
  os << '(' << v.x << ',' << v.y << ')';
  return os;
}

template <typename real>
inline bool
operator==(const Vec<2, real>& u, const Vec<2, real>& v) 
{
  return u.x == v.x && u.y == v.y
}

using Vec2f = Vec2<float>;

}

#endif 