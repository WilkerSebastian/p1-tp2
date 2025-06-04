#ifndef __Concepts_h
#define __Concepts_h

#include <vector>
#include <concepts>

namespace tcii::p1
{ // begin namespace tcii::p1

template <typename T, typename U>
concept SameAs = std::same_as<T, U>;

template <typename R>
concept IsReal = std::floating_point<R>;

template <typename C, typename E>
concept IsIndexable = requires (C x, size_t i)
{
  { x[i] } -> std::convertible_to<E>;
};

template <template <size_t, typename> typename P, typename R, size_t D>
concept IsPoint = IsReal<R> &&
  SameAs<R, typename P<D, R>::value_type> &&
  IsIndexable<P<D, R>, R> &&
requires (P<D, R> p, R r)
{
  { p + p } -> SameAs<P<D, R>>;
  { p - p } -> SameAs<P<D, R>>;
  { p * r } -> SameAs<P<D, R>>;
  { r * p } -> SameAs<P<D, R>>;
};

template <typename A,
  template <size_t, typename> typename P,
  typename R,
  size_t D>
concept IsPointArray = IsPoint<P, R, D> &&
  IsIndexable<A, P<D, R>> &&
  std::move_constructible<A> &&
requires (A a, size_t i)
{
  { a.size() } -> SameAs<size_t>;
};

} // end namespace tcii::p1

#endif // __Concepts_h
