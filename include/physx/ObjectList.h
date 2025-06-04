#ifndef __ObjectList_h
#define __ObjectList_h

#include <cassert>

namespace tcii::physx
{ // begin namespace tcii::physx

template <typename T>
inline void
destroy(T* object)
{
  delete object;
}

template <typename T> class ObjectList;
template <typename T> class ObjectListIterator;

template <typename T>
class ObjectListNode
{
  ObjectList<T>* _list;
  T* _next;
  T* _prev;

  friend ObjectList<T>;

}; // ObjectListNode

template <typename T>
class ObjectListIterator
{
public:
  ObjectListIterator(const T* node, const ObjectList<T>* list):
#ifdef DEBUG
    _list{list},
#endif // DEBUG
    _node{node}
  {
    (void)list;
  }

  auto& operator *() const
  {
#ifdef DEBUG
    if (!_node)
      throw std::logic_error{"List iterator not dereferencable"};
#endif // DEBUG
    return *_node;
  }

  auto& operator ++()
  {
#ifdef DEBUG
    if (!_node)
      throw std::logic_error{"List iterator not incrementable"};
#endif // DEBUG
    _node = ObjectList<T>::next(_node);
    return *this;
  }

  auto operator ++(int)
  {
    auto temp = *this;

    ++*this;
    return temp;
  }

  auto& operator --()
  {
#ifdef DEBUG
    if (!_node)
      throw std::logic_error{"List iterator not decrementable"};
#endif // DEBUG
    _node = ObjectList<T>::prev(_node);
    return *this;
  }

  auto operator --(int)
  {
    auto temp = *this;

    --*this;
    return temp;
  }

  auto operator ->()
  {
    return &(operator *());
  }

  auto operator ==(const ObjectListIterator& other) const
  {
#ifdef _DEBUG
    if (_list != other._list)
      throw std::logic_error{"List iterators incompatible"};
#endif // _DEBUG
    return _node == other._node;
  }

  auto operator !=(const ObjectListIterator& other) const
  {
    return !operator ==(other);
  }

private:
  const T* _node;
#if _DEBUG
  const ObjectList<T>* _list;
#endif // DEBUG

}; // ObjectListIterator

template <typename T>
class ObjectList
{
public:
  ~ObjectList()
  {
    clear();
  }

  ObjectList() = default;

  ObjectList(const ObjectList&) = delete;
  ObjectList& operator =(const ObjectList&) = delete;

  auto empty() const
  {
    return _head == nullptr;
  }

  auto size() const
  {
    return _count;
  }

  void clear();
  bool add(T* node);
  bool remove(T* node);

  ObjectListIterator<T> begin() const
  {
    return {_head, this};
  }

  auto end() const
  {
    return ObjectListIterator<T>{nullptr, this};
  }

private:
  T* _head{};
  T* _tail{};
  unsigned _count{};

  static auto& next(const T* node)
  {
    return ((ObjectListNode<T>*)node)->_next;
  }

  static auto& prev(const T* node)
  {
    return ((ObjectListNode<T>*)node)->_prev;
  }

  static auto& list(const T* node)
  {
    return ((ObjectListNode<T>*)node)->_list;
  }

  friend ObjectListIterator<T>;

}; // ObjectList

template <typename T>
void
ObjectList<T>::clear()
{
  while (auto temp = _head)
  {
    _head = next(temp);
    destroy(temp);
  }
  _tail = nullptr, _count = 0;
}

template <typename T>
bool
ObjectList<T>::add(T* node)
{
  assert(node != nullptr);
  // If the node already belongs to this list,
  // then we return false.
  if (list(node) == this)
    return false;
  // If the node belongs to another list, then
  // we remove it from its owner.
  if (list(node))
    list(node)->remove(node);
  // Next, we insert the node in this list and
  // set this as its owner.
  if (!_head)
    _head = node;
  else
    next(_tail) = node;
  _tail = node;
  prev(node) = _tail;
  next(node) = nullptr;
  ++_count;
  list(node) = this;
  return true;
}

template <typename T>
bool
ObjectList<T>::remove(T* node)
{
  assert(node != nullptr);
  if (list(node) != this)
    return false;
  if (prev(node))
    next(prev(node)) = next(node);
  if (next(node))
    prev(next(node)) = prev(node);
  if (node == _head)
    _head = next(node);
  if (node == _tail)
    _tail = prev(node);
  --_count;
  list(node) = nullptr;
  return true;
}

} // end namespace tcii::physx

#endif // __ObjectList_h
