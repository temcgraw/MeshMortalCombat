#pragma once
#include <glm/glm.hpp>
#include <vector>

template <class T>
struct vector3d
{
   std::vector<T> mVec;
   glm::ivec3 mSize;
   glm::ivec3 mStride;

   vector3d();
   vector3d(const glm::ivec3& size);
   vector3d(const glm::ivec3& size, T val);

   int index(const glm::ivec3& coord) const { return coord.x + coord.y * mStride.y + coord.z * mStride.z; }
   bool valid_index(const glm::ivec3& coord) const
   {
      if (coord.x < 0 || coord.y < 0 || coord.z < 0) return false;
      if (coord.x >= mSize.x || coord.y >= mSize.y || coord.z >= mSize.z) return false;
      return true;
   }

   T get(const glm::ivec3& coord) const;
   T get(int i, int j, int k) const;
   T get_clamp(int i, int j, int k) const;
   T& getRef(const glm::ivec3& coord);
   T& getRef(int i, int j, int k);

   void set(int i, int j, int k, const T& val);
   void set(const glm::ivec3& coord, const T& val);
   void fill(const T& val);

   T* data() { return mVec.data(); }

   void resize(const glm::ivec3& size);
   glm::ivec3 size() const { return mSize; }
   int size1d() const { return mVec.size(); }
};

template <class T>
vector3d<T>::vector3d() :mSize(0), mStride(0)
{

}

template <class T>
vector3d<T>::vector3d(const glm::ivec3& size) :mSize(size), mVec(size.x* size.y* size.z)
{
   mStride = glm::ivec3(1, mSize.x, mSize.x * mSize.y);
}

template <class T>
vector3d<T>::vector3d(const glm::ivec3& size, T val) :mSize(size), mVec(size.x* size.y* size.z, val)
{
   mStride = glm::ivec3(1, mSize.x, mSize.x * mSize.y);
}

template <class T>
void vector3d<T>::resize(const glm::ivec3& size)
{
   mSize = size;
   mVec.resize(size.x * size.y * size.z);
   mStride = glm::ivec3(1, mSize.x, mSize.x * mSize.y);
}

template <class T>
T vector3d<T>::get(int i, int j, int k) const
{
   int ix = index(glm::ivec3(i, j, k));
   return mVec[ix];
}

template <class T>
T vector3d<T>::get(const glm::ivec3& coord) const
{
   int ix = index(coord);
   return mVec[ix];
}

// added reference getter
template <class T>
T& vector3d<T>::getRef(int i, int j, int k)
{
    int ix = index(glm::ivec3(i, j, k));
    return mVec[ix];
}

template <class T>
T& vector3d<T>::getRef(const glm::ivec3& coord)
{
    int ix = index(coord);
    return mVec[ix];
}


template <class T>
T vector3d<T>::get_clamp(int i, int j, int k) const
{
   glm::ivec3 coord = glm::clamp(glm::ivec3(i, j, k), glm::ivec3(0), mSize - glm::ivec3(1));
   return get(coord);
}

template <class T>
void vector3d<T>::set(const glm::ivec3& coord, const T& val)
{
   int ix = index(coord);
   mVec[ix] = val;
}

template <class T>
void vector3d<T>::set(int i, int j, int k, const T& val)
{
   int ix = index(glm::ivec3(i, j, k));
   mVec[ix] = val;
}

template <class T>
void vector3d<T>::fill(const T& val)
{
   for (T& v : mVec)
   {
      v = val;
   }
}