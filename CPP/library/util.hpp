#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <istream>
#include <ostream>
#include <cmath>
#include <vector>
#include <memory>
#include <exception>
#include <stdexcept>

template <typename T>
inline void SWAP(T *a, T *b)
{
    T temp = *a;
    *a = *b;
    *b = temp;
}

template <typename T>
inline T SQR(T a)
{
    return a * a;
}

#endif