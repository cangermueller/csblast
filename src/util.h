#ifndef UTIL_H
#define UTIL_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Collection of commonly used inline utility functions.

#include <cstdlib>
#include <cctype>
#include <climits>
#include <cmath>

#include "my_exception.h"

// Returns the base 2 logarithm of num.
inline double log2(double num)
{
    return log(num)/log(2);
}

// Round to the nearest integer.
inline int iround(double x)
{
    return static_cast<int>(floor(x+0.5));
}

// Returns pointer to first non-white-space character in str OR to NULL if none found
inline const char* strscn(const char* str)
{
    if (!str) return NULL;
    const char* ptr=str;
    while (*ptr!='\0' && isspace(*ptr)) ptr++;
    return (*ptr=='\0') ? NULL : ptr;
}

// Returns leftmost integer in ptr and sets the pointer to first char after
// the integer. If no integer is found, returns INT_MIN and sets ptr to NULL
inline int strtoi(const char*& ptr)
{
    const char* ptr0=ptr;
    if (!ptr) return INT_MIN;
    while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
    if (*ptr=='\0') {
        ptr=NULL;
        return INT_MIN;
    }
    int i = (ptr>ptr0 && *(ptr-1)=='-') ? atoi(ptr-1) : atoi(ptr);
    while (*ptr>='0' && *ptr<='9') ptr++;
    return i;
}

// Same as strtoi, but interpretes '*' as default.
inline int strtoi_asterix(const char*& ptr, int deflt=INT_MAX)
{
    const char* ptr0=ptr;
    if (!ptr) return INT_MIN;
    while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9') && *ptr!='*') ptr++;
    if (*ptr=='\0') {
        ptr=NULL;
        return INT_MIN;
    }
    if (*ptr=='*') {
        ptr++;
        return deflt;
    }
    int i = (ptr>ptr0 && *(ptr-1)=='-') ? atoi(ptr-1) : atoi(ptr);
    while (*ptr>='0' && *ptr<='9') ptr++;
    return i;
}

// Normalizes a float array to one and returns the array sum. In case the array sums to
// zero and default_array is provided the array is set to default_array, otherwise an
// exception is thrown.
inline float normalize_to_one(float* array, int length, const float* default_array=NULL)
{
    float sum = 0.0f;
    for (int i = 0; i < length; ++i) sum += array[i];
    if (sum != 0.0f) {
        float fac = 1.0f / sum;
        for (int  i = 0; i < length; ++i) array[i] *= fac;
    } else if (default_array) {
        for (int i = 0; i < length; ++i) array[i] = default_array[i];
    } else {
        throw MyException("Unable to normalize array to one. Array sum is zero!");
    }
    return sum;
}

// Reset all entries in a float array to given value or zero if none provided.
inline void reset(float* array, int length, float value = 0.0f)
{
    for (int i = 0; i < length; ++i) array[i] = value;
}

#endif
