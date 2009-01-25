#ifndef CS_UTIL_H
#define CS_UTIL_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Collection of commonly used inline utility functions.

#include <cstdlib>
#include <cctype>
#include <cmath>

#include <limits>
#include <string>
#include <vector>

namespace cs
{

// Returns the base 2 logarithm of num.
inline float log2(float x)
{
    return 1.442695041 * log(x);
}

// Returns the base 10 logarithm of num.
inline float log10(float x)
{
    return 0.434294481 * log(x);
}

// Round to the nearest integer.
inline int iround(double x)
{
    return static_cast<int>(floor(x + 0.5));
}

// // Returns pointer to first non-white-space character in str OR to NULL if none found
// inline const char* strscn(const char* str)
// {
//     if (!str) return NULL;
//     const char* ptr=str;
//     while (*ptr!='\0' && isspace(*ptr)) ptr++;
//     return (*ptr=='\0') ? NULL : ptr;
// }

// // Returns leftmost integer in ptr and sets the pointer to first char after
// // the integer. If no integer is found, returns INT_MIN and sets ptr to NULL
// inline int strtoi(const char*& ptr)
// {
//     const char* ptr0=ptr;
//     if (!ptr) return std::numeric_limits<int>::min();
//     while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
//     if (*ptr=='\0') {
//         ptr=NULL;
//         return std::numeric_limits<int>::min();
//     }
//     int i = (ptr>ptr0 && *(ptr-1)=='-') ? atoi(ptr-1) : atoi(ptr);
//     while (*ptr>='0' && *ptr<='9') ptr++;
//     return i;
// }

// // Same as strtoi, but interpretes '*' as default.
// inline int strtoi_asterix(const char*& ptr, int deflt = std::numeric_limits<int>::max())
// {
//     const char* ptr0=ptr;
//     if (!ptr) return std::numeric_limits<int>::min();
//     while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9') && *ptr!='*') ptr++;
//     if (*ptr=='\0') {
//         ptr=NULL;
//         return std::numeric_limits<int>::min();
//     }
//     if (*ptr=='*') {
//         ptr++;
//         return deflt;
//     }
//     int i = (ptr>ptr0 && *(ptr-1)=='-') ? atoi(ptr-1) : atoi(ptr);
//     while (*ptr>='0' && *ptr<='9') ptr++;
//     return i;
// }

// Normalize a float array such that it sums to one
// If it sums to 0 then assign def_array elements to array (optional)
inline float normalize_to_one(float* array, int length, const float* default_array = NULL)
{
    float sum = 0.0f;
    for (int i = 0; i < length; ++i) sum += array[i];
    if (sum != 0.0f) {
        float fac = 1.0f / sum;
        for (int i = 0; i < length; ++i) array[i] *= fac;
    } else if (default_array) {
        for (int i = 0; i < length; ++i) array[i] = default_array[i];
    }
    return sum;
}

// Reset all entries in a float array to given value or zero if none provided.
inline void reset(float* array, int length, float value = 0.0f)
{
    for (int i = 0; i < length; ++i) array[i] = value;
}

// Splits a string into multiple strings with character delimiter removed.
inline void split(const std::string& s, char c, std::vector<std::string>& v)
{
    std::string::size_type i = 0;
    std::string::size_type j = s.find(c);

    while (j != std::string::npos) {
        v.push_back(s.substr(i, j-i));
        i = ++j;
        j = s.find(c, j);

        if (j == std::string::npos)
            v.push_back(s.substr(i, s.length()));
    }
}

}  // cs

#endif
