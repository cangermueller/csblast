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
#include <cstdarg>

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace cs
{

// Returns the base 2 logarithm of x.
inline float log2(float x)
{
    return 1.442695041 * log(x);
}

// Round to the nearest integer.
inline int iround(double x)
{
    return static_cast<int>(floor(x + 0.5));
}

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

// Stringifies the range elemeents delimited by given character.
template<typename Fwd>
inline std::string stringify_range(Fwd first, Fwd last, char delim = ',')
{
    std::ostringstream out;
    out << "{";
    while (first != last) {
        out << *first;
        if (++first != last)
            out << delim << ' ';
    }
    out << "}";
    return out.str();
}

// Stringifies the provided container delimited by given character.
template<typename C>
inline std::string stringify_container(const C& c, char delim = ',')
{
    return stringify_range(c.begin(), c.end(), delim);
}

// sprintf-like helper that returns a string.
inline std::string strprintf(const char* str, ...)
{
    char *buffer = new char[1000];
    va_list ap;
    va_start(ap, str);
    vsprintf(buffer, str, ap);
    va_end(ap);
    std::string rv(buffer);
    delete [] buffer;
    return rv;
}

// Returns the file extension of the given file name.
inline std::string get_file_ext(const std::string& s)
{
    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos) {
        return(s.substr(i+1, s.length() - i));
    }

    return("");
}

// Returns the last component of the filename. If ext is false
// the file extension is removed.
inline std::string get_file_basename(const std::string& s, bool ext = true)
{
    char sep = '/';

#ifdef _WIN32
    sep = '\\';
#endif

    size_t i = s.rfind(sep, s.length());
    if (i != std::string::npos) {
        if (ext) {
            return(s.substr(i+1, s.length() - i));
        } else {
            size_t j = s.rfind('.', s.length());
            if (j != std::string::npos && j > i) {
                return(s.substr(i+1, j-i));
            } else {
                return(s.substr(i+1, s.length() - i));
            }
        }
    }

    return("");
}

// Returns all components of the filename except the last one.
inline std::string get_file_dirname(const std::string& s)
{
    char sep = '/';

#ifdef _WIN32
    sep = '\\';
#endif

    size_t i = s.rfind(sep, s.length( ));
    if (i != std::string::npos) {
        return(s.substr(0, i));
    }

    return("");
}


}  // cs

#endif
