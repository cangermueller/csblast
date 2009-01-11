#ifndef UTIL_H
#define UTIL_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Collection of commonly used inline utility functions.

#include <cstdlib>
#include <cctype>      // isspace
#include <climits>     // INT_MIN

// Returns the base 2 logarithm of num
inline double log2(double num) { return log(num)/log(2); }

// Round to the nearest integer
inline int iround(double x) { return static_cast<int>(floor(x+0.5)); }

// Returns pointer to first non-white-space character in str OR to 0 if none found
inline const char* strscn(const char* str)
{
    if (!str) return 0;
    const char* ptr=str;
    while (*ptr!='\0' && isspace(*ptr)) ptr++;
    return (*ptr=='\0') ? 0 : ptr;
}

// Returns leftmost integer in ptr and sets the pointer to first char after
// the integer. If no integer is found, returns INT_MIN and sets ptr to 0
inline int strtoi(const char*& ptr)
{
    const char* ptr0=ptr;
    if (!ptr) return INT_MIN;
    while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
    if (*ptr=='\0') {
        ptr=0;
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
        ptr=0;
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

#endif
