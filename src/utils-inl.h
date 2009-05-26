// Copyright 2009, Andreas Biegert
// Collection of commonly used inline utility functions.

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <cassert>
#include <cctype>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <stdint.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "globals.h"

namespace cs {

// Returns the base 2 logarithm of x.
inline float log2(float x) {
  return 1.442695041f * log(x);
}

// Returns the base 2 logarithm of x.
inline double log2(double x) {
  return 1.442695041 * log(x);
}

// Round to the nearest integer.
inline int iround(double x) {
  return static_cast<int>(floor(x + 0.5));
}

// This function returns log2 with a max abolute deviation of +/- 1.5E-5
// It takes 1.42E-8 s  whereas log2(x) takes 9.5E-7 s. It is hence 9.4 times faster.
// It makes use of the representation of 4-byte floating point numbers:
// seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// s is the sign,
// the following 8 bits, eee eee e, give the exponent + 127 (in hex: 0x7f).
// The following 23 bits, m, give the mantisse, the binary digits behind the
// decimal point.
// In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
// The expression (((*(int *)&x) & 0x7f800000 ) >>23 )-0x7f is the exponent eeeeeeee,
// i.e. the largest integer that is smaller than log2(x) (e.g. -1 for 0.9). *(int *)&x
// is an integer which
// contains the bytes as the floating point variable x is represented in memory.
// Check:  assert( sizeof(f) == sizeof(int) );
// Check:  assert( sizeof(f) == 4 );
inline float fast_log2(float x) {
  if (x <= 0.0f) return log2(x);

  static float lg2[1025];   // lg2[i] = log2[1+x/1024]
  static float diff[1025];  // diff[i]= (lg2[i+1]-lg2[i])/8096 (for interpolation)
  static bool initialized;

  if (!initialized) {
      assert( sizeof(x) == kIntSize );
      assert( sizeof(x) == 4 );

      float prev = 0.0f;
      lg2[0] = 0.0f;
      for (int i = 1; i <= 1024; ++i) {
        lg2[i] = log(float(1024+i))*1.442695041-10.0f;
        diff[i-1] = (lg2[i]-prev)*1.2352E-4;
        prev = lg2[i];
      }
      initialized=true;
  }

  int a = (((*((int *)&x)) & 0x7F800000) >>23 )-0x7f;
  int b =  ((*((int *)&x)) & 0x007FE000) >>13;
  int c =  ((*((int *)&x)) & 0x00001FFF);

  return a + lg2[b] + diff[b]*(float)(c);
}

// Fast 2^x
// ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// Relative deviation < 2.3E-7
// Speed: 2.3E-8s per call! (exp(): 8.5E-8, pow(): 1.7E-7)
//                        seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
inline float fast_pow2(float x) {
  if (x == -INFINITY || x > FLT_MAX_EXP || x < FLT_MIN_EXP)
    return pow(2.0f, x);

  // store address of float as pointer to long
  int *px = (int*)(&x);
  // temporary value for truncation: x-0.5 is added to a large integer (3<<22)
  float tx = (x-0.5f) + (3<<22);
  int lx = *((int*)&tx) - 0x4b400000;   // integer value of x
  float dx = x-(float)(lx);             // float remainder of x
  x = 1.0f + dx*(0.693153f              // polynomial apporoximation of 2^x
           + dx*(0.240153f              // for x in the range [0, 1]
           + dx*(0.0558282f
           + dx*(0.00898898f
           + dx* 0.00187682f ))));

  *px += (lx<<23);                      // add integer power of 2 to exponent
  return x;
}

// /////////////////////////////////////////////////////////////////////////////////////
// // fast 2^x
// // ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// // Relative deviation < ??  (< 2.3E-7 with 5'th order polynomial)
// // Speed: ???E-8s (??E-8s) per call!
// // Internal representation of double number according to IEEE 754:
// //   1bit sign, 11 bits exponent, 52 bits mantissa: seee eeee eeee mmmm mmmm mmmm mmmm mmmm mmmm mmmm mmmm mmmm mmmm mmmm mmmm mmmm
// // In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm... * 2^(eeeeeeeeee-1023)
// /////////////////////////////////////////////////////////////////////////////////////
// inline double fast_pow2(double x) {
//   if (x == -INFINITY || x > DBL_MAX_EXP || x < DBL_MIN_EXP)
//     return pow(2.0, x);
//   uint64_t *px = (uint64_t*)(&x);       // store address of float as pointer to long int
//   double tx = (x - 0.5) + (static_cast<uint64_t>(3) << 51);        // temporary value for truncation: x-0.5 is added to a large integer (3<<51)
//                                         // 3<<51 = (1.1bin)*252 = (1.1bin)*2^(1075-1023)
//                                         // which, in internal bits, is written 0x4338000000000000 (since 0100 0011 0011 bin = 1075)
//   const uint64_t kTMP = (static_cast<uint64_t>(0x43380000) << 32) | 0x00000000;
//   uint64_t lx = *((uint64_t*)&tx) - kTMP;   // integer value of x
//   //  uint64_t lx = *((uint64_t*)&tx) -             0x4338000000000000ull;   // integer value of x
//   double dx = x-(double)(lx);           // float remainder of x
//   // x = 1.0f + dx*(0.693019d              // polynomial apporoximation of 2^x
//   //          + dx*(0.241404d              // for x in the range [0, 1]
//   //          + dx*(0.0520749d
//   //          + dx* 0.0134929d )));
//   x = 1.0 + dx*(0.693153               // polynomial apporoximation of 2^x
//            + dx*(0.240153              // for x in the range [0, 1]
//            + dx*(0.0558282
//            + dx*(0.00898898
//            + dx* 0.00187682 ))));

//   *px += ((static_cast<uint64_t>(lx) << 32) | 0x00000000);
//   //  *px += (static_cast<uint64_t>(lx) << 52);                      // add integer power of 2 to exponent
//   return x;
// }

// Fast 2^x
// ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// Relative deviation < 2.3E-7
// Speed: 2.3E-8s per call! (exp(): 8.5E-8, pow(): 1.7E-7)
//                        seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
inline double fast_pow2(double d) {
  if (d == -INFINITY || d > FLT_MAX_EXP || d < FLT_MIN_EXP)
    return pow(2.0, d);

  float x = d;
  // store address of float as pointer to long
  int *px = (int*)(&x);
  // temporary value for truncation: x-0.5 is added to a large integer (3<<22)
  float tx = (x-0.5f) + (3<<22);
  int lx = *((int*)&tx) - 0x4b400000;   // integer value of x
  float dx = x-(float)(lx);             // float remainder of x
  x = 1.0f + dx*(0.693153f              // polynomial apporoximation of 2^x
           + dx*(0.240153f              // for x in the range [0, 1]
           + dx*(0.0558282f
           + dx*(0.00898898f
           + dx* 0.00187682f ))));

  *px += (lx<<23);                      // add integer power of 2 to exponent
  return x;
}

// Fast 2^x
// ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// Relative deviation < 2.3E-7
// Speed: 2.3E-8s per call! (exp(): 8.5E-8, pow(): 1.7E-7)
//                        seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
// inline double fast_pow2(double d) {
//   if (d == -INFINITY || d > DBL_MAX_EXP || d < DBL_MIN_EXP)
//     return pow(2.0, d);

//   double y = 1.0;
//   while (d > FLT_MAX_EXP) {
//     d -= FLT_MAX_EXP;
//     y *= FLT_MAX;
//   }
//   while (d < FLT_MIN_EXP) {
//     d += FLT_MAX_EXP;
//     y /= FLT_MAX;
//   }
//   float x = d;

//   // store address of float as pointer to long
//   int *px = (int*)(&x);
//   // temporary value for truncation: x-0.5 is added to a large integer (3<<22)
//   float tx = (x-0.5f) + (3<<22);
//   int lx = *((int*)&tx) - 0x4b400000;   // integer value of x
//   float dx = x-(float)(lx);             // float remainder of x
//   x = 1.0f + dx*(0.693153f              // polynomial apporoximation of 2^x
//            + dx*(0.240153f              // for x in the range [0, 1]
//            + dx*(0.0558282f
//            + dx*(0.00898898f
//            + dx* 0.00187682f ))));

//   *px += (lx<<23);                      // add integer power of 2 to exponent
//   return y * x;
// }

// Normalize a float array such that it sums to one. If it sums to 0 then assign
// def_array elements to array (optional)
inline float normalize_to_one(float* array, int length,
                              const float* default_array = NULL) {
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

// Normalize a double array such that it sums to one. If it sums to 0 then assign
// def_array elements to array (optional)
inline double normalize_to_one(double* array, int length,
                               const double* default_array = NULL) {
  double sum = 0.0;
  for (int i = 0; i < length; ++i) sum += array[i];
  if (sum != 0.0) {
    double fac = 1.0 / sum;
    for (int i = 0; i < length; ++i) array[i] *= fac;
  } else if (default_array) {
    for (int i = 0; i < length; ++i) array[i] = default_array[i];
  }
  return sum;
}

// Reset all entries in a float array to given value or zero if none provided.
inline void Reset(float* array, int length, float value = 0.0f) {
  for (int i = 0; i < length; ++i) array[i] = value;
}

// Splits a string into multiple strings with character delimiter removed.
inline void split(const std::string& s, char c, std::vector<std::string>* v) {
  std::string::size_type i = 0;
  std::string::size_type j = s.find(c);

  while (j != std::string::npos) {
    v->push_back(s.substr(i, j-i));
    i = ++j;
    j = s.find(c, j);

    if (j == std::string::npos)
      v->push_back(s.substr(i, s.length()));
  }
}

// Stringifies the range elemeents delimited by given character.
template<typename Fwd>
inline std::string stringify_range(Fwd first, Fwd last, char delim = ',') {
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
inline std::string stringify_container(const C& c, char delim = ',') {
  return stringify_range(c.begin(), c.end(), delim);
}

// sprintf-like helper that returns a string.
inline std::string strprintf(const char* str, ...) {
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
inline std::string get_file_ext(const std::string& s) {
  size_t i = s.rfind('.', s.length());
  if (i != std::string::npos) {
    return(s.substr(i+1, s.length() - i));
  }

  return "";
}

// Returns the last component of the filename. If ext is false
// the file extension is removed.
inline std::string get_file_basename(const std::string& s, bool ext = true) {
  size_t i = s.rfind(kDirSep, s.length());
  if (i != std::string::npos) {
    if (ext) return(s.substr(i+1, s.length() - i));

    size_t j = s.rfind('.', s.length());
    if (j != std::string::npos && j > i)
      return(s.substr(i+1, j-i));
    else
      return(s.substr(i+1, s.length() - i));
  } else if (!ext) {
    size_t j = s.rfind('.', s.length());
    if (j != std::string::npos)
      return(s.substr(0, j+1));
    else
      return s;
  }
  return s;
}

// Returns all components of the filename except the last one.
inline std::string get_file_dirname(const std::string& s) {
  size_t i = s.rfind(kDirSep, s.length( ));
  if (i != std::string::npos) {
    return(s.substr(0, i));
  }

  return "";
}

// Removes the newline and other control characters at the end of a string (if
//  present) and returns the new length of the string (-1 if str is NULL)
inline int chomp(char* str) {
  if (!str) return -1;
  int l = 0;
  for (l = strlen(str) - 1; l >= 0 && str[l] < 32; --l)
    /* do nothing */;
  str[++l] = '\0';
  return l;
}

// Emulates the ifstream::getline method; similar to fgets(str,maxlen,FILE*),
// but removes the newline at the end and returns NULL if at end of file or read
// error
inline char* fgetline(char* str, int maxlen, FILE* file) {
  if (!fgets(str, maxlen, file)) return NULL;
  if (chomp(str) + 1 >= maxlen)  // if line is cut after maxlen characters...
    while (fgetc(file) != '\n')  // ... read in rest of line
      /* do nothing */;
  return(str);
}

// Returns leftmost integer in ptr and sets the pointer to first char after
// the integer. If no integer is found, returns INT_MIN and sets ptr to NULL
inline int strtoi(const char*& ptr) {
  int i;
  const char* ptr0 = ptr;
  if (!ptr) return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9')) ++ptr;
  if (*ptr == '\0') {
    ptr = 0;
    return INT_MIN;
  }
  if (ptr > ptr0 && *(ptr-1) == '-') i = -atoi(ptr); else i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9') ++ptr;
  return i;
}

// Same as strtoi, but interpretes '*' as default
inline int strtoi_ast(const char*& ptr, int deflt = INT_MAX) {
  int i;
  const char* ptr0 = ptr;
  if (!ptr) return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9') && *ptr != '*') ++ptr;
  if (*ptr == '\0') {
    ptr = 0;
    return INT_MIN;
  }
  if (*ptr == '*') {
    ++ptr;
    return deflt;
  }
  if (ptr > ptr0 &&  *(ptr-1) == '-') i = -atoi(ptr); else i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9') ++ptr;
  return i;
}

// Returns pointer to first non-white-space character in str OR to NULL if none
// found
inline const char* strscn(const char* str) {
  if (!str) return NULL;
  const char* ptr = str;
  while (*ptr != '\0' && isspace(*ptr)) ++ptr;
  return (*ptr == '\0') ? NULL : ptr;
}

}  // namespace cs

#endif  // SRC_UTILS_H_
