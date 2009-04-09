// Copyright 2009, Andreas Biegert
// Collection of commonly used inline utility functions.

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <cctype>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cmath>
#include <cstdarg>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cs {

// Returns the base 2 logarithm of x.
inline float log2(float x) {
  return 1.442695041 * log(x);
}

// Round to the nearest integer.
inline int iround(double x) {
  return static_cast<int>(floor(x + 0.5));
}

// Normalize a float array such that it sums to one. If it sums to 0 then assign
//  def_array elements to array (optional)
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

// Reset all entries in a float array to given value or zero if none provided.
inline void reset(float* array, int length, float value = 0.0f) {
  for (int i = 0; i < length; ++i) array[i] = value;
}

// Splits a string into multiple strings with character delimiter removed.
inline void split(const std::string& s, char c, std::vector<std::string>& v) {
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
  char sep = '/';

#ifdef _WIN32
  sep = '\\';
#endif

  size_t i = s.rfind(sep, s.length());
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
  char sep = '/';

#ifdef _WIN32
  sep = '\\';
#endif

  size_t i = s.rfind(sep, s.length( ));
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
