// Copyright 2009, Andreas Biegert

#ifndef SRC_GLOBALS_H_
#define SRC_GLOBALS_H_

namespace cs {

// Constants

#ifdef DEBUG
const bool kDebug = true;
#else
const bool kDebug = false;
#endif  // DEBUG

const int KB = 1024;
const int MB = KB * KB;
const int GB = KB * KB * KB;
const int kMaxInt = 0x7FFFFFFF;
const int kMinInt = -kMaxInt - 1;

const int kCharSize     = sizeof(char);    // NOLINT
const int kShortSize    = sizeof(short);   // NOLINT
const int kIntSize      = sizeof(int);     // NOLINT
const int kFloatSize    = sizeof(float);   // NOLINT
const int kDoubleSize   = sizeof(double);  // NOLINT
const int kPointerSize  = sizeof(void*);   // NOLINT


// A macro to disallow the evil copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName)      \
  TypeName(const TypeName&);                    \
  void operator=(const TypeName&)


// A macro to disallow all the implicit constructors, namely the
// default constructor, copy constructor and operator= functions.
//
// This should be used in the private: declarations for a class
// that wants to prevent anyone from instantiating it. This is
// especially useful for classes containing only static methods.
#define DISALLOW_IMPLICIT_CONSTRUCTORS(TypeName) \
  TypeName();                                    \
  DISALLOW_COPY_AND_ASSIGN(TypeName)

}  // namespace cs

#endif  // SRC_GLOBALS_H_



