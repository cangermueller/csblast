// Copyright 2009, Andreas Biegert

#ifndef SRC_LOG_H_
#define SRC_LOG_H_

#include <cstdio>
#include <sys/time.h>

#include <iostream>
#include <sstream>
#include <string>

enum LogLevel { ERROR, WARNING, INFO, DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4 };

class Log {
 public:
  Log();
  virtual ~Log();
  std::ostringstream& get(LogLevel level = INFO);

  static LogLevel& reporting_level();
  static std::string to_string(LogLevel log_level);
  static LogLevel from_string(const std::string& log_level);
  static LogLevel from_integer(int log_level);

 protected:
  std::ostringstream os;
  LogLevel level;

 private:
  Log(const Log&);
  Log& operator =(const Log&);
};



#ifndef LOG_MAX_LEVEL
#define LOG_MAX_LEVEL DEBUG4
#endif

#define LOG(level)                              \
  if (level > LOG_MAX_LEVEL) ;                  \
  else if (level > Log::reporting_level()) ;    \
  else Log().get(level)

#endif  // SRC_LOG_H_
