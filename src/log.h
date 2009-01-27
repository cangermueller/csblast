#ifndef CS_LOG_H
#define CS_LOG_H

#include <cstdio>
#include <sys/time.h>

#include <iostream>
#include <sstream>
#include <string>

namespace cs
{

enum LogLevel {ERROR, WARNING, INFO, DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4};

class Log
{
  public:
    Log();
    virtual ~Log();
    std::ostringstream& get(LogLevel level = INFO);

    static LogLevel& reporting_level();
    static std::string to_string(LogLevel log_level);
    static LogLevel from_string(const std::string& log_level);

  protected:
    static const int kMargin = 22;
    static const int kIndent = 4;

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
    if (level > LOG_MAX_LEVEL) ;                \
    else if (level > Log::reporting_level()) ;   \
    else Log().get(level)

}  // cs

#endif
