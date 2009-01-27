/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "log.h"

#include <cstdio>
#include <sys/time.h>

#include <iostream>
#include <sstream>
#include <string>

namespace
{

inline std::string now_time()
{
    char buffer[11];
    time_t now = time(0);
    strftime(buffer, sizeof(buffer), "%X", localtime(&now));
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[100] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
    return result;
}

}  // namespace


namespace cs
{

Log::Log()
{}

std::ostringstream& Log::get(LogLevel level)
{
    log_level = level;
    os << now_time() << " " << to_string(level) << std::string(9 - to_string(level).length(), ' ');
    os << std::string(level > DEBUG ? (level - DEBUG) * 4 : 0, ' ');
    return os;
}

Log::~Log()
{
    if (os.str().find('\n') == std::string::npos)  // no newline in log output
        std::cerr << os.str() << std::endl;
    else {                                         // multiline log output -> insert left margin
        std::string s(os.str());
        const std::string margin(std::string(kLeftMargin, ' ') +
                                 std::string(log_level > DEBUG ? (log_level - DEBUG) * 4 : 0, ' '));

        if (s[s.length()-1] == '\n') s.erase(s.length()-1, 1);  // remove trailing newline

        std::string::size_type i = 0;
        std::string::size_type j = s.find('\n');

        while (j != std::string::npos) {
            std::cerr << s.substr(i, j-i) << std::endl << margin;
            i = ++j;
            j = s.find('\n', j);

            if (j == std::string::npos)
                std::cerr << s.substr(i, s.length()) << std::endl;
        }
    }
    std::cerr.flush();
}

LogLevel& Log::reporting_level()
{
    static LogLevel level = DEBUG4;
    return level;
}

std::string Log::to_string(LogLevel level)
{
    static const char* const buffer[] = {"ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
    return buffer[level];
}

LogLevel Log::from_string(const std::string& level)
{
    if (level == "DEBUG4")
        return DEBUG4;
    if (level == "DEBUG3")
        return DEBUG3;
    if (level == "DEBUG2")
        return DEBUG2;
    if (level == "DEBUG1")
        return DEBUG1;
    if (level == "DEBUG")
        return DEBUG;
    if (level == "INFO")
        return INFO;
    if (level == "WARNING")
        return WARNING;
    if (level == "ERROR")
        return ERROR;
    Log().get(WARNING) << "Unknown logging level '" << level << "'. Using INFO level as default.";
    return INFO;
}

}  // cs
