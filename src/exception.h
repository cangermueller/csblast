#ifndef CS_EXCEPTION_H
#define CS_EXCEPTION_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//This is a simple but useful class for construction of standard exceptions
//The constructor is implemented to take variadic arguments
//in the style of printf for simple message construction

#include <exception>
#include <string>
#include <cstdarg>

namespace cs
{

class Exception : public std::exception
{
public:
    Exception(const std::string &m) : msg(m), c(0) {}

    Exception(const char *str, ...) : c(0)
    {
        char *buffer = new char[1000];
        va_list ap;
        va_start(ap, str);
        vsprintf(buffer, str, ap);
        va_end(ap);
        msg = buffer;
        delete [] buffer;
    }

    Exception(const int type, const char *str, ...) : c(type)
    {
        char *buffer = new char[1000];
        va_list ap;
        va_start(ap, str);
        vsprintf(buffer, str, ap);
        va_end(ap);
        msg = buffer;
        delete [] buffer;
    }

    virtual ~Exception() throw() {};

    virtual const char* what() const throw()
    { return msg.c_str(); }

    virtual int code() const
    { return c; }

private:
    std::string msg;
    int c;
};

}  // cs

#endif
