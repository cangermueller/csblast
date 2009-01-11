#ifndef SMART_PTR_H
#define SMART_PTR_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// The SmartPtr class template stores a pointer to a dynamically allocated
// object, typically with a C++ new-expression. The object pointed to is
// guaranteed to be deleted when the last SmartPtr pointing to it is destroyed
// or reset.

template <class T>
class SmartPtr
{
public:
    SmartPtr(T *p) : p_(p), c_(new long(1)) {}
    ~SmartPtr() { if(!--*c_) { delete c_; delete p_; } }
    SmartPtr(const SmartPtr &init) : p_(init.p_), c_(init.c_) { ++*c_; }

    SmartPtr &operator =(const SmartPtr &rhs)
    {
        if (this != &rhs) {
            if (!--*c_) { delete c_; delete p_; }
            p_ = rhs.p_;
            ++*(c_ = rhs.c_);
        }
        return *this;
    }
    T &operator *() const { return *p_; }
    T &operator ->() const { return p_; }

private:
    T *p_;
    long *c_;
};

#endif
