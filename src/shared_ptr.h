#ifndef SMART_PTR_H
#define SMART_PTR_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// The shared_ptr class template stores a pointer to a dynamically allocated
// object, typically with a C++ new-expression. The object pointed to is
// guaranteed to be deleted when the last shared_ptr pointing to it is destroyed
// or reset.

template <class T>
class shared_ptr
{
  public:
    shared_ptr(T* p) : p_(p), c_(new long(1)) {}
    ~shared_ptr() { if(!--*c_) { delete c_; delete p_; } }
    shared_ptr(const shared_ptr &init) : p_(init.p_), c_(init.c_) { ++*c_; }

    shared_ptr &operator =(const shared_ptr &rhs)
    {
        if (this != &rhs) {
            if (!--*c_) { delete c_; delete p_; }
            p_ = rhs.p_;
            ++*(c_ = rhs.c_);
        }
        return *this;
    }
    T& operator *() const { return *p_; }
    T* operator ->() const { return p_; }

    // Template function for implicit conversion (see More Effectove C++ page 176 for details).
    template <class newType>
    operator shared_ptr<newType>() { return shared_ptr<newType>(p_); }

  private:
    T *p_;
    long *c_;
};

#endif
