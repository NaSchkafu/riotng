/*

  File: refcounter.h, 2004/08/18

  Copyright (C) Ingo Eble,    ingoeble@web.de

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef REFCOUNTER_H_INCLUDED
#define REFCOUNTER_H_INCLUDED

/*

  +---------+
  | :Handle |
  +-+-----+-+            +----+
  | :RC |------------->| :T |
  +-----+\             +----+
  \               |
  \           +---------+
  ---------->| counter |
  +---------+
*/
namespace riot
{

  template<class T>
    class ReferenceCounter
  {
  public:

    /*
      Constructors
    */

  ReferenceCounter() : objptr_((T*)(0)), counter_(new unsigned int(1)) {}
  ReferenceCounter(T *ptr) : objptr_(ptr), counter_(new unsigned int(1)) {}
  ReferenceCounter(T *ptr, unsigned int *cnt) : objptr_(ptr), counter_(cnt) { ++(*counter_); }
  ReferenceCounter(const ReferenceCounter& x)
    : objptr_(x.objptr_), counter_(x.counter_)
    {
      ++(*counter_);
    }

    ~ReferenceCounter()
    {
      if( --(*counter_) == 0 )
      {
        delete objptr_;
        delete counter_;
      }
    }
    ReferenceCounter& operator = (const ReferenceCounter& x)
      {
        ++(*(x.counter_));

        if( --(*counter_) == 0 )
        {
          delete objptr_;
          delete counter_;
        }

        objptr_  = x.objptr_;
        counter_ = x.counter_;

        return *this;
      }

    /*
      Operators for comparison
    */

    bool operator >  (unsigned int i) const { return *counter_ >  i; }
    bool operator <  (unsigned int i) const { return *counter_ <  i; }
    bool operator == (unsigned int i) const { return *counter_ == i; }

    /*
      Access to object
    */

    T* operator -> () { return  objptr_; }
    T& operator *  () { return *objptr_; }
    const T* operator -> () const { return  objptr_; }
    const T& operator *  () const { return *objptr_; }

  private:

    T            *objptr_;
    unsigned int *counter_;
  };
}

#endif

/*

  End of File: refcounter.h

*/
