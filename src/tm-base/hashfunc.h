/*

  File: hashfunc.h, 2004/08/18

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

#ifndef HASHFUNC_H_INCLUDED
#define HASHFUNC_H_INCLUDED

#include "primetable.h"
#include "hashtable.h"
#include "monom.h"

namespace riot
{

/*-----------------------------------------------------------------------------+
  |        Here we define Hash functions for keys of type Monom.                |
  +-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------+
  | The implementation of the following hash-function is taken from the         |
  | "spock"-Package, a package to calculate with sparse polynomials. See        |
  | http://www.mathematik.tu-muenchen.de/~kaplan/ca/spock/                      |
  |                                                                             |
  | Copyright (c) 1997-98  spock-developers and                                 |
  |                        Technische Universitaet Muenchen, Munich (Germany)   |
  +-----------------------------------------------------------------------------*/

  class MonomHasher1
  {
  private:

    /*
      To make it a Singleton.
    */

    MonomHasher1() {}
    MonomHasher1(const MonomHasher1&);
    ~MonomHasher1() {}
    MonomHasher1& operator = (const MonomHasher1&);

  public:

    static MonomHasher1& getFunction()
    {
      static MonomHasher1 singleton;
      return singleton;
    }

    long unsigned int operator ()(const Monom& m, long unsigned int size) const
    {
      long unsigned int h = 0;

      for(unsigned int i = m.min_var_code(); i <= m.max_var_code(); i++)
      {
        h = ( h*h + i * m.expnt_of_var( i ) ) % size;
      }

      return h % size;
    }
  };

/*-----------------------------------------------------------------------------+
  | Another implementation (selfmade)                                           |
  +-----------------------------------------------------------------------------*/

  class MonomHasher2
  {
  private:

    /*
      To make it a Singleton.
    */

  MonomHasher2() : base_( PrimeNumbers::getTable().find_GreatestPrimeSmallerThan(127) )
    {
      PrimeNumbers::getTable().lock( base_ );
    }
    MonomHasher2(const MonomHasher2&);
    ~MonomHasher2() {}
    MonomHasher2& operator = (const MonomHasher2&);

  public:

    static MonomHasher2& getFunction()
    {
      static MonomHasher2 singleton;
      return singleton;
    }

    long unsigned int operator ()(const Monom& m, long unsigned int size) const
    {
      long unsigned int h = m.expnt_of_var(8);

      for(unsigned int i = 7; i > 0; i--)
      {
        h *= base_;
        h += m.expnt_of_var( i );
        h %= size;
      }

      return h % size;
    }

  private:

    /*
      Data.
    */

    long unsigned int base_;
  };
}

#endif

/*

  End of File: hashfunc.h

*/
