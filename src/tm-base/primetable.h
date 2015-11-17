/*

  File: primetable.h, 2004/08/18

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

#ifndef PRIMETABLE_H_INCLUDED
#define PRIMETABLE_H_INCLUDED

#include <vector>

namespace riot
{

/*-----------------------------------------------------------------------------+
  |                  Table of Prime numbers as 'Singleton'                      |
  +-----------------------------------------------------------------------------*/

  class PrimeNumbers
  {
  private:

    /*
      To make it a Singleton.
    */

    PrimeNumbers() {}
    PrimeNumbers(unsigned long);
    PrimeNumbers(const PrimeNumbers&);
    ~PrimeNumbers() {}

    PrimeNumbers& operator = (const PrimeNumbers&);

  public:

    static PrimeNumbers& getTable()
    {
      static PrimeNumbers singleton(1000000);
      return singleton;
    }

    void lock(unsigned long l) { locked_prime_ = l; }

    unsigned long find_SmallestPrimeGreaterThan(const double&) const;
    unsigned long find_GreatestPrimeSmallerThan(const double&) const;

  private:

    /*
      Data.
    */

    std::vector<unsigned long> table_;
    unsigned long              locked_prime_;
  };
}

#endif

/*

  End of File: primetable.h

*/
