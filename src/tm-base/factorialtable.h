/*

  File: factorialtable.h, 2004/08/18

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

#ifndef FACTORIALTABLE_H_INCLUDED
#define FACTORIALTABLE_H_INCLUDED

#include <vector>
#include "adaptintval.h"

using Adapt::Interval;
namespace riot
{

  class Factorials
  {
  private:

    /*
      To make it a Singleton.
    */

    Factorials();
    Factorials(const Factorials&);
    ~Factorials() {}
    Factorials& operator = (const Factorials&);

  public:

    static Factorials& getTable()
    {
      static Factorials singleton;
      return singleton;
    }

    Interval factorial    (unsigned int i) const { return facs_[i];     }
    Interval odd_factorial(unsigned int i) const { return odd_facs_[i]; }

  private:

    /*
      Data.
    */

    std::vector<Interval> facs_, odd_facs_;
  };
}

#endif

/*

  End of File: factorialtable.h

*/
