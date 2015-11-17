/*

  File: factorialtable.cpp, 2004/08/18

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

#include "factorialtable.h"
namespace riot
{

  Factorials::Factorials() : facs_( 200, Interval(1.0) ), odd_facs_( 200, Interval(1.0) )
  {
    for(int i = 2; i < 200; i++ )     facs_[  i] = i * facs_[i-1];
    for(int i = 3; i < 199; i+=2) odd_facs_[i+1] = (odd_facs_[i] = i * odd_facs_[i-1]);
  }
}

/*

  End of File: factorialtable.cpp

*/
