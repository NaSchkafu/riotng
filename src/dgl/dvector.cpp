/*

 File: dvector.cpp, 2005/01/21

 Copyright (C) Ingo Eble,    IngoEble@web.de

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

#include "dvector.h"

double Vector::max() const
{
  if( row_ > 0 )
    {
      double m = _M_[0];
      for(unsigned int i = 1; i < row_; i++) m = (m < _M_[i])?_M_[i]:m;
      return m;
    }
  else return 0.0;
}

double Vector::min() const
{
  if( row_ > 0 )
    {
      double m = _M_[0];
      for(unsigned int i = 1; i < row_; i++) m = (m > _M_[i])?_M_[i]:m;
      return m;
    }
  else return 0.0;
}

std::ostream& operator << (std::ostream& os, const Vector& A)
{
  for(unsigned int i = 0; i < A.row_; i++)
    {
      os << A._M_[i] << " ";
    }
  os << std::endl;
  
  return os;
}

/*

  End of File: vector.cpp

*/
