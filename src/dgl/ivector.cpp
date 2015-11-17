/*

 File: ivector.cpp, 2003/09/29

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

#include "ivector.h"

Vector IVector::sup() const
{
  Vector res(row_);
  for(unsigned int i = 0; i < row_; i++) res[i] = _M_[i].sup();
  return res;
}
Vector IVector::inf() const
{
  Vector res(row_);
  for(unsigned int i = 0; i < row_; i++) res[i] = _M_[i].inf();
  return res;
}
Vector IVector::mid() const
{
  Vector res(row_);
  for(unsigned int i = 0; i < row_; i++) res[i] = _M_[i].mid();
  return res;
}

IVector abs(const IVector& v)
{
  IVector res(v.row_);
  for(unsigned int i = 0; i < v.row_; i++) res[i] = v._M_[i].abs();
  return res;
}

std::ostream& operator << (std::ostream& os, const IVector& A)
{
  for(unsigned int i = 0; i < A.row_; i++)
    {
      os << A._M_[i] << " ";
    }
  os << std::endl;
  
  return os;
}


/*

Other operators.

*/

IVector operator * (const Matrix& A, const IVector& x)
{
  IVector res(A.size(1));

  if( A.size(2) == x.size() )
    {
      for(unsigned int i = 0; i <= res.Ub(); i++) 
	for(unsigned int j = 0; j <= x.Ub(); j++) 
	  res[i] += A(i,j) * x[j];
    }
  else
    std::cout << "Matrix * IVector: *** Dimensions must agree ***" << std::endl;

  return res;
}

IVector operator * (const IMatrix& A, const IVector& x)
{
  IVector res(A.size(1));

  if( A.size(2) == x.size() )
    {
      for(unsigned int i = 0; i <= res.Ub(); i++) 
	for(unsigned int j = 0; j <= x.Ub(); j++) 
	  res[i] += A(i,j) * x[j];
    }
  else
    std::cout << "IMatrix * IVector: *** Dimensions must agree ***" << std::endl;

  return res;
}


/*

  End of File: ivector.cpp

*/
