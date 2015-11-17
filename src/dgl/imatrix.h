/*

 File: imatrix.h, 2003/09/23

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

#ifndef IMATRIX_H_INCLUDED
#define IMATRIX_H_INCLUDED

#include "matrix.h"
#include "adaptintval.h"

using Adapt::Interval;

class IMatrix
{
 public:

  IMatrix() : row_(0), col_(0), _M_(0) {}
  IMatrix(int row, int col) : row_(row), col_(col), _M_(0)
    {
      _M_ = new Interval[row*col];

      for(unsigned int i = 0; i < row*col; i++) _M_[i] = Adapt::ZERO_INTERVAL();
    }
  IMatrix(const IMatrix& A) : row_(A.row_), col_(A.col_), _M_(0)
    {
      _M_ = new Interval[A.row_*A.col_];
      
      for(int i = 0; i < A.row_*A.col_; i++) _M_[i] = A._M_[i];
    }
  ~IMatrix() { delete [] _M_; }

  IMatrix& operator = (const IMatrix& A)
    {
      if( &A == this ) return *this;

      row_ = A.row_;
      col_ = A.col_;

      delete [] _M_;

      _M_ = new Interval[row_*col_];
      
      for(int i = 0; i < row_*col_; i++) _M_[i] = A._M_[i];

      return *this;
    }

  Interval& operator () (int i, int j) { return _M_[i*row_+j]; }
  const Interval& operator () (int i, int j) const { return _M_[i*row_+j]; }

  int Lb(int i) const { return 0; }
  int Ub(int i) const 
    { 
      switch(i)
	{
	case 1: return row_-1; 
	case 2: return col_-1; 
	}
      return 0;
    }
  unsigned int size(unsigned int i) const
    {
      switch(i)
	{
	case 1: return row_; 
	case 2: return col_; 
	}
      return 0;
    }

  friend std::ostream& operator << (std::ostream& os, const IMatrix& A)
    {
      for(int i = 0; i < A.row_; i++)
	{
	  for(int j = 0; j < A.col_; j++) os << A._M_[i*A.row_+j] << " ";
	  os << std::endl;
	}

      return os;
    }

 private:

  int row_, col_;
  Interval *_M_;
};

void imatinv(const Matrix&, IMatrix&, int&);

#endif

/*

End of File: imatrix.h

*/
