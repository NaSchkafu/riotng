/*

 File: matrix.h, 2005/02/20

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

#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>

class Matrix
{
 public:

  Matrix() : row_(0), col_(0), _M_(0) {}
  Matrix(unsigned int row, unsigned int col) : row_(row), col_(col), _M_(0)
    {
      _M_ = new double[row*col];

      for(unsigned int i = 0; i < row*col; i++) _M_[i] = 0.0;
    }
  Matrix(const Matrix& A) : row_(A.row_), col_(A.col_), _M_(0)
    {
      _M_ = new double[A.row_*A.col_];
      
      for(unsigned int i = 0; i < A.row_*A.col_; i++) _M_[i] = A._M_[i];
    }
  ~Matrix() { delete [] _M_; }

  Matrix& operator = (const Matrix& A)
    {
      if( &A == this ) return *this;

      row_ = A.row_;
      col_ = A.col_;

      delete [] _M_;

      _M_ = new double[row_*col_];
      
      for(unsigned int i = 0; i < row_*col_; i++) _M_[i] = A._M_[i];

      return *this;
    }

  Matrix& operator = (const double& d)
    {
      for(unsigned int i = 0; i < row_*col_; i++) _M_[i] = d;
      return *this;
    }

  double&       operator () (unsigned int i, unsigned int j)       { return _M_[i*col_+j]; }
  const double& operator () (unsigned int i, unsigned int j) const { return _M_[i*col_+j]; }

  unsigned int Lb(unsigned int i) const { return 0; }
  unsigned int Ub(unsigned int i) const 
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

  double max() const;
  double min() const;

  friend std::ostream& operator << (std::ostream& os, const Matrix& A)
    {
      for(unsigned int i = 0; i < A.row_; i++)
	{
	  for(unsigned int j = 0; j < A.col_; j++) os << A._M_[i*A.col_+j] << " ";
	  os << std::endl;
	}

      return os;
    }

 private:

  unsigned int row_, col_;
  double *_M_;
};

inline double max(const Matrix& A) { return A.max(); }
inline double min(const Matrix& A) { return A.min(); }

/*
  Computes the determinant if an (n x n)-Matrix.
*/

double det(const Matrix&,int&);

/*
  Computes absolute values of elements.
*/

Matrix abs(const Matrix&);

#endif

/*

End of File: matrix.h

*/
