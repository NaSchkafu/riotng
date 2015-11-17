/*

 File: vector.h, 2005/01/21

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

#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED

#include <iostream>

class Vector
{
 public:

  Vector() : row_(0), _M_(0) {}
  Vector(unsigned int row) : row_(row), _M_(0)
    {
      _M_ = new double[row];

      for(unsigned int i = 0; i < row; i++) _M_[i] = 0.0;
    }
  Vector(const Vector& A) : row_(A.row_), _M_(0)
    {
      _M_ = new double[A.row_];
      
      for(unsigned int i = 0; i < A.row_; i++) _M_[i] = A._M_[i];
    }
  ~Vector() { delete [] _M_; }

  Vector& operator = (const Vector& A)
    {
      if( &A == this ) return *this;

      row_ = A.row_;

      delete [] _M_;

      _M_ = new double[row_];
      
      for(unsigned int i = 0; i < row_; i++) _M_[i] = A._M_[i];
      
      return *this;
    }

  Vector& operator = (const double& d)
    {
      for(unsigned int i = 0; i < row_; i++) _M_[i] = d;
      return *this;
    }

  double&       operator [] (unsigned int i)       { return _M_[i]; }
  const double& operator [] (unsigned int i) const { return _M_[i]; }

  unsigned int Lb()   const { return      0; }
  unsigned int Ub()   const { return row_-1; }
  unsigned int size() const { return   row_; }

  double max() const;
  double min() const;

  /*
    Friend functions.
  */

  friend std::ostream& operator << (std::ostream&, const Vector&);

 private:

  unsigned int row_;
  double *_M_;
};

inline double max(const Vector& v) { return v.max(); }

#endif

/*

End of File: vector.h

*/
