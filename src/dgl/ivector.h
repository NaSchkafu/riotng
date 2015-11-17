/*

 File: ivector.h, 2004/08/05

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

#ifndef IVECTOR_H_INCLUDED
#define IVECTOR_H_INCLUDED

#include "adaptintval.h"

#include "dvector.h"
#include "matrix.h"
#include "imatrix.h"

using Adapt::Interval;

class IVector
{
 public:

  IVector() : row_(0), _M_(0) {}
  IVector(unsigned int row) : row_(row), _M_(0)
    {
      _M_ = new Interval[row_];

      for(unsigned int i = 0; i < row; i++) _M_[i] = Adapt::ZERO_INTERVAL();
    }
  IVector(const IVector& A) : row_(A.row_), _M_(0)
    {
      _M_ = new Interval[row_];
      
      for(unsigned int i = 0; i < row_; i++) _M_[i] = A._M_[i];
    }
  ~IVector() { delete [] _M_; }

  IVector& operator = (const IVector& A)
    {
      if( &A == this ) return *this;

      row_ = A.row_;

      delete [] _M_;

      _M_ = new Interval[row_];
      
      for(unsigned int i = 0; i < row_; i++) _M_[i] = A._M_[i];

      return *this;
    }

  IVector& operator += (const IVector& A)
    {
      if( row_ == A.row_ )
	{
	  for(unsigned int i = 0; i < row_; i++) _M_[i] += A._M_[i];
	}

      return *this;
    }
  IVector& operator -= (const IVector& A)
    {
      if( row_ == A.row_ )
	{
	  for(unsigned int i = 0; i < row_; i++) _M_[i] -= A._M_[i];
	}

      return *this;
    }

  Interval&       operator [] (unsigned int i)       { return _M_[i]; }
  const Interval& operator [] (unsigned int i) const { return _M_[i]; }

  unsigned int Lb()   const { return      0; }
  unsigned int Ub()   const { return row_-1; }
  unsigned int size() const { return   row_; }

  Vector sup() const;
  Vector inf() const;
  Vector mid() const;

  /*

  Friend functions.

  */

  friend IVector abs(const IVector&);

  friend std::ostream& operator << (std::ostream&, const IVector&);

 private:

  unsigned int row_;
  Interval *_M_;
};

inline Vector sup(const IVector& v) { return v.sup(); }
inline Vector inf(const IVector& v) { return v.inf(); }
inline Vector mid(const IVector& v) { return v.mid(); }

IVector operator * (const  Matrix&,const IVector&);
IVector operator * (const IMatrix&,const IVector&);

#endif

/*

End of File: ivector.h

*/
