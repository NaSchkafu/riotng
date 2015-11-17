/*

 File: vvtaylormodel.h, 2003/10/08

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

#ifndef VVTAYLORMODEL_H_INCLUDED
#define VVTAYLORMODEL_H_INCLUDED

#include "taylormodel.h"
#include "dvector.h"
#include "ivector.h"
#include "matrix.h"
#include "imatrix.h"

namespace riot {
using Adapt::Interval;

class VVTaylorModel
{
 public:

  VVTaylorModel() : row_(0), _M_(0) {}
  VVTaylorModel(unsigned int dim) : row_(dim), _M_(0) 
    {
      _M_ = new TaylorModel[dim];

      for(unsigned int i = 0; i < dim; i++) _M_[i] = TaylorModel::ZERO_TM();
    }

  /*
    Check of dimensions for arguments is still missing!
  */

  VVTaylorModel(unsigned int dim,
		const std::vector<std::string>& vname, 
		const std::vector<double>&      vdpoint, 
		const std::vector<Interval>&    vdomain ) 
    : row_(dim), _M_(0)
    {
      _M_ = new TaylorModel[dim];

      for(unsigned int i = 0; i < dim; i++) _M_[i] = TaylorModel(vname[i],vdpoint[i],vdomain[i]);
    }
  VVTaylorModel(unsigned int dim, 
		const std::vector<std::string>& vname, 
		const Vector&                   vdpoint, 
		const IVector&                  vdomain ) 
    : row_(dim), _M_(0)
    {
      _M_ = new TaylorModel[dim];

      for(unsigned int i = 0; i < dim; i++) _M_[i] = TaylorModel(vname[i],vdpoint[i],vdomain[i]);
    }

  VVTaylorModel(const VVTaylorModel& A) : row_(A.row_)
    {
      _M_ = new TaylorModel[A.row_];

      for(unsigned int i = 0; i < A.row_; i++) _M_[i] = A._M_[i];
    }

  VVTaylorModel& operator = (const VVTaylorModel& A)
    {
      if( &A == this ) return *this;

      row_ = A.row_;

      delete [] _M_;

      _M_ = new TaylorModel[A.row_];

      if( A._M_ ) for(unsigned int i = 0; i < A.row_; i++) _M_[i] = A._M_[i];

      return *this;
    }

  ~VVTaylorModel() { delete [] _M_; }

  /*
    Read access.
  */

  unsigned int size() const { return row_; }
  IVector interval_part() const;
  Vector  constant_part() const;
  IVector eval()          const;
  const TaylorModel& operator [] (unsigned int i) const { return _M_[i]; }

  /*
    Write access.
  */

  TaylorModel& operator [] (unsigned int i) { return _M_[i]; }

  VVTaylorModel  operator - () const;

  VVTaylorModel& operator += (const VVTaylorModel&);
  VVTaylorModel& operator -= (const VVTaylorModel&);

  VVTaylorModel& operator += (const TaylorModel&);
  VVTaylorModel& operator -= (const TaylorModel&);
  VVTaylorModel& operator *= (const TaylorModel&);
  VVTaylorModel& operator /= (const TaylorModel&);

  VVTaylorModel& operator += (const Interval&);
  VVTaylorModel& operator -= (const Interval&);
  VVTaylorModel& operator *= (const Interval&);
  VVTaylorModel& operator /= (const Interval&);

  VVTaylorModel& operator += (const double&);
  VVTaylorModel& operator -= (const double&);
  VVTaylorModel& operator *= (const double&);
  VVTaylorModel& operator /= (const double&);

  void add_to_interval_part(const IVector&); 
  void replace_interval_part_with(const IVector&);

  /*
    Friend functions.
  */
  
  friend VVTaylorModel integrate (const VVTaylorModel&,unsigned int);
  friend VVTaylorModel derivate  (const VVTaylorModel&,unsigned int);
  friend VVTaylorModel substitute(const VVTaylorModel&,unsigned int,const Interval&);

  friend VVTaylorModel operator + (const VVTaylorModel&, const Vector&);
  friend VVTaylorModel operator - (const VVTaylorModel&, const Vector&);

  friend VVTaylorModel operator + (const Vector&, const VVTaylorModel&);
  friend VVTaylorModel operator - (const Vector&, const VVTaylorModel&); 

  friend VVTaylorModel operator * (const  Matrix&, const VVTaylorModel&);
  friend VVTaylorModel operator * (const IMatrix&, const VVTaylorModel&);

  friend std::ostream& operator <<  (std::ostream&, const VVTaylorModel&);
  friend std::string&  operator <<  (std::string& , const VVTaylorModel&);

  
 private:
  
  unsigned int row_;
  TaylorModel *_M_;
};

const VVTaylorModel operator + (const VVTaylorModel&, const VVTaylorModel&);
const VVTaylorModel operator - (const VVTaylorModel&, const VVTaylorModel&);

const VVTaylorModel operator + (const VVTaylorModel&, const TaylorModel&);
const VVTaylorModel operator - (const VVTaylorModel&, const TaylorModel&);
const VVTaylorModel operator * (const VVTaylorModel&, const TaylorModel&);
const VVTaylorModel operator / (const VVTaylorModel&, const TaylorModel&);

const VVTaylorModel operator + (const TaylorModel&, const VVTaylorModel&);
const VVTaylorModel operator - (const TaylorModel&, const VVTaylorModel&);
const VVTaylorModel operator * (const TaylorModel&, const VVTaylorModel&);

const VVTaylorModel operator + (const Interval&, const VVTaylorModel&);
const VVTaylorModel operator - (const Interval&, const VVTaylorModel&);
const VVTaylorModel operator * (const Interval&, const VVTaylorModel&);

const VVTaylorModel operator + (const VVTaylorModel&, const Interval&);
const VVTaylorModel operator - (const VVTaylorModel&, const Interval&);
const VVTaylorModel operator * (const VVTaylorModel&, const Interval&);
const VVTaylorModel operator / (const VVTaylorModel&, const Interval&);

const VVTaylorModel operator + (const double&, const VVTaylorModel&);
const VVTaylorModel operator - (const double&, const VVTaylorModel&);
const VVTaylorModel operator * (const double&, const VVTaylorModel&);

const VVTaylorModel operator + (const VVTaylorModel&, const double&);
const VVTaylorModel operator - (const VVTaylorModel&, const double&);
const VVTaylorModel operator * (const VVTaylorModel&, const double&);
const VVTaylorModel operator / (const VVTaylorModel&, const double&);
}
#endif

/*

  End of File: vvtaylormodel.h

*/


