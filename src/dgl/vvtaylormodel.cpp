/*

 File: vvtaylormodel.cpp, 2003/10/08

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

#include "vvtaylormodel.h"

/*

The read only functions.

*/
namespace riot {
IVector VVTaylorModel::interval_part() const
{
  IVector res(row_);

  for(unsigned int i = 0; i < row_; i++) res[i] = _M_[i].interval_part();

  return res;
}
Vector VVTaylorModel::constant_part() const
{
  Vector res(row_);

  for(unsigned int i = 0; i < row_; i++) res[i] = _M_[i].coefficient_of( Monom() );
  
  return res;
}
IVector VVTaylorModel::eval() const
{
  IVector res(row_);

  for(unsigned int i = 0; i < row_; i++) res[i] = _M_[i].eval();

  return res;
}

/*

The operators.

*/

VVTaylorModel VVTaylorModel::operator - () const
{
  VVTaylorModel res;
  res.row_ = row_;
  res._M_  = new TaylorModel[row_];

  for(unsigned int i = 0; i < row_; i++) res._M_[i] = - _M_[i];

  return res;
}

VVTaylorModel& VVTaylorModel::operator += (const VVTaylorModel& A)
{
  if( row_ == A.row_ )
    for(unsigned int i = 0; i < row_; i++) _M_[i] += A._M_[i];
  else 
    std::cout << "VVTaylorModel += VVTaylorModel: *** Wrong dimensions! ***" << std::endl;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator -= (const VVTaylorModel& A)
{
  if( row_ == A.row_ )
    for(unsigned int i = 0; i < row_; i++) _M_[i] -= A._M_[i];
  else 
    std::cout << "VVTaylorModel -= VVTaylorModel: *** Wrong dimensions! ***" << std::endl;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator += (const TaylorModel& s)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] += s;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator -= (const TaylorModel& s)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] -= s;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator *= (const TaylorModel& s)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] *= s;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator /= (const TaylorModel& s)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] /= s;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator += (const Interval& a)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] += a;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator -= (const Interval& a)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] -= a;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator *= (const Interval& a)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] *= a;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator /= (const Interval& a)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] /= a;

  return *this;
}

VVTaylorModel& VVTaylorModel::operator += (const double& d)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] += d;

  return *this;
}
VVTaylorModel& VVTaylorModel::operator -= (const double& d)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] -= d;

  return *this;
}
VVTaylorModel& VVTaylorModel::operator *= (const double& d)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] *= d;

  return *this;
}
VVTaylorModel& VVTaylorModel::operator /= (const double& d)
{
  for(unsigned int i = 0; i < row_; i++) _M_[i] /= d;

  return *this;
}

void VVTaylorModel::add_to_interval_part(const IVector& A)
{
  if( row_ == A.size() )
    for(unsigned int i = 0; i < row_; i++) _M_[i].add_to_interval_part( A[i] );
  else
    std::cout << "VVTaylorModel::add_to_interval_part: *** Wrong dimensions ***" << std::endl;
}

void VVTaylorModel::replace_interval_part_with(const IVector& A)
{
  if( row_ == A.size() )
    for(unsigned int i = 0; i < row_; i++) _M_[i].replace_interval_part_with( A[i] );
  else
    std::cout << "VVTaylorModel::replace_interval_part_with: *** Wrong dimensions ***" << std::endl;
}

//void subtract_polynomial_part_of(const VVTaylorModel&);
//void remove_from_polynomial(const Monom&);

/*
  
Friend functions.

*/
  
VVTaylorModel integrate (const VVTaylorModel& A, unsigned int code)
{
  VVTaylorModel res;
  res.row_ = A.row_;
  res._M_  = new TaylorModel[A.row_];

  for(unsigned int i = 0; i < A.row_; i++) res._M_[i] = integrate( A._M_[i], code );

  return res;
}

VVTaylorModel derivate  (const VVTaylorModel& A, unsigned int code)
{
  VVTaylorModel res;
  res.row_ = A.row_;
  res._M_  = new TaylorModel[A.row_];

  for(unsigned int i = 0; i < A.row_; i++) res._M_[i] = derivate( A._M_[i], code );

  return res;
}

VVTaylorModel substitute(const VVTaylorModel& A, unsigned int code, const Interval& domain)
{
  VVTaylorModel res;
  res.row_ = A.row_;
  res._M_  = new TaylorModel[A.row_];

  for(unsigned int i = 0; i < A.row_; i++) res._M_[i] = substitute( A._M_[i], code, domain );

  return res;
}

VVTaylorModel operator + (const VVTaylorModel& A, const Vector& v)
{
  VVTaylorModel res(A);

  if( A.row_ == v.size() )
    for(unsigned int i = 0; i < res.row_; i++) res._M_[i] += v[i];
  else
    std::cout << "VVTaylorModel + Vector: *** Wrong dimensions ***" << std::endl;

  return res;
}

VVTaylorModel operator - (const VVTaylorModel& A, const Vector& v)
{
  VVTaylorModel res(A);

  if( A.row_ == v.size() )
    for(unsigned int i = 0; i < res.row_; i++) res._M_[i] -= v[i];
  else
    std::cout << "VVTaylorModel - Vector: *** Wrong dimensions ***" << std::endl;

  return res;
}

VVTaylorModel operator + (const Vector& v, const VVTaylorModel& A)
{
  VVTaylorModel res(A);

  if( A.row_ == v.size() )
    for(unsigned int i = 0; i < res.row_; i++) res._M_[i] += v[i];
  else
    std::cout << "Vector + VVTaylorModel: *** Wrong dimensions ***" << std::endl;

  return res;
}

VVTaylorModel operator - (const Vector& v, const VVTaylorModel& A)
{
  VVTaylorModel res(-A);

  if( A.row_ == v.size() )
    for(unsigned int i = 0; i < res.row_; i++) res._M_[i] += v[i];
  else
    std::cout << "Vector - VVTaylorModel: *** Wrong dimensions ***" << std::endl;

  return res;
}

VVTaylorModel operator * (const  Matrix& M, const VVTaylorModel& A)
{
  VVTaylorModel res;
  res.row_ = M.size(1); //Number of rows.
  res._M_  = new TaylorModel[res.row_];

  if( A.row_ == M.size(2) )
    for(unsigned int i = 0; i < res.row_; i++)
      {
	res._M_[i] = TaylorModel(TaylorModel::ZERO_TM());
	for(unsigned int j = 0; j < A.row_; j++) res._M_[i] += M(i,j) * A._M_[j];
      }
  else
    std::cout << "Matrix * VVTaylorModel: *** Wrong dimensions ***" << std::endl;

  return res;
}

VVTaylorModel operator * (const IMatrix& M, const VVTaylorModel& A)
{
  VVTaylorModel res;
  res.row_ = M.size(1); //Number of rows.
  res._M_  = new TaylorModel[res.row_];

  if( A.row_ == M.size(2) )
    for(unsigned int i = 0; i < res.row_; i++)
      {
	res._M_[i] = TaylorModel::ZERO_TM();
	for(unsigned int j = 0; j < A.row_; j++) res._M_[i] += M(i,j) * A._M_[j];
      }
  else
    std::cout << "IMatrix * VVTaylorModel: *** Wrong dimensions ***" << std::endl;

  return res;
}

std::ostream& operator <<  (std::ostream& os, const VVTaylorModel& A)
{
  for(unsigned int i = 0; i < A.row_; i++)
    os << i << ". Komponente:\n\n" << A._M_[i] << std::endl;

  return os;
}

std::string&  operator <<  (std::string& s, const VVTaylorModel& A)
{
  std::cout << "string << VVTaylorModel: *** Noch nicht implementiert ***" << std::endl;
  return s;
}

/*

Sonstige Operatoren.

*/

const VVTaylorModel operator + (const VVTaylorModel& A, const VVTaylorModel& B)
{
  VVTaylorModel C(A);
  return C += B;
}

const VVTaylorModel operator - (const VVTaylorModel& A, const VVTaylorModel& B)
{
  VVTaylorModel C(A);
  return C -= B;
}

const VVTaylorModel operator + (const VVTaylorModel& A, const TaylorModel& B)
{
  VVTaylorModel C(A);
  return C += B;
}

const VVTaylorModel operator - (const VVTaylorModel& A, const TaylorModel& B)
{
  VVTaylorModel C(A);
  return C -= B;
}

const VVTaylorModel operator * (const VVTaylorModel& A, const TaylorModel& B)
{
  VVTaylorModel C(A);
  return C *= B;
}

const VVTaylorModel operator / (const VVTaylorModel& A, const TaylorModel& B)
{
  VVTaylorModel C(A);
  return C /= B;
}

const VVTaylorModel operator + (const TaylorModel& A, const VVTaylorModel& B)
{
  VVTaylorModel C(B);
  return C += A;
}

const VVTaylorModel operator - (const TaylorModel& A, const VVTaylorModel& B)
{
  VVTaylorModel C(-B);
  return C += A;
}

const VVTaylorModel operator * (const TaylorModel& A, const VVTaylorModel& B)
{
  VVTaylorModel C(B);
  return C *= A;
}

const VVTaylorModel operator + (const Interval& I, const VVTaylorModel& A)
{
  VVTaylorModel C(A);
  return C += I;
}

const VVTaylorModel operator - (const Interval& I, const VVTaylorModel& A)
{
  VVTaylorModel C(-A);
  return C += I;
}

const VVTaylorModel operator * (const Interval& I, const VVTaylorModel& A)
{
  VVTaylorModel C(A);
  return C *= I;
}

const VVTaylorModel operator + (const VVTaylorModel& A, const Interval& I)
{
  VVTaylorModel C(A);
  return C += I;
}

const VVTaylorModel operator - (const VVTaylorModel& A, const Interval& I)
{
  VVTaylorModel C(A);
  return C -= I;
}

const VVTaylorModel operator * (const VVTaylorModel& A, const Interval& I)
{
  VVTaylorModel C(A);
  return C *= I;
}

const VVTaylorModel operator / (const VVTaylorModel& A, const Interval& I)
{
  VVTaylorModel C(A);
  return C /= I;
}

const VVTaylorModel operator + (const double& d, const VVTaylorModel& A)
{
  VVTaylorModel C(A);
  return C += d;
}

const VVTaylorModel operator - (const double& d, const VVTaylorModel& A)
{
  VVTaylorModel C(-A);
  return C += d;
}

const VVTaylorModel operator * (const double& d, const VVTaylorModel& A)
{
  VVTaylorModel C(A);
  return C *= d;
}

const VVTaylorModel operator + (const VVTaylorModel& A, const double& d)
{
  VVTaylorModel C(A);
  return C += d;
}

const VVTaylorModel operator - (const VVTaylorModel& A, const double& d)
{
  VVTaylorModel C(A);
  return C -= d;
}

const VVTaylorModel operator * (const VVTaylorModel& A, const double& d)
{
  VVTaylorModel C(A);
  return C *= d;
}

const VVTaylorModel operator / (const VVTaylorModel& A, const double& d)
{
  VVTaylorModel C(A);
  return C /= d;
}

/*

  End of File: vvtaylormodel.h

*/
}
