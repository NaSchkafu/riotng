/*

  File: taylormodel.cpp, 2004/08/18

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

#include "taylormodel.h"

/*
  Definition of friend functions.
*/
namespace riot
{

  TaylorModel sqr (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( sqr (*(s._M_RC_)) );
    return result;
  }

  TaylorModel sqrt(const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( sqrt(*(s._M_RC_)) );
    return result;
  }

  TaylorModel invsqrt(const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( invsqrt(*(s._M_RC_)) );
    return result;
  }

  TaylorModel exp (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( exp (*(s._M_RC_)) );
    return result;
  }

  TaylorModel log (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( log (*(s._M_RC_)) );
    return result;
  }

  TaylorModel sin (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( sin (*(s._M_RC_)) );
    return result;
  }

  TaylorModel cos (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( cos (*(s._M_RC_)) );
    return result;
  }

  TaylorModel tan  (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( tan (*(s._M_RC_)) );
    return result;
  }

  TaylorModel cot  (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( cot (*(s._M_RC_)) );
    return result;
  }

  TaylorModel asin (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( asin (*(s._M_RC_)) );
    return result;
  }

  TaylorModel acos (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( acos (*(s._M_RC_)) );
    return result;
  }

  TaylorModel atan (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( atan (*(s._M_RC_)) );
    return result;
  }

  TaylorModel acot (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( acot (*(s._M_RC_)) );
    return result;
  }

  TaylorModel sinh (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( sinh (*(s._M_RC_)) );
    return result;
  }

  TaylorModel cosh (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( cosh (*(s._M_RC_)) );
    return result;
  }

  TaylorModel tanh (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( tanh (*(s._M_RC_)) );
    return result;
  }

  TaylorModel coth (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( coth (*(s._M_RC_)) );
    return result;
  }

  TaylorModel asinh(const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( asinh (*(s._M_RC_)) );
    return result;
  }

  TaylorModel acosh(const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( acosh (*(s._M_RC_)) );
    return result;
  }

  TaylorModel atanh(const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( atanh (*(s._M_RC_)) );
    return result;
  }

  TaylorModel acoth(const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( acoth (*(s._M_RC_)) );
    return result;
  }

  TaylorModel pow  (const TaylorModel& s,const TaylorModel& t)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( pow (*(s._M_RC_),*(t._M_RC_)) );
    return result;
  }

  TaylorModel power(const TaylorModel& s,int n)
  {
    if     ( n == 0 ) return TaylorModel::ONE_TM();
    else if( n == 1 )
    {
      TaylorModel result(s);
      return result;
    }
    else if( n == 2 ) return sqr(s);
    else
    {
      TaylorModel result;
      result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( power (*(s._M_RC_),n) );
      return result;
    }
  }

  TaylorModel invert    (const TaylorModel& s)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( invert(*(s._M_RC_)) );
    return result;
  }

  TaylorModel integrate (const TaylorModel& s, unsigned int var_code)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( integrate (*(s._M_RC_),var_code) );
    return result;
  }

  TaylorModel derivate  (const TaylorModel& s, unsigned int var_code)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( derivate  (*(s._M_RC_),var_code) );
    return result;
  }

  TaylorModel substitute(const TaylorModel& s,unsigned int code,const Interval& iv)
  {
    TaylorModel result;
    result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( substitute (*(s._M_RC_),code,iv) );
    return result;
  }

  std::ostream& operator <<  (std::ostream& os, const TaylorModel& s)
  {
    return os << *(s._M_RC_);
  }

  std::string&  operator <<  (std::string& st, const TaylorModel& s)
  {
    return st << *(s._M_RC_);
  }


/*
  Binary operators: TaylorModel o TaylorModel, o \in {+,-,*,/}.
*/

  const TaylorModel operator + (const TaylorModel& s, const TaylorModel& t)
  {
    TaylorModel result(s);
    return (result += t);
  }

  const TaylorModel operator - (const TaylorModel& s, const TaylorModel& t)
  {
    TaylorModel result(s);
    return (result -= t);
  }

  const TaylorModel operator * (const TaylorModel& s, const TaylorModel& t)
  {
    TaylorModel result(s);
    return (result *= t);
  }

  const TaylorModel operator / (const TaylorModel& s, const TaylorModel& t)
  {
    TaylorModel result(s);
    return (result /= t);
  }

/*
  Binary operators: Interval o TaylorModel, o \in {+,-,*,/}.
*/

  const TaylorModel operator + (const Interval& i, const TaylorModel& t)
  {
    TaylorModel result(t);
    return (result += i);
  }

  const TaylorModel operator - (const Interval& i, const TaylorModel& t)
  {
    TaylorModel result(-t);
    return (result += i);
  }

  const TaylorModel operator * (const Interval& i, const TaylorModel& t)
  {
    TaylorModel result(t);
    return (result *= i);
  }

  const TaylorModel operator / (const Interval& i, const TaylorModel& t)
  {
    return (i * invert(t));
  }

/*
  Binary operators: TaylorModel o Interval, o \in {+,-,*,/}.
*/

  const TaylorModel operator + (const TaylorModel& t, const Interval& i)
  {
    TaylorModel result(t);
    return (result += i);
  }

  const TaylorModel operator - (const TaylorModel& t, const Interval& i)
  {
    TaylorModel result(t);
    return (result -= i);
  }

  const TaylorModel operator * (const TaylorModel& t, const Interval& i)
  {
    TaylorModel result(t);
    return (result *= i);
  }

  const TaylorModel operator / (const TaylorModel& t, const Interval& i)
  {
    TaylorModel result(t);
    return (result /= i);
  }

/*
  Binary operators: double o TaylorModel, o \in {+,-,*,/}.
*/

  const TaylorModel operator + (const double& d, const TaylorModel& t)
  {
    TaylorModel result(t);
    return (result += d);
  }

  const TaylorModel operator - (const double& d, const TaylorModel& t)
  {
    TaylorModel result(-t);
    return (result += d);
  }

  const TaylorModel operator * (const double& d, const TaylorModel& t)
  {
    TaylorModel result(t);
    return (result *= d);
  }

  const TaylorModel operator / (const double& d, const TaylorModel& t)
  {
    return (d * invert(t));
  }

/*
  Binary operators: TaylorModel o double, o \in {+,-,*,/}.
*/

  const TaylorModel operator + (const TaylorModel& t, const double& d)
  {
    TaylorModel result(t);
    return (result += d);
  }

  const TaylorModel operator - (const TaylorModel& t, const double& d)
  {
    TaylorModel result(t);
    return (result -= d);
  }

  const TaylorModel operator * (const TaylorModel& t, const double& d)
  {
    TaylorModel result(t);
    return (result *= d);
  }

  const TaylorModel operator / (const TaylorModel& t, const double& d)
  {
    TaylorModel result(t);
    return (result /= d);
  }


/*

  End of File: taylormodel.cpp

*/
}
