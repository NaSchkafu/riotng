// -*-c++-*
/*

 File: adaptintval.h, 2003/05/22

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

#ifndef ADAPTINTVAL_H_INCLUDED
#define ADAPTINTVAL_H_INCLUDED

#include <sstream>

#ifdef FILIB_VERSION

#define FILIB_NAMESPACES

#include "interval/interval.hpp"
typedef filib::interval<double> _INTERVAL_TYPE_;

#endif
//#define CXSC_VERSION USE_IV_CXSC
#ifdef CXSC_VERSION

#include "interval.hpp"
#include "imath.hpp"
typedef cxsc::interval _INTERVAL_TYPE_;

#endif

/*
  Adapt the chosen interval library to a unique interface:

  Construction:
  
  Interval( double )                              :    creates point interval
  Interval( double , double )                     :    creates interval with lower and upper bound

  Access to components:

  inf(), inf( Interval )                          :    returns lower bound                 
  sup(), sup( Interval )                          :    returns upper bound                 

  Questions to Intervals:

  is_point(), is_point( Interval )                :    
  is_empty(), is_empty( Interval )                :    

  Functions:

  mid(), mid( Interval )                          :    returns mid point 
  diam(), diam( Interval )                        :    returns diameter
  abs(), abs( Interval )                          :

  Set operations:

  intersect( Interval ), intersect( Interval , Interval )   
  Interval & Interval                             :    intersection of two intervals
  hull( Interval ), hull( Interval , Interval )
  hull( Interval , double )
  hull( double , Interval )
  hull( double , double )   
  double   | Interval                               
  Interval | double                              
  double   | double                               :    union of two (point-)intervals
  in( double , Interval )                         :    checks if double is in Interval
  subset( Interval ) 
  subset( Interval , Interval )                   :
  superset( Interval ) 
  superset( Interval , Interval )                 :

  Relations:

  Interval == Interval                            :
  Interval != Interval                            :
  disjoint( Interval , Interval )                 :

  In-/Output:

  ostream << Interval                             :
  istream >> Interval                             :

  Arithmetik Operations:

  Interval o= {Interval,double}
  Interval o Interval
  Interval o double
  double   o Interval, o \in {+,-,*,/}            :

  Intrinsic Functions:

  sqr  ( Interval )
  sqrt ( Interval )
  exp  ( Interval )
  ln   ( Interval )
  sin  ( Interval )
  cos  ( Interval )
  tan  ( Interval )
  cot  ( Interval )
  asin ( Interval )
  acos ( Interval )
  atan ( Interval )
  acot ( Interval )
  sinh ( Interval )
  cosh ( Interval )
  tanh ( Interval )
  coth ( Interval )
  asinh( Interval )
  acosh( Interval )
  atanh( Interval )
  acoth( Interval )
  power( Interval , int )
  pow  ( Interval , Interval )

*/

namespace Adapt
{
  void setup();

  class Interval
    {
      _INTERVAL_TYPE_ x_;

    public:

      /*
	Ctors.
      */

      Interval(const double& d = 0.0) : x_(d,d) { setup(); }
      Interval(const double& l, const double& r) : x_(l,r) { setup(); }
      Interval(const Interval& i) : x_(i.x_) {}
      Interval(const _INTERVAL_TYPE_& x) : x_(x) { setup(); }

      Interval& operator = (const Interval& i) { x_ = i.x_; return *this; }
      
      /*
	Access.
      */

      double inf() const 
	{
#ifdef FILIB_VERSION
	  return x_.inf();
#endif
#ifdef CXSC_VERSION
	  return cxsc::_double( cxsc::Inf(x_) );
#endif
	}
      double sup() const 
	{
#ifdef FILIB_VERSION
	  return x_.sup();
#endif
#ifdef CXSC_VERSION
	  return cxsc::_double( cxsc::Sup(x_) );
#endif
	}

      /*
	State.
      */

      bool is_point() const
	{
#ifdef FILIB_VERSION
	  return x_.isPoint();
#endif
#ifdef CXSC_VERSION
	  return cxsc::Inf(x_) == cxsc::Sup(x_);
#endif
	}
      bool is_empty() const
	{
#ifdef FILIB_VERSION
	  return x_.isEmpty();
#endif
#ifdef CXSC_VERSION
	  return cxsc::IsEmpty(x_);
#endif
	}

      /*
	Useful functions.
      */

      double mid() const
	{
#ifdef FILIB_VERSION
	  return x_.mid();
#endif
#ifdef CXSC_VERSION
	  return cxsc::_double( cxsc::mid(x_) );
#endif
	}
      double diam() const
	{
#ifdef FILIB_VERSION
	  return x_.diam();
#endif
#ifdef CXSC_VERSION
	  return cxsc::_double( cxsc::diam(x_) );
#endif
	}
      Interval abs() const
	{
#ifdef FILIB_VERSION
	  return x_.abs();
#endif
#ifdef CXSC_VERSION
	  return cxsc::abs(x_);
#endif
	}

      /*
	Set operations.
      */

      Interval& intersect(const Interval& x)
	{
#ifdef FILIB_VERSION
	  x_  = x_.intersect(x.x_);
#endif
#ifdef CXSC_VERSION
	  x_ &= x.x_;
#endif
	  return *this;
	}
      Interval& operator &= (const Interval& x)
	{
#ifdef FILIB_VERSION
	  x_  = x_.intersect(x.x_);
#endif
#ifdef CXSC_VERSION
	  x_ &= x.x_;
#endif
	  return *this;
	}
      Interval& hull(const Interval& x)
	{
#ifdef FILIB_VERSION
	  x_  = x_.hull(x.x_);
#endif
#ifdef CXSC_VERSION
	  x_ |= x.x_;
#endif
	  return *this;
	}
      Interval& operator |= (const Interval& x)
	{
#ifdef FILIB_VERSION
	  x_  = x_.hull(x.x_);
#endif
#ifdef CXSC_VERSION
	  x_ |= x.x_;
#endif
	  return *this;
	}
      Interval& hull(const double& d)
	{
#ifdef FILIB_VERSION
	  x_  = x_.hull(d);
#endif
#ifdef CXSC_VERSION
	  x_ |= d;
#endif
	  return *this;
	}
      Interval& operator |= (const double& d)
	{
#ifdef FILIB_VERSION
	  x_  = x_.hull(d);
#endif
#ifdef CXSC_VERSION
	  x_ |= d;
#endif
	  return *this;
	}
      bool disjoint(const Interval& x) const
	{
#ifdef FILIB_VERSION
	  return x_.disjoint(x.x_);
#endif
#ifdef CXSC_VERSION
	  return ( (( Inf(x_) > Inf(x.x_) )? Inf(x_) : Inf(x.x_))   > (( Sup(x_) < Sup(x.x_) )? Sup(x_) : Sup(x.x_)) );
#endif
	}
      bool contains(const double& d) const
	{
#ifdef FILIB_VERSION
	  return x_.contains(d);
#endif
#ifdef CXSC_VERSION
	  return d <= x_;
#endif
	}
      bool subset(const Interval& x) const
	{
#ifdef FILIB_VERSION
	  return x_.subset(x.x_);
#endif
#ifdef CXSC_VERSION
	  return x_ <= x.x_;
#endif
	}
      bool superset(const Interval& x) const
	{
#ifdef FILIB_VERSION
	  return x_.superset(x.x_);
#endif
#ifdef CXSC_VERSION
	  return x_ >= x.x_;
#endif
	}
      Interval blow(const double& eps) const
	{
#ifdef FILIB_VERSION
	  return x_.blow(eps);
#endif
#ifdef CXSC_VERSION
	  return (1+eps) * x_ - eps * x_;
#endif
	}
      bool interior(const Interval& x) const
	{
#ifdef FILIB_VERSION
	  return x_.interior(x.x_);
#endif
#ifdef CXSC_VERSION
	  return x_ < x.x_;
#endif
	}


      /*
	Arithmetic operations.
      */

      Interval  operator -  ()                   const { return -x_; }
      Interval& operator += (const Interval& x)
	{
	  x_ += x.x_;
	  return *this;
	}
      Interval& operator -= (const Interval& x)
	{
	  x_ -= x.x_;
	  return *this;
	}
      Interval& operator *= (const Interval& x)
	{
	  x_ *= x.x_;
	  return *this;
	}
      Interval& operator /= (const Interval& x)
	{
	  x_ /= x.x_;
	  return *this;
	}
      Interval& operator += (const double& d)
	{
	  x_ += d;
	  return *this;
	}
      Interval& operator -= (const double& d)
	{
	  x_ -= d;
	  return *this;
	}
      Interval& operator *= (const double& d)
	{
	  x_ *= d;
	  return *this;
	}
      Interval& operator /= (const double& d)
	{
	  x_ /= d;
	  return *this;
	}

      /*
	Friend functions.
      */

      friend bool operator == (const Interval& x, const Interval& y) { return x.x_ == y.x_; }
      friend bool operator != (const Interval& x, const Interval& y) { return x.x_ != y.x_; }

      friend std::ostream& operator << (std::ostream& os, const Interval& x) 
	{
	  return (os << x.x_); 
	}
      friend std::istream& operator >> (std::istream& is,       Interval& x) { return (is >> x.x_); }

      friend std::string&  operator << (std::string& s,   const Interval& x) 
	{ std::ostringstream out;    out << x.x_; return (s = out.str()); }

      friend std::string&  operator >> (std::string& s,         Interval& x) 
	{ std::istringstream in(s);  in  >> x.x_; return (s = in .str()); }

      friend void     operator >> (const std::string& s,   Interval& x) 
	{ std::istringstream in(s);  in  >> x.x_;                         }

      friend Interval sqr  (const Interval& x) 
	{
#ifdef FILIB_VERSION
	  return filib::sqr(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: sqr(x.x_);
#endif
	}
      friend Interval sqrt (const Interval& x) 
	{
#ifdef FILIB_VERSION
	  return filib::sqrt(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: sqrt(x.x_);
#endif
	}

      friend Interval exp  (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::exp(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: exp(x.x_);
#endif
	}

      friend Interval ln   (const Interval& x) 
	{ 
#ifdef FILIB_VERSION
	  return filib::log(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: ln (x.x_);
#endif
	}
      friend Interval sin  (const Interval& x)                    
	{
#ifdef FILIB_VERSION
	  return filib::sin(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: sin(x.x_);
#endif
	}
      friend Interval cos  (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::cos(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: cos(x.x_);
#endif
	}
      friend Interval tan  (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::tan(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: tan(x.x_);
#endif
	}
      friend Interval cot  (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::cot(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: cot(x.x_);
#endif
	}
      friend Interval asin (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::asin(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: asin(x.x_);
#endif
	}
      friend Interval acos (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::acos(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: acos(x.x_);
#endif
	}
      friend Interval atan (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::atan(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: atan(x.x_);
#endif
	}
      friend Interval acot (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::acot(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: acot(x.x_);
#endif
	}
      friend Interval sinh (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::sinh(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: sinh(x.x_);
#endif
	}
      friend Interval cosh (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::cosh(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: cosh(x.x_);
#endif
	}
      friend Interval tanh (const Interval& x) 
	{
#ifdef FILIB_VERSION
	  return filib::tanh(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: tanh(x.x_);
#endif
	}
      friend Interval coth (const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::coth(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: coth(x.x_);
#endif
	}
      friend Interval asinh(const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::asinh(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: asinh(x.x_);
#endif
	}
      friend Interval acosh(const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::acosh(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: acosh(x.x_);
#endif
	}
      friend Interval atanh(const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::atanh(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: atanh(x.x_);
#endif
	}
      friend Interval acoth(const Interval& x)
	{
#ifdef FILIB_VERSION
	  return filib::acoth(x.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: acoth(x.x_);
#endif
	}
      friend Interval power(const Interval& x, int n)
	{
#ifdef FILIB_VERSION      
	  return filib::power(x.x_,n);  
#endif
#ifdef CXSC_VERSION
	  if     ( n <  0         ) return (Interval(1.0) /= power(x,-n));
	  else if( n == 0         ) return Interval(1.0);
	  else if( n == 1         ) return x;
	  else if( n == 2         ) return cxsc::sqr(x.x_);
	  else                      return cxsc::power( x.x_, n );
#endif
	}
      friend Interval pow  (const Interval& x, const Interval& y) 
	{
#ifdef FILIB_VERSION
	  return filib::pow(x.x_,y.x_);
#endif
#ifdef CXSC_VERSION
	  return cxsc:: pow(x.x_,y.x_);
#endif
	}

    };

  inline double inf(const Interval& x) { return x.inf(); }
  inline double sup(const Interval& x) { return x.sup(); }

  inline bool is_point(const Interval& x) { return x.is_point(); }
  inline bool is_empty(const Interval& x) { return x.is_empty(); }

  inline double   mid (const Interval& x) { return x.mid();  }
  inline double   diam(const Interval& x) { return x.diam(); }
  inline Interval abs (const Interval& x) { return x.abs();  }

  inline bool disjoint(const Interval& x, const Interval& y) { return x.disjoint(y); }

  inline Interval intersect  (const Interval& x, const Interval& y) 
    { 
      Interval z(x);
      return z.intersect(y); 
    }
  inline Interval operator & (const Interval& x, const Interval& y) 
    { 
      Interval z(x);
      return z.intersect(y); 
    }
  inline Interval hull       (const Interval& x, const Interval& y) 
    {
      Interval z(x);
      return z.hull(y);      
    }
  inline Interval operator | (const Interval& x, const Interval& y) 
    {
      Interval z(x);
      return z.hull(y);      
    }
  inline Interval hull       (const double&   d, const Interval& y) 
    {
      Interval z(d);
      return z.hull(y);      
    }
  inline Interval hull       (const Interval& x, const double&   d) 
    {
      Interval z(x);
      return z.hull(d); 
     }
  inline Interval hull       (const double&   d, const double&   e) 
    {
      Interval z(d);
      return z.hull(e);      
    }
  inline Interval operator | (const double&   d, const Interval& y) 
    {
      Interval z(d);
      return z.hull(y);      
    }
  inline Interval operator | (const Interval& x, const double&   d) 
    {
      Interval z(x);
      return z.hull(d);    
    }

  inline bool in           (const double&   d, const Interval& x) { return x.contains(d); }
  inline bool subset       (const Interval& x, const Interval& y) { return x.subset(y);   }
  inline bool operator <=  (const Interval& x, const Interval& y) { return x.subset(y);   }
  inline bool superset     (const Interval& x, const Interval& y) { return x.superset(y); }
  inline bool operator >=  (const Interval& x, const Interval& y) { return x.superset(y); }

  inline Interval blow(const Interval& x, const double& eps) { return x.blow(eps); }

  inline bool interior(const Interval& x, const Interval& y) { return x.interior(y); }

  inline Interval operator + (const Interval& x, const Interval& y) 
    {
      Interval z(x);
      return z += y; 
    }
  inline Interval operator - (const Interval& x, const Interval& y) 
    {
      Interval z(x);
      return z -= y; 
    }
  inline Interval operator * (const Interval& x, const Interval& y) 
    {
      Interval z(x);
      return z *= y; 
    }
  inline Interval operator / (const Interval& x, const Interval& y) 
    {
      Interval z(x);
      return z /= y; 
    }

  inline Interval subtract_bounds(const Interval& x, const Interval& y)
    {
#ifdef FILIB_VERSION
      double 
	lb = filib::fp_traits<double>::downward_minus( x.inf(), y.inf(), true ),
	ub = filib::fp_traits<double>::upward_minus( x.sup(), y.sup(), false );
#endif
#ifdef CXSC_VERSION
      double 
	lb = cxsc::_double( cxsc::subdown( x.inf(), y.inf() ) ),
	ub = cxsc::_double( cxsc::subup  ( x.sup(), y.sup() ) );
#endif
      return Interval(lb,ub);
    }

  inline Interval operator + (const double& d, const Interval& y) 
    {
      Interval z(d);
      return z += y; 
    }
  inline Interval operator - (const double& d, const Interval& y) 
    {
      Interval z(d);
      return z -= y; 
    }
  inline Interval operator * (const double& d, const Interval& y) 
    {
      Interval z(d);
      return z *= y; 
    }
  inline Interval operator / (const double& d, const Interval& y) 
    {
      Interval z(d);
      return z /= y; 
    }

  inline Interval operator + (const Interval& x, const double& d) 
    {
      Interval z(x);
      return z += d; 
    }
  inline Interval operator - (const Interval& x, const double& d) 
    {
      Interval z(x);
      return z -= d; 
    }
  inline Interval operator * (const Interval& x, const double& d) 
    { 
      Interval z(x);
      return z *= d; 
    }
  inline Interval operator / (const Interval& x, const double& d) 
    {
      Interval z(x);
      return z /= d; 
    }

  /*
    Useful constants.
  */
  
  inline const Interval& ZERO_INTERVAL()
    { 
      static const Interval zero = Interval(0.0);
      return zero;
    }

  inline const Interval& ONE_INTERVAL()
    { 
      static const Interval one = Interval(1.0);
      return one;
    }
  
  inline const Interval& PI()
    { 
      static const Interval pi = acos(Interval(-1.0));
      return pi;
    }

  inline const Interval& HALFPI()
    { 
      static const Interval hp = PI() / 2.0;
      return hp;
    }

  inline const Interval& QUARTERPI()
    {
      static const Interval qp = PI() / 4.0;
      return qp;
    }

}

#endif

/*

End of File: adaptintval.h

*/
