// -*-c++-*-
/*

  File: taylormodel.h, 2004/08/28

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

#ifndef TAYLORMODEL_H_INCLUDED
#define TAYLORMODEL_H_INCLUDED

#include "refcounter.h"
#include "taylormodel_impl.h"

#include "monom.h"

using Adapt::Interval;

namespace riot
{
  class TaylorModel
  {
  public:

    /*
      Constructors
    */

    TaylorModel() {}
    TaylorModel(const std::string& vname, const double& dpoint, const Interval& domain)
      : _M_RC_( new TaylorModel_Impl(vname,dpoint,domain) )
      {}

    /*
      Read access
    */

    const Interval&    interval_part()                const { return _M_RC_->interval_part();      }
    double             coefficient_of(const Monom& m) const { return _M_RC_->coefficient_of(m);    }
    std::vector<Monom> monom_of_order(unsigned int n) const { return _M_RC_->monom_of_order(n);    }
    Interval           eval()                         const { return _M_RC_->eval();               }
    bool               is_zero()                      const { return _M_RC_->is_zero();            }
    bool               polynomial_is_zero()           const { return _M_RC_->polynomial_is_zero(); }

    /*
      Write access
    */

    TaylorModel operator - () const
      {
        TaylorModel t;
        t._M_RC_ = ReferenceCounter<TaylorModel_Impl>(-(*_M_RC_));
        return t;
      }
    TaylorModel& operator += (const TaylorModel& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) += *(s._M_RC_);

        return *this;
      }
    TaylorModel& operator -= (const TaylorModel& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) -= *(s._M_RC_);

        return *this;
      }
    TaylorModel& operator *= (const TaylorModel& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) *= *(s._M_RC_);

        return *this;
      }
    TaylorModel& operator /= (const TaylorModel& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) /= *(s._M_RC_);

        return *this;
      }

    TaylorModel& subtract_polynomial_part_of(const TaylorModel& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_).subtract_polynomial_part_of( *(s._M_RC_) );

        return *this;
      }

    TaylorModel& operator += (const Interval& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) += s;

        return *this;
      }
    TaylorModel& operator -= (const Interval& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) -= s;

        return *this;
      }
    TaylorModel& operator *= (const Interval& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) *= s;

        return *this;
      }
    TaylorModel& operator /= (const Interval& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) /= s;

        return *this;
      }

    TaylorModel& operator += (const double& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) += s;

        return *this;
      }
    TaylorModel& operator -= (const double& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) -= s;

        return *this;
      }
    TaylorModel& operator *= (const double& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) *= s;

        return *this;
      }
    TaylorModel& operator /= (const double& s)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_) /= s;

        return *this;
      }

    void add_to_interval_part(const Interval& I)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_).add_to_interval_part(I);
      }
    void replace_interval_part_with(const Interval& I)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_).replace_interval_part_with(I);
      }
    void remove_from_polynomial(const Monom& m)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        (*_M_RC_).remove_from_polynomial(m);
      }
    double read_coefficient_and_remove(const Monom& m)
      {
        if( _M_RC_ > 1 )
          _M_RC_ = ReferenceCounter<TaylorModel_Impl>(new TaylorModel_Impl(*_M_RC_));

        return (*_M_RC_).read_coefficient_and_remove(m);
      }
    void change_var_data(unsigned int code, const double& dp, const Interval& dm)
      {
        (*_M_RC_).change_var_data(code,dp,dm);
      }

    /*
      Friend functions
    */

    friend TaylorModel sqr    (const TaylorModel&);
    friend TaylorModel sqrt   (const TaylorModel&);
    friend TaylorModel invsqrt(const TaylorModel&);
    friend TaylorModel exp    (const TaylorModel&);
    friend TaylorModel log    (const TaylorModel&);
    friend TaylorModel sin    (const TaylorModel&);
    friend TaylorModel cos    (const TaylorModel&);
    friend TaylorModel tan    (const TaylorModel&);
    friend TaylorModel cot    (const TaylorModel&);
    friend TaylorModel asin   (const TaylorModel&);
    friend TaylorModel acos   (const TaylorModel&);
    friend TaylorModel atan   (const TaylorModel&);
    friend TaylorModel acot   (const TaylorModel&);
    friend TaylorModel sinh   (const TaylorModel&);
    friend TaylorModel cosh   (const TaylorModel&);
    friend TaylorModel tanh   (const TaylorModel&);
    friend TaylorModel coth   (const TaylorModel&);
    friend TaylorModel asinh  (const TaylorModel&);
    friend TaylorModel acosh  (const TaylorModel&);
    friend TaylorModel atanh  (const TaylorModel&);
    friend TaylorModel acoth  (const TaylorModel&);
    friend TaylorModel pow    (const TaylorModel&,const TaylorModel&);
    friend TaylorModel power  (const TaylorModel&,int);

    friend TaylorModel invert    (const TaylorModel&);
    friend TaylorModel integrate (const TaylorModel&,unsigned int);
    friend TaylorModel derivate  (const TaylorModel&,unsigned int);
    friend TaylorModel substitute(const TaylorModel&,unsigned int,const Interval&);

    friend std::ostream& operator <<  (std::ostream&, const TaylorModel&);
    friend std::string&  operator <<  (std::string& , const TaylorModel&);

    /*
      Static constants
    */

    static TaylorModel& ZERO_TM()
      {
        static TaylorModel zero( Const_TM(0.0) );
        return zero;
      }

    static TaylorModel& ONE_TM()
      {
        static TaylorModel one( Const_TM(1.0) );
        return one;
      }

    /*
      Functions to read and change the static members
    */

    static unsigned int order                 ()
      {
        return TaylorModel_Impl::order();
      }
    static unsigned int set_order             (unsigned int n)
      {
        return TaylorModel_Impl::set_order(n);
      }
    static double       set_sparsity_tol      (const double& tol)
      {
        return TaylorModel_Impl::set_sparsity_tol(tol);
      }
    static void         set_polynomial_range_evaluation(const PolynomialRange& pb)
      {
        TaylorModel_Impl::set_polynomial_range_evaluation(pb);
      }
    static void         set_degree_check(const DegreeBase& db)
      {
        TaylorModel_Impl::set_degree_check(db);
      }
    static void         set_mode              (unsigned int m)
      {
        TaylorModel_Impl::set_mode(m);
      }
    static void         mode_message          (unsigned int m)
      {
        TaylorModel_Impl::mode_message(m);
      }

  private:

    static TaylorModel Const_TM(const double& d)
      {
        TaylorModel result;
        result._M_RC_ = ReferenceCounter<TaylorModel_Impl>( TaylorModel_Impl::Const_TM(d) );

        return result;
      }

    /*
      Data
    */

    ReferenceCounter<TaylorModel_Impl> _M_RC_;
  };

  const TaylorModel operator + (const TaylorModel&, const TaylorModel&);
  const TaylorModel operator - (const TaylorModel&, const TaylorModel&);
  const TaylorModel operator * (const TaylorModel&, const TaylorModel&);
  const TaylorModel operator / (const TaylorModel&, const TaylorModel&);

  const TaylorModel operator + (const Interval&, const TaylorModel&);
  const TaylorModel operator - (const Interval&, const TaylorModel&);
  const TaylorModel operator * (const Interval&, const TaylorModel&);
  const TaylorModel operator / (const Interval&, const TaylorModel&);

  const TaylorModel operator + (const TaylorModel&, const Interval&);
  const TaylorModel operator - (const TaylorModel&, const Interval&);
  const TaylorModel operator * (const TaylorModel&, const Interval&);
  const TaylorModel operator / (const TaylorModel&, const Interval&);

  const TaylorModel operator + (const double&, const TaylorModel&);
  const TaylorModel operator - (const double&, const TaylorModel&);
  const TaylorModel operator * (const double&, const TaylorModel&);
  const TaylorModel operator / (const double&, const TaylorModel&);

  const TaylorModel operator + (const TaylorModel&, const double&);
  const TaylorModel operator - (const TaylorModel&, const double&);
  const TaylorModel operator * (const TaylorModel&, const double&);
  const TaylorModel operator / (const TaylorModel&, const double&);

}


#endif

/*

  End of File: taylormodel.h

*/
