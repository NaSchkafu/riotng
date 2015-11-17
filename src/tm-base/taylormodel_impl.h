// -*-c++-*-
/*

  File: taylormodel_impl.h, 2004/12/30

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

#ifndef TAYLORMODEL_IMPL_H_INCLUDED
#define TAYLORMODEL_IMPL_H_INCLUDED

#include <cstdlib>

#include "adaptintval.h"
#include "vardatatable.h"
#include "refcounter.h"
#include "hashtable.h"
#include "hashfunc.h"
#include "monom.h"
#include "degreeeval.h"

namespace riot
{

  class PolynomialRange;
  using Adapt::Interval;

/*-----------------------------------------------------------------------------+

  It follows the implementation of type Taylor Model based on hashtables.

  Notice: The zero polynomial has no entries in the hashtable!

*/

  enum {
    NORMAL_MODE    = 0,
    REMAINDER_REC  = 1,
    REMAINDER_PLAY = 2,
    NO_REMAINDER   = 4
  };

  enum {
    SEQUENCE_BREAK = 0
  };

  class TaylorModel_Impl
  {
  public:

    /*
      Introduce better readable names
    */

    typedef HashTable<Monom,double,MonomHasher2>                   polynomial_type;
    typedef HashTable<Monom,double,MonomHasher2>::iterator         polynomial_iterator;
    typedef HashTable<Monom,double,MonomHasher2>::const_iterator   const_polynomial_iterator;

    typedef HashTable<Monom,Interval,MonomHasher2>                 ipolynomial_type;
    typedef HashTable<Monom,Interval,MonomHasher2>::iterator       ipolynomial_iterator;
    typedef HashTable<Monom,Interval,MonomHasher2>::const_iterator const_ipolynomial_iterator;

    typedef std::vector<ReferenceCounter<VarDataListNode> >        additional_var_data;

    /*

      Constructors

    */

    TaylorModel_Impl() : polynomial_(1,MonomHasher2::getFunction()), rest_(0.0), data_(0) {}
    TaylorModel_Impl(const std::string& vname, const double& dpoint, const Interval& domain)
      : polynomial_(3,MonomHasher2::getFunction()), rest_(0.0)
      {
        unsigned int var_code = vtab_->insert( vname );

        polynomial_.insert( Monom( var_code, 1 ), //Monom to insert.
                            1.0,                  //The coefficient.
                            DefaultAction<Monom,double>::getObject() );

        data_ = additional_var_data( var_code );

        VarDataListNode *data_ptr = VarData::getTable().insert( var_code, dpoint, domain );
        data_[ var_code - 1 ] = ReferenceCounter<VarDataListNode>( data_ptr, data_ptr->cnt );
      }
    TaylorModel_Impl(const TaylorModel_Impl& s)
      : polynomial_(s.polynomial_), rest_( s.rest_ ), data_( s.data_ )
      {}

    ~TaylorModel_Impl() {}

    TaylorModel_Impl& operator = (const TaylorModel_Impl& s)
      {
        polynomial_ = s.polynomial_;
        rest_       = s.rest_;
        data_       = s.data_;

        return *this;
      }

    /*
      Read access
    */

    double             coefficient_of(const Monom&) const;
    std::vector<Monom> monom_of_order(unsigned int) const;
    Interval           eval()               const;
    const Interval&    interval_part()      const { return rest_;                                               }
    bool               is_zero()            const { return polynomial_.number_of_entries() == 0 && rest_ == 0; }
    bool               polynomial_is_zero() const { return polynomial_.number_of_entries() == 0;                }

    /*
      Write access
    */

    TaylorModel_Impl* operator - () const;

    TaylorModel_Impl& operator += (const TaylorModel_Impl& s)
      {
        if( Mode & NO_REMAINDER ) return add_without_seizing_rounding_errors( s );
        else                      return add_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator -= (const TaylorModel_Impl& s)
      {
        if( Mode & NO_REMAINDER ) return subtract_without_seizing_rounding_errors( s );
        else                      return subtract_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator *= (const TaylorModel_Impl& s)
      {
        if( Mode & NO_REMAINDER ) return multiply_without_seizing_rounding_errors( s );
        else                      return multiply_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator /= (const TaylorModel_Impl&);

    TaylorModel_Impl& subtract_polynomial_part_of(const TaylorModel_Impl&);

    TaylorModel_Impl& operator += (const Interval& s)
      {
        if( Mode & NO_REMAINDER ) return add_without_seizing_rounding_errors( s );
        else                      return add_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator -= (const Interval& s)
      {
        if( Mode & NO_REMAINDER ) return subtract_without_seizing_rounding_errors( s );
        else                      return subtract_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator *= (const Interval& s)
      {
        if( Mode & NO_REMAINDER ) return multiply_without_seizing_rounding_errors( s );
        else                      return multiply_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator /= (const Interval& s)
      {
        if( Mode & NO_REMAINDER ) return divide_without_seizing_rounding_errors( s );
        else                      return divide_with_seizing_rounding_errors   ( s );
      }

    TaylorModel_Impl& operator += (const double& s)
      {
        if( Mode & NO_REMAINDER ) return add_without_seizing_rounding_errors( s );
        else                      return add_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator -= (const double& s)
      {
        if( Mode & NO_REMAINDER ) return subtract_without_seizing_rounding_errors( s );
        else                      return subtract_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator *= (const double& s)
      {
        if( Mode & NO_REMAINDER ) return multiply_without_seizing_rounding_errors( s );
        else                      return multiply_with_seizing_rounding_errors   ( s );
      }
    TaylorModel_Impl& operator /= (const double& s)
      {
        if( Mode & NO_REMAINDER ) return divide_without_seizing_rounding_errors( s );
        else                      return divide_with_seizing_rounding_errors   ( s );
      }

    void   add_to_interval_part(const Interval& I)       { rest_ += I; }
    void   replace_interval_part_with(const Interval& I) { rest_  = I; }
    void   remove_from_polynomial(const Monom& m)        { polynomial_.erase(m); }
    double read_coefficient_and_remove(const Monom&);
    void   change_var_data(unsigned int code, const double& dp, const Interval& dm)
      {
        if( code <= data_.size() )
        {
          data_[ code - 1 ]->devel_point_ = dp;
          data_[ code - 1 ]->domain_      = dm;
        }
      }

    /*
      Iterating the polynomial
    */

    polynomial_iterator       begin()       { return polynomial_.begin(); }
    const_polynomial_iterator begin() const { return polynomial_.begin(); }
    polynomial_iterator       end()         { return polynomial_.end();   }
    const_polynomial_iterator end()   const { return polynomial_.end();   }

    /*
      Friend functions
    */

    friend TaylorModel_Impl* sqr    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* sqrt   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* invsqrt(const TaylorModel_Impl&);
    friend TaylorModel_Impl* exp    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* log    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* sin    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* cos    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* tan    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* cot    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* asin   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* acos   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* atan   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* acot   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* sinh   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* cosh   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* tanh   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* coth   (const TaylorModel_Impl&);
    friend TaylorModel_Impl* asinh  (const TaylorModel_Impl&);
    friend TaylorModel_Impl* acosh  (const TaylorModel_Impl&);
    friend TaylorModel_Impl* atanh  (const TaylorModel_Impl&);
    friend TaylorModel_Impl* acoth  (const TaylorModel_Impl&);
    friend TaylorModel_Impl* pow    (const TaylorModel_Impl&,const TaylorModel_Impl&);
    friend TaylorModel_Impl* power  (const TaylorModel_Impl&,int);

    friend TaylorModel_Impl* invert    (const TaylorModel_Impl&);
    friend TaylorModel_Impl* integrate (const TaylorModel_Impl&,unsigned int);
    friend TaylorModel_Impl* derivate  (const TaylorModel_Impl&,unsigned int);
    friend TaylorModel_Impl* substitute(const TaylorModel_Impl&,unsigned int,const Interval&);

    friend std::ostream& operator <<  (std::ostream&, const TaylorModel_Impl&);
    friend std::string&  operator <<  (std::string& , const TaylorModel_Impl&);

    /*
      Static functions following
    */

    //
    // Function to create a constant polynomial with zero remainder interval
    //
    static TaylorModel_Impl* Const_TM(const double&);

    //
    // Functions to read and/or change the static members
    //
    static unsigned int order()
      {
        return order_;
      }
    static unsigned int set_order(unsigned int n)
      {
        unsigned int old = order_;
        order_ = n;
        return old;
      }
    static double set_sparsity_tol(const double& tol)
      {
        if( tol > 1e-2 ) return sparsity_tol_;

        double   old  = sparsity_tol_;
        sparsity_tol_ = tol;

        return old;
      }
    static void set_polynomial_range_evaluation(const PolynomialRange& pb)
      {
        poly_eval_ = &pb;
      }
    static void set_degree_check(const DegreeBase& db)
      {
        degree_ = &db;
      }
    static void set_mode    (unsigned int);
    static void mode_message(unsigned int);

  private:

    /*
      Functions for internal use
    */

    //
    // A constructor to create a TaylorModel_Impl-Object with an empty polynom
    // from an existing TaylorModel_Impl-Object
    //
    TaylorModel_Impl(unsigned long poly_size, const Interval& rest, const additional_var_data& data)
      : polynomial_( poly_size, MonomHasher2::getFunction() ), rest_(rest), data_(data)
      {}

    static bool check_var_data(additional_var_data& x, const additional_var_data& y);

    //
    // Two functions each creating an array of counters. The array implementation
    // makes it possible to break the evaluation sequence and start the sequence from beginning
    //
    static unsigned int* init_play_cntr  (unsigned int);
    static unsigned int* init_record_cntr(unsigned int);

    //
    // Arithmetic operators for internal use
    //
    TaylorModel_Impl& add_without_seizing_rounding_errors     (const TaylorModel_Impl&);
    TaylorModel_Impl& add_with_seizing_rounding_errors        (const TaylorModel_Impl&);
    TaylorModel_Impl& subtract_without_seizing_rounding_errors(const TaylorModel_Impl&);
    TaylorModel_Impl& subtract_with_seizing_rounding_errors   (const TaylorModel_Impl&);
    TaylorModel_Impl& multiply_without_seizing_rounding_errors(const TaylorModel_Impl&);
    TaylorModel_Impl& multiply_with_seizing_rounding_errors   (const TaylorModel_Impl&);

    TaylorModel_Impl& add_without_seizing_rounding_errors     (const Interval&);
    TaylorModel_Impl& add_with_seizing_rounding_errors        (const Interval&);
    TaylorModel_Impl& subtract_without_seizing_rounding_errors(const Interval&);
    TaylorModel_Impl& subtract_with_seizing_rounding_errors   (const Interval&);
    TaylorModel_Impl& multiply_without_seizing_rounding_errors(const Interval&);
    TaylorModel_Impl& multiply_with_seizing_rounding_errors   (const Interval&);
    TaylorModel_Impl& divide_without_seizing_rounding_errors  (const Interval&);
    TaylorModel_Impl& divide_with_seizing_rounding_errors     (const Interval&);

    TaylorModel_Impl& add_without_seizing_rounding_errors     (const double&);
    TaylorModel_Impl& add_with_seizing_rounding_errors        (const double&);
    TaylorModel_Impl& subtract_without_seizing_rounding_errors(const double&);
    TaylorModel_Impl& subtract_with_seizing_rounding_errors   (const double&);
    TaylorModel_Impl& multiply_without_seizing_rounding_errors(const double&);
    TaylorModel_Impl& multiply_with_seizing_rounding_errors   (const double&);
    TaylorModel_Impl& divide_without_seizing_rounding_errors  (const double&);
    TaylorModel_Impl& divide_with_seizing_rounding_errors     (const double&);

    /*
      Data
    */

    polynomial_type     polynomial_; //Hashtable.
    Interval            rest_;
    additional_var_data data_;       //Contains the development point and the domain
    //interval of each variable.

    /*
      Data for all Taylor-Models
    */

    static unsigned int           order_;
    static double                 sparsity_tol_;
    static const PolynomialRange *poly_eval_;
    static const DegreeBase      *degree_;

    static unsigned int           Mode;

    static Variables             *vtab_;
  };
}

#endif

/*

  End of File: taylormodel_impl.h

*/
