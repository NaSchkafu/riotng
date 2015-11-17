/*

 File: polyeval.h, 2004/11/15

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

#ifndef POLYEVAL_H_INCLUDED
#define POLYEVAL_H_INCLUDED

#include "adaptintval.h"
#include "taylormodel_impl.h"

using Adapt::Interval;

namespace riot 
{
  

class PolynomialRange
{
 public:

  // Introduce better readable names

  typedef TaylorModel_Impl::polynomial_type             polynomial_type;
  typedef TaylorModel_Impl::polynomial_iterator         polynomial_iterator;
  typedef TaylorModel_Impl::const_polynomial_iterator   const_polynomial_iterator;

  typedef TaylorModel_Impl::ipolynomial_type            ipolynomial_type;
  typedef TaylorModel_Impl::ipolynomial_iterator        ipolynomial_iterator;
  typedef TaylorModel_Impl::const_ipolynomial_iterator  const_ipolynomial_iterator;

  typedef TaylorModel_Impl::additional_var_data         additional_var_data;
  
  // Virtual functions

  virtual ~PolynomialRange() {}

  virtual Interval operator ()(const polynomial_type* ,const additional_var_data*) const = 0;
  virtual Interval operator ()(const ipolynomial_type*,const additional_var_data*) const = 0;
};

/*
  The IntervalEvaluation is a naive evaluation of the polynomial. Every term gets separatly enclosed,
  then all enclosures will be added.
*/

class IntervalEvaluation : public PolynomialRange
{
 private:
  
  /*
    To make it a Singleton.
  */

  IntervalEvaluation() {}
  IntervalEvaluation(const IntervalEvaluation&);
  ~IntervalEvaluation() {}

  IntervalEvaluation& operator = (const IntervalEvaluation&);

 public:

  static IntervalEvaluation& getObject()
    {
      static IntervalEvaluation singleton;
      return singleton;
    }

  Interval operator ()(const polynomial_type* ,const additional_var_data*) const;
  Interval operator ()(const ipolynomial_type*,const additional_var_data*) const;
};

/*
  The MeanValueForm is the evaluation of the polynomial according to the well known
  mean value form.
*/

class MeanValueForm : public PolynomialRange
{
 private:
  
  /*
    To make it a Singleton.
  */

  MeanValueForm() {}
  MeanValueForm(const MeanValueForm&);
  ~MeanValueForm() {}

  MeanValueForm& operator = (const MeanValueForm&);

 public:

  static MeanValueForm& getObject()
    {
      static MeanValueForm singleton;
      return singleton;
    }

  Interval operator ()(const polynomial_type* ,const additional_var_data*) const;
  Interval operator ()(const ipolynomial_type*,const additional_var_data*) const;
};

/*
  The LDB is the Linear Dominated Bounder.
*/

class LDB : public PolynomialRange
{
 private:
  
  /*
    To make it a Singleton.
  */

  LDB() {}
  LDB(const LDB&);
  ~LDB() {}

  LDB& operator = (const LDB&);

 public:

  static LDB& getObject()
    {
      static LDB singleton;
      return singleton;
    }

  Interval operator ()(const polynomial_type* ,const additional_var_data*) const;
  Interval operator ()(const ipolynomial_type*,const additional_var_data*) const;
};

/*
  The BernsteinForm computes bounds of the polynomial using the
  Bernstein form.
*/

class BernsteinForm : public PolynomialRange
{
 private:
  
  /*
    To make it a Singleton.
  */

  BernsteinForm() {}
  BernsteinForm(const BernsteinForm&);
  ~BernsteinForm() {}

  BernsteinForm& operator = (const BernsteinForm&);

 public:

  static BernsteinForm& getObject()
    {
      static BernsteinForm singleton;
      return singleton;
    }

  Interval operator ()(const polynomial_type* ,const additional_var_data*) const;
  Interval operator ()(const ipolynomial_type*,const additional_var_data*) const;
};

#endif
}

/*

End of File: polyeval.h

*/
