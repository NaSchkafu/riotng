/*

  File: polyeval.cpp, 2004/11/15

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

#include "polyeval.h"
#include <cmath>
namespace riot
{

  Interval IntervalEvaluation::operator ()(const IntervalEvaluation::polynomial_type     *poly_ptr,
                                           const IntervalEvaluation::additional_var_data *data_ptr  ) const
  {
    Interval bound(0.0);

    std::vector<Interval> arguments( data_ptr->size() );
    for(unsigned int i = 0; i < data_ptr->size(); i++)
      if( (*data_ptr)[ i ].operator ->() != 0 )
        arguments[ i ] = ((*data_ptr)[ i ])->domain_ - ((*data_ptr)[ i ])->devel_point_;

    IntervalEvaluation::const_polynomial_iterator
      curr = poly_ptr->begin(),
      last = poly_ptr->end();

    while( curr != last ) //Walk trough the polynomial.
    {
      Interval prod(1.0);

      //Evaluate the current monomial.
      for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
      {
        int e = (*curr).key().expnt_of_var(i);
        if( e != 0 )
        {
          prod *= power( arguments[i-1], e );
        }
      }

      //Multiply with coefficient.
      prod *= (*curr).value();   // <---- Multiplication with double value.

      //Add enclosure to bound interval.
      bound += prod;

      ++curr;
    }

    return bound;
  }
  Interval IntervalEvaluation::operator ()(const IntervalEvaluation::ipolynomial_type    *poly_ptr,
                                           const IntervalEvaluation::additional_var_data *data_ptr  ) const
  {
    Interval bound(0.0);

    std::vector<Interval> arguments( data_ptr->size() );
    for(unsigned int i = 0; i < data_ptr->size(); i++)
      if( (*data_ptr)[ i ].operator ->() != 0 )
        arguments[ i ] = ((*data_ptr)[ i ])->domain_ - ((*data_ptr)[ i ])->devel_point_;

    IntervalEvaluation::const_ipolynomial_iterator
      curr = poly_ptr->begin(),
      last = poly_ptr->end();

    while( curr != last ) //Walk trough the polynomial.
    {
      Interval prod(1.0);

      //Evaluate the current monomial.
      for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
      {
        int e = (*curr).key().expnt_of_var(i);
        if( e != 0 )
        {
          prod *= power( arguments[i-1], e );
        }
      }

      //Multiply with coefficient.
      prod *= (*curr).value();   // <---- Multiplication with interval value.

      //Add enclosure to bound interval.
      bound += prod;

      ++curr;
    }

    return bound;
  }

  Interval MeanValueForm::operator ()(const MeanValueForm::polynomial_type     *poly_ptr,
                                      const MeanValueForm::additional_var_data *data_ptr  ) const
  {
    Interval bound1(0.0), bound2(0.0);
    std::vector<Interval> bound_of_derivatives( data_ptr->size(), Interval(0.0) );

    std::vector<Interval> arguments( data_ptr->size() ), arguments_mid( data_ptr->size() );
    for(unsigned int i = 0; i < data_ptr->size(); i++)
      if( (*data_ptr)[ i ].operator ->() != 0 )
      {
        arguments    [ i ] = ((*data_ptr)[ i ])->domain_ - ((*data_ptr)[ i ])->devel_point_;
        arguments_mid[ i ] = Interval( arguments[ i ].mid() );
      }

    MeanValueForm::const_polynomial_iterator
      curr = poly_ptr->begin(),
      last = poly_ptr->end();

    while( curr != last ) //Walk trough the polynomial.
    {
      Interval prod1(1.0), prod2(1.0);

      //Evaluate the current monomial at the midpoint of the domain.
      for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
      {
        int e = (*curr).key().expnt_of_var(i);
        if( e != 0 )
        {
          prod1 *= power( arguments    [i-1], e );
          prod2 *= power( arguments_mid[i-1], e );
        }
      }

      //Multiply with coefficient.
      prod1 *= (*curr).value();   // <---- Multiply with double value.
      prod2 *= (*curr).value();   // <---- Multiply with double value.

      //Add enclosure to bound interval.
      bound1 += prod1;
      bound2 += prod2;

      //Derivate the current monomial and evaluate it.
      for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
      {
        unsigned int e_i = (*curr).key().expnt_of_var(i);
        if( e_i != 0 )
        {
          prod2 = e_i * power( arguments[i-1], e_i-1 ); //Value of derivative.

          unsigned int j;
          for(j = (*curr).key().min_var_code(); j < i; j++)
          {
            unsigned int e_j = (*curr).key().expnt_of_var(j);
            if( e_j != 0 )
            {
              prod2 *= power( arguments[j-1], e_j );
            }
          }
          for(j = i+1; j <= (*curr).key().max_var_code(); j++)
          {
            unsigned int e_j = (*curr).key().expnt_of_var(j);
            if( e_j != 0 )
            {
              prod2 *= power( arguments[j-1], e_j );
            }
          }

          //Multiply with coefficient.
          prod2 *= (*curr).value();   // <---- Multiply with double value.

          //Save enclosure of derivative.
          bound_of_derivatives[ i - 1 ] += prod2;
        }
      }

      ++curr;
    }

    for(unsigned int i = 0; i < bound_of_derivatives.size(); i++)
      if( bound_of_derivatives[i] != Interval(0.0) )
        bound2 += bound_of_derivatives[i] * (((*data_ptr)[i])->domain_ - (((*data_ptr)[i])->domain_).mid() );

    return bound1 & bound2; //Intersection of naive and mean value form evaluation.
  }
  Interval MeanValueForm::operator ()(const MeanValueForm::ipolynomial_type    *poly_ptr,
                                      const MeanValueForm::additional_var_data *data_ptr  ) const
  {
    Interval bound1(0.0), bound2(0.0);
    std::vector<Interval> bound_of_derivatives( data_ptr->size(), Interval(0.0) );

    std::vector<Interval> arguments( data_ptr->size() ), arguments_mid( data_ptr->size() );
    for(unsigned int i = 0; i < data_ptr->size(); i++)
      if( (*data_ptr)[ i ].operator ->() != 0 )
      {
        arguments    [ i ] = ((*data_ptr)[ i ])->domain_ - ((*data_ptr)[ i ])->devel_point_;
        arguments_mid[ i ] = Interval( arguments[ i ].mid() );
      }

    MeanValueForm::const_ipolynomial_iterator
      curr = poly_ptr->begin(),
      last = poly_ptr->end();

    while( curr != last ) //Walk trough the polynomial.
    {
      Interval prod1(1.0), prod2(1.0);

      //Evaluate the current monomial at the midpoint of the domain.
      for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
      {
        unsigned int e = (*curr).key().expnt_of_var(i);
        if( e != 0 )
        {
          prod1 *= power( arguments    [i-1], e );
          prod2 *= power( arguments_mid[i-1], e );
        }
      }

      //Multiply with coefficient.
      prod1 *= (*curr).value();   // <---- Multiply with double value.
      prod2 *= (*curr).value();   // <---- Multiply with double value.

      //Add enclosure to bound interval.
      bound1 += prod1;
      bound2 += prod2;

      //Derivate the current monomial and evaluate it.
      for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
      {
        unsigned int e_i = (*curr).key().expnt_of_var(i);
        if( e_i != 0 )
        {
          prod2 = e_i * power( arguments[i-1], e_i-1 ); //Value of derivative.

          unsigned int j;
          for(j = (*curr).key().min_var_code(); j < i; j++)
          {
            unsigned int e_j = (*curr).key().expnt_of_var(j);
            if( e_j != 0 )
            {
              prod2 *= power( arguments[j-1], e_j );
            }
          }
          for(j = i+1; j <= (*curr).key().max_var_code(); j++)
          {
            unsigned int e_j = (*curr).key().expnt_of_var(j);
            if( e_j != 0 )
            {
              prod2 *= power( arguments[j-1], e_j );
            }
          }

          //Multiply with coefficient.
          prod2 *= (*curr).value();   // <---- Multiply with double value.

          //Save enclosure of derivative.
          bound_of_derivatives[ i - 1 ] += prod2;
        }
      }

      ++curr;
    }

    for(unsigned int i = 0; i < bound_of_derivatives.size(); i++)
      if( bound_of_derivatives[i] != Interval(0.0) )
        bound2 += bound_of_derivatives[i] * (((*data_ptr)[i])->domain_ - (((*data_ptr)[i])->domain_).mid() );

    return bound1 & bound2; //Intersection of naive and mean value form evaluation.
  }

  Interval LDB::operator ()(const LDB::polynomial_type     *poly_ptr ,
                            const LDB::additional_var_data *data_ptr   ) const
  {
    unsigned int data_size = data_ptr->size();

    /*
      Abspalten der linearen Terme. Die Konstante wird den hoeheren Termen zugerechnet.
      Gleichzeitig werden die hoeheren Terme eingeschlossen.
    */

    Interval rest(0.0);

    std::vector<Interval> arguments( data_ptr->size() );
    for(unsigned int i = 0; i < data_size; i++)
      if( (*data_ptr)[ i ].operator ->() != 0 )
        arguments[ i ] = ((*data_ptr)[ i ])->domain_ - ((*data_ptr)[ i ])->devel_point_;

    std::vector<double>       linear_coeffs(data_size,0.0);
    std::vector<unsigned int> var_codes; var_codes.reserve(data_size); //Codes of Vars with linear term.

    LDB::polynomial_type::const_iterator
      curr = poly_ptr->begin(),
      last = poly_ptr->end();

    while( curr != last )
    {
      if( curr->key().expnt_sum() == 0 )      //Konstanter Term.
      {
        rest += curr->value();
      }
      else if( curr->key().expnt_sum() == 1 ) //Linearer Term.
      {
        unsigned int i = curr->key().min_var_code();

        var_codes.push_back(i-1);
        linear_coeffs[i-1] = curr->value();
      }
      else //Nonlinear term with more than one variable.
      {
        Interval prod(1.0);

        //Evaluate the current monomial.
        for(unsigned int i = curr->key().min_var_code(); i <= curr->key().max_var_code(); i++)
        {
          unsigned int e = curr->key().expnt_of_var(i);
          if( e != 0 )
          {
            prod *= power( arguments[i-1], e );
          }
        }

        //Multiply with coefficient.
        prod *= curr->value();

        //Add enclosure to bound interval.
        rest += prod;
      }

      ++curr;
    }

    if( var_codes.size() == 0 ) return rest; //No linear terms in polynomial. So return with enclosure of
    //higher order terms.

    /*
      Lokale Kopien der 'Domains' anlegen, da diese im weiteren Verlauf unter Umstaenden
      veraendert werden.
    */

    std::vector<Interval> domain_of_max(data_size,Interval(0.0));
    std::vector<Interval> domain_of_min(data_size,Interval(0.0));

    for(unsigned int i = 0; i < data_size; i++)
    {
      if( (*data_ptr)[i].operator ->() ) //There is information.
      {
        domain_of_max[i] = domain_of_min[i] = ((*data_ptr)[i])->domain_;
      }
    }

    /*
      Jetzt beginnt der Algorithmus des LDB range bounders.
    */

    double delta = rest.diam();
    bool resizing = false;
    for(unsigned int i = 0; i < var_codes.size(); i++)
    {
      unsigned int k = var_codes[i];

      Interval abs_b  = Interval( std::abs(linear_coeffs[k]) );
      Interval domain = ((*data_ptr)[k])->domain_;
      Interval upper = abs_b * domain.diam();
      if( upper.inf() > delta )
      {
        Interval tmp = delta / abs_b;
        if( linear_coeffs[k] > 0 )
        {
          domain_of_min[k] = domain_of_min[k].inf() + Interval(0,tmp.sup());
          domain_of_max[k] = domain_of_max[k].sup() - Interval(0,tmp.sup());
        }
        else
        {
          domain_of_min[k] = domain_of_min[k].sup() - Interval(0,tmp.sup());
          domain_of_max[k] = domain_of_max[k].inf() + Interval(0,tmp.sup());
        }
        resizing = true;
      }
    }

    for(unsigned int i = 0; i < data_size; i++)
    {
      if( (*data_ptr)[i].operator ->() ) //There is information.
      {
        domain_of_min[i] -= ((*data_ptr)[i])->devel_point_;
        domain_of_max[i] -= ((*data_ptr)[i])->devel_point_;
      }
    }

    if( resizing ) //Die Domains wurden verändert.
    {
      Interval min(0.0),max(0.0);

      curr = poly_ptr->begin();

      while( curr != last ) //Walk trough the polynomial.
      {
        Interval prod1(1.0),prod2(1.0);

        //Evaluate the current monomial.
        for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
        {
          unsigned int e = (*curr).key().expnt_of_var(i);
          if( e != 0 )
          {
            prod1 *= power( domain_of_min[ i - 1 ], e );
            prod2 *= power( domain_of_max[ i - 1 ], e );
          }
        }

        //Multiply with coefficient.
        prod1 *= (*curr).value();   // <---- Multiply with double value.
        prod2 *= (*curr).value();   // <---- Multiply with double value.

        //Add enclosure to interval.
        min += prod1;
        max += prod2;

        ++curr;
      }

      return Interval( min.inf(), max.sup() );
    }
    else
    {
      for(unsigned int i = 0; i < var_codes.size(); i++)
      {
        unsigned int k = var_codes[i];
        rest += linear_coeffs[k] * domain_of_min[k];
      }
      return rest;
    }
  }
  Interval LDB::operator ()(const LDB::ipolynomial_type    *poly_ptr,
                            const LDB::additional_var_data *data_ptr  ) const
  {
    unsigned int data_size = data_ptr->size();

    /*
      Abspalten der linearen Terme. Die Konstante wird den hoeheren Termen zugerechnet.
      Gleichzeitig werden die hoeheren Terme eingeschlossen.
    */

    Interval rest(0.0);
    std::vector<Interval> arguments( data_ptr->size() );
    for(unsigned int i = 0; i < data_size; i++)
      if( (*data_ptr)[ i ].operator ->() != 0 )
        arguments[ i ] = ((*data_ptr)[ i ])->domain_ - ((*data_ptr)[ i ])->devel_point_;

    std::vector<Interval> linear_coeffs(data_size, Interval(0.0));
    std::vector<unsigned int> var_codes; var_codes.reserve(data_size); //Codes of Vars with linear term.

    LDB::ipolynomial_type::const_iterator
      curr = poly_ptr->begin(),
      last = poly_ptr->end();

    while( curr != last )
    {
      if( curr->key().expnt_sum() == 0 )      //Konstanter Term.
      {
        rest += curr->value();
      }
      else if( curr->key().expnt_sum() == 1 ) //Linearer Term.
      {
        unsigned int i = curr->key().min_var_code();

        var_codes.push_back(i-1);
        linear_coeffs[i-1] = curr->value();
      }
      else //Nonlinear term.
      {
        Interval prod(1.0);

        //Evaluate the current monomial.
        for(unsigned int i = curr->key().min_var_code(); i <= curr->key().max_var_code(); i++)
        {
          unsigned int e = curr->key().expnt_of_var(i);
          if( e != 0 )
          {
            prod *= power( arguments[i-1], e );
          }
        }
        //Multiply with coefficient.
        prod *= curr->value();

        //Add enclosure to bound interval.
        rest += prod;
      }

      ++curr;
    }

    if( var_codes.size() == 0 ) return rest; //No linear terms in polynomial. So return with enclosure of
    //higher order terms.

    /*
      Lokale Kopien der 'Domains' anlegen, da diese im weiteren Verlauf unter Umstaenden
      veraendert werden.
    */

    std::vector<Interval> domain_of_max(data_size,Interval(0.0));
    std::vector<Interval> domain_of_min(data_size,Interval(0.0));

    for(unsigned int i = 0; i < data_size; i++)
    {
      if( (*data_ptr)[i].operator ->() ) //There is information.
      {
        domain_of_max[i] = domain_of_min[i] = ((*data_ptr)[i])->domain_;
      }
    }

    /*
      Jetzt beginnt der Algorithmus des LDB range bounders.
    */

    Interval idelta(rest.diam());
    for(unsigned int i = 0; i < var_codes.size(); i++)
    {
      unsigned int k = var_codes[i];
      Interval b = linear_coeffs[k];
      idelta += b.diam() * abs( ((*data_ptr)[k])->domain_ - ((*data_ptr)[k])->devel_point_ );
    }
    double delta = idelta.sup();

    bool resizing = false;
    for(unsigned int i = 0; i < var_codes.size(); i++)
    {
      unsigned int k = var_codes[i];

      double   mid_b  = linear_coeffs[k].mid();
      Interval abs_b  = Interval( std::abs(mid_b) );
      Interval domain = ((*data_ptr)[k])->domain_;
      Interval upper  = abs_b * domain.diam();
      if( upper.inf() > delta )
      {
        Interval tmp = delta / abs_b;
        if( mid_b > 0 )
        {
          domain_of_min[k] = domain_of_min[k].inf() + Interval(0,tmp.sup());
          domain_of_max[k] = domain_of_max[k].sup() - Interval(0,tmp.sup());
        }
        else
        {
          domain_of_min[k] = domain_of_min[k].sup() - Interval(0,tmp.sup());
          domain_of_max[k] = domain_of_max[k].inf() + Interval(0,tmp.sup());
        }
        resizing = true;
      }
    }

    for(unsigned int i = 0; i < data_size; i++)
    {
      if( (*data_ptr)[i].operator ->() ) //There is information.
      {
        domain_of_min[i] -= ((*data_ptr)[i])->devel_point_;
        domain_of_max[i] -= ((*data_ptr)[i])->devel_point_;
      }
    }

    if( resizing ) //Die Domains wurden verändert.
    {
      Interval min(0.0),max(0.0);

      curr = poly_ptr->begin();

      while( curr != last ) //Walk trough the polynomial.
      {
        Interval prod1(1.0),prod2(1.0);

        //Evaluate the current monomial.
        for(unsigned int i = (*curr).key().min_var_code(); i <= (*curr).key().max_var_code(); i++)
        {
          unsigned int e = (*curr).key().expnt_of_var(i);
          if( e != 0 )
          {
            prod1 *= power( domain_of_min[ i - 1 ], e );
            prod2 *= power( domain_of_max[ i - 1 ], e );
          }
        }

        //Multiply with coefficient.
        prod1 *= (*curr).value();   // <---- Multiply with double value.
        prod2 *= (*curr).value();   // <---- Multiply with double value.

        //Add enclosure to interval.
        min += prod1;
        max += prod2;

        ++curr;
      }

      return Interval( min.inf(), max.sup() );
    }
    else
    {
      for(unsigned int i = 0; i < var_codes.size(); i++)
      {
        unsigned int k = var_codes[i];
        rest += linear_coeffs[k] * domain_of_min[k];
      }
      return rest;
    }
  }

  Interval BernsteinForm::operator ()(const BernsteinForm::polynomial_type     *poly_ptr,
                                      const BernsteinForm::additional_var_data *data_ptr  ) const
  {
    std::cerr<<"BernsteinForm: Noch nicht implementiert"<<std::endl;
    Interval bound(0.0);
    return bound;
  }
  Interval BernsteinForm::operator ()(const BernsteinForm::ipolynomial_type    *poly_ptr,
                                      const BernsteinForm::additional_var_data *data_ptr  ) const
  {
    std::cerr<<"BernsteinForm: Noch nicht implementiert"<<std::endl;
    Interval bound(0.0);
    return bound;
  }
}

/*

  End of File: polyeval.cpp

*/
