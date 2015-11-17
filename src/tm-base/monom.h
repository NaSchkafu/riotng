/*

  File: monom.h, 2004/08/18

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

#ifndef MONOM_H_INCLUDED
#define MONOM_H_INCLUDED

#include <string>
#include <iostream>
#include "vartable.h"

/*

  Monom verwaltet derzeit die Information über MAXIMAL 8 Variablen.
  Erweiterung sind leicht möglich.

*/
namespace riot
{
  class Monom
  {
  public:

    /*
      Constructors
    */

  Monom() : block_1to4_(0), block_5to8_(0), b_sum_(0), left_(0), right_(0)
    {
      b_array_[0] = (unsigned char*)&block_1to4_;
      b_array_[1] = (unsigned char*)&block_5to8_;
    }
    Monom(unsigned int, unsigned int = 1);
    Monom(const std::string&, unsigned int = 1);
  Monom(const Monom& x)
    : block_1to4_(x.block_1to4_), block_5to8_(x.block_5to8_), b_sum_(x.b_sum_), left_(x.left_), right_(x.right_)
    {
      b_array_[0] = (unsigned char*)&block_1to4_;
      b_array_[1] = (unsigned char*)&block_5to8_;
    }
    ~Monom() {}

    Monom& operator = (const Monom&);

    /*
      State
    */

    bool is_const() const             { return b_sum_ == 0; }
    unsigned int min_var_code() const { return left_;  }
    unsigned int max_var_code() const { return right_; }
    unsigned int expnt_of_var(unsigned int i) const
    {
      if( right_ > 0 )
      {
        --i;
        return (b_array_[ i / 4 ])[ i % 4 ]; // 'i / 4' gibt an in welchem 4er-Block die Var. ihren Platz hat ...
        // und 'i % 4' gibt die Stelle innerhalb des 4er-Blocks an.
      }
      else return 0;
    }
    unsigned int expnt_sum() const { return b_sum_; }

    /*
      Arithmetic operation
    */

    Monom& operator *= (const Monom&);

    unsigned int integrate_regarding_to(unsigned int);
    unsigned int derivate_regarding_to (unsigned int);
    unsigned int substitute            (unsigned int);

    bool operator == (const Monom& x) const
    {
      return (block_1to4_ == x.block_1to4_ && block_5to8_ == x.block_5to8_);
    }
    bool operator != (const Monom& x) const
    {
      return (block_1to4_ != x.block_1to4_ || block_5to8_ != x.block_5to8_);
    }
    bool operator <  (const Monom& x) const
    {
      return (block_5to8_ < x.block_5to8_) || (block_5to8_ == x.block_5to8_ && block_1to4_ < x.block_1to4_);
    }
    bool operator >  (const Monom& x) const
    {
      return (block_5to8_ > x.block_5to8_) || (block_5to8_ == x.block_5to8_ && block_1to4_ > x.block_1to4_);
    }

    /*
      Friend functions
    */

    friend std::ostream& operator << (std::ostream&,const Monom&);
    friend std::string&  operator << (std::string& ,const Monom&);

  private:

    unsigned long int block_1to4_, block_5to8_;
    unsigned int      b_sum_;
    unsigned char    *b_array_[2];

    unsigned int left_, right_;

    static Variables *vtab_;
  };

  inline Monom operator * (const Monom& x, const Monom& y)
  {
    Monom z(x);
    return z *= y;
  }
}

#endif //Header

/*

  End of File: monom.h

*/
