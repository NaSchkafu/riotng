/*

  File: monom.cpp, 2004/08/18

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

#include "monom.h"

/*
  Set the pointer to the Table of Variables.
*/
namespace riot
{

  Variables *Monom::vtab_ = &Variables::getTable();

/*
  Constructor.
*/

  Monom::Monom(unsigned int code, unsigned int e)
    : block_1to4_(0), block_5to8_(0), b_sum_(0), left_(0), right_(0)
  {
    b_array_[0] = (unsigned char*)&block_1to4_;
    b_array_[1] = (unsigned char*)&block_5to8_;

    if( 9 > code && code > 0 && e > 0 ) //<------- Beschränkung auf 8 Var.
    {
      left_ = right_ = code;
      --code;
      b_sum_ = e;

      (b_array_[ code / 4 ])[ code % 4 ] = e; // 'code / 4' gibt an in welchem 4er-Block die Var. ihren Platz hat ...
                                              // und 'code % 4' gibt die Stelle innerhalb des 4er-Blocks an.
    }
  }

  Monom::Monom(const std::string& vname, unsigned int e)
    : block_1to4_(0), block_5to8_(0), b_sum_(0), left_(0), right_(0)
  {
    b_array_[0] = (unsigned char*)&block_1to4_;
    b_array_[1] = (unsigned char*)&block_5to8_;

    int code = vtab_->insert( vname );  //Var.code holen.

    if( 9 > code && code > 0 && e > 0 ) //<------- Beschränkung auf 8 Var.
    {
      left_ = right_ = code;
      --code;
      b_sum_ = e;

      (b_array_[ code / 4 ])[ code % 4 ] = e; // 'code / 4' gibt an in welchem 4er-Block die Var. ihren Platz hat ...
                                              // und 'code % 4' gibt die Stelle innerhalb des 4er-Blocks an.
    }
  }

/*
  The assignment operator.
*/

  Monom& Monom::operator = (const Monom& x)
  {
    block_1to4_ = x.block_1to4_;
    block_5to8_ = x.block_5to8_;
    b_array_[0] = (unsigned char*)&block_1to4_;
    b_array_[1] = (unsigned char*)&block_5to8_;
    b_sum_      = x.b_sum_;
    left_       = x.left_;
    right_      = x.right_;

    return *this;
  }

/*
  The monom multiplication.
*/

  Monom& Monom::operator *= (const Monom& m)
  {
    block_1to4_ += m.block_1to4_;
    block_5to8_ += m.block_5to8_;

    if( m.b_sum_ == 0 ) return *this;
    else if( b_sum_ == 0 )
    {
      left_  = m.left_;
      right_ = m.right_;
    }
    else
    {
      left_  = (left_  < m.left_ )?left_ :m.left_ ;
      right_ = (right_ > m.right_)?right_:m.right_;
    }

    b_sum_ += m.b_sum_;

    return *this;
  }

/*
  Integration regarding to given variable. Returns new exponent.
*/

  unsigned int Monom::integrate_regarding_to(unsigned int code)
  {
    if( 9 > code && code > 0 )
    {
      if( code < left_  ) left_  = code;
      if( code > right_ )
        if( right_ == 0 ) left_  = right_ = code;
        else              right_ = code;

      --code;
      ++b_sum_;
      return ( (b_array_[ code / 4 ])[ code % 4 ] += 1 ); //Return new exponent.
    }
    else return code;
  }

/*
  Derivate regarding to given variable. Returns old exponent.
*/

  unsigned int Monom::derivate_regarding_to (unsigned int code)
  {
    if( 9 > code && code > 0 )
    {
      if( code < left_ || code > right_ )
      {
        return (left_ = right_ = b_sum_ = block_1to4_ = block_5to8_ = 0);
      }
      else // left_ <= code <= right_
      {
        --code;
        unsigned int exponent = (b_array_[ code / 4 ])[ code % 4 ];
        --b_sum_;

        if( exponent == 0 || ( exponent == 1 && left_ == right_ ) )
        {
          left_ = right_ = b_sum_ = block_1to4_ = block_5to8_ = 0;
        }
        else if( exponent == 1 ) //and left_ != right_.
        {
          (b_array_[ code / 4 ])[ code % 4 ] -= 1;

          ++code;

          if( code == left_ )
          {
            unsigned int new_left = left_;

            //Search for next variable with exponent != 0.
            while( (b_array_[ new_left / 4 ])[ new_left % 4 ] == 0 ) { ++new_left; }

            left_ = ++new_left;     //Set new lower bound.
          }
          else if( code == right_ )
          {
            unsigned int new_right = right_-2;

            //Search for next variable with exponent != 0.
            while( (b_array_[ new_right / 4 ])[ new_right % 4 ] == 0 ) { --new_right; }

            right_ = ++new_right;   //Set new upper bound.
          }
        }
        else (b_array_[ code / 4 ])[ code % 4 ] -= 1;

        return exponent;
      }
    }
    else return code;
  }

/*
  Substitute a given variable. Returns exponent.
*/

  unsigned int Monom::substitute(unsigned int code)
  {
    if( 9 > code && code > 0 )
    {
      if     ( code < left_ || right_ < code ) return 0; //No variable with code 'code'.
      else if( left_ == right_ ) //code = left_ = right_.
        left_ = right_ = 0;
      else
      {
        if( code == left_ )
        {
          unsigned int new_left = left_;

          //Search for next variable with exponent != 0.
          while( (b_array_[ new_left / 4 ])[ new_left % 4 ] == 0 ) { ++new_left; }

          left_ = ++new_left;     //Set new lower bound.
        }
        else if( code == right_ )
        {
          unsigned int new_right = right_-2;

          //Search for next variable with exponent != 0.
          while( (b_array_[ new_right / 4 ])[ new_right % 4 ] == 0 ) { --new_right; }

          right_ = ++new_right;   //Set new upper bound.
        }
      }

      --code;
      unsigned int exponent = (b_array_[ code / 4 ])[ code % 4 ];
      if( exponent > 0 ) { (b_array_[ code / 4 ])[ code % 4 ] = 0; b_sum_ -= exponent; }
      return exponent; //Return exponent.
    }
    else return code;
  }

/*
  Output operator.
*/

  std::ostream& operator << (std::ostream& os,const Monom& x)
  {
    if( x.b_sum_ == 0 ) return os;
    else
    {
      unsigned int e;
      for(int i = 0; i < 8; i++)
      {
        e = (x.b_array_[ i / 4 ])[ i % 4 ];
        if( e > 0 )
        {
          os << "*" << (Monom::vtab_->look_up(i+1));

          if( e > 1 ) os << "^" << e;
        }
      }

      return os;
    }
  }

  std::string&  operator << (std::string&  st,const Monom& x)
  {
    std::cerr << "string << Monom: Not yet implemented" << std::endl;
    return st;
  }
}

/*

  End of File: monom.cpp

*/
