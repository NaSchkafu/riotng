/*

 File: degree.h, 2004/12/28

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

#ifndef DEGREE_H_INCLUDED
#define DEGREE_H_INCLUDED

#include "degreeeval.h"

namespace riot{

class TotalDegree_TimeSeparate : public DegreeBase, public Observer //The time variable plays an extra part.
{
 public:

  static TotalDegree_TimeSeparate& getObject()
    {
      static TotalDegree_TimeSeparate singleton;
      return singleton;
    }

  void update(const Subject*, const Action* a) 
    {
      if     ( a == &NewVariable::getObject() ) //New variable.
	{
	  if( code_ == 0 && rep_ == NewVariable::getObject().identifier_ ) code_ = NewVariable::getObject().code_;
	}
      else if( a == &IdentifierChanged::getObject() ) //Change of identifier. Maybe useful in the future.
	{
	  if     ( rep_ == IdentifierChanged::getObject().old_identifier_ ) rep_  = IdentifierChanged::getObject().new_identifier_;
	  else if( rep_ == IdentifierChanged::getObject().new_identifier_ ) code_ = IdentifierChanged::getObject().code_;
	}
      else //Confirmation of announcement.
	{
	  code_ = Variables::getTable().look_up( rep_ );
	}
    }

  void change_time_identifier(const std::string& s) 
    { 
      rep_  = s;
      code_ = Variables::getTable().look_up(s);
    }

  bool check(const Monom* m, unsigned int max_deg) const
    {
      unsigned int time_expnt = m->expnt_of_var(code_);//Could be zero.
      return ( m->expnt_sum() <= max_deg + time_expnt ) and ( time_expnt <= max_deg );
    }

 private:

  std::string  rep_; //String representation of time variable.
  unsigned int code_;//Table code of time variable.

  TotalDegree_TimeSeparate() : rep_("t"), code_(0) 
    {
      Variables::getTable().register_Observer(this);
    }
  TotalDegree_TimeSeparate(const TotalDegree_TimeSeparate&);
  ~TotalDegree_TimeSeparate()
    {
      Variables::getTable().remove_Observer(this);
    }

  TotalDegree_TimeSeparate& operator = (const TotalDegree_TimeSeparate&);
};

#endif

/*

  End of File: degree.h

*/
}
