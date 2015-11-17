// -*-c++-*-
/*

  File: vartable.h, 2004/08/18

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

#ifndef VARTABLE_H_INCLUDED
#define VARTABLE_H_INCLUDED

#include "observer.h"
#include "actions.h"

#include <vector>
#include <string>
namespace riot
{

/*
  The Variables is implementet as Singleton.
  The storage of the data is useful only for less than 100 entries.
*/

  class Variables : public Subject
  {
  private:

  Variables() : Subject(), num_(0) { Tab_.reserve(100); }
    Variables(const Variables&);
    ~Variables() {}

    Variables& operator = (const Variables&);

  public:

    static Variables& getTable()
    {
      static Variables singleton;
      return singleton;
    }

    unsigned int insert(const std::string& s)
    {
      if( s != "" )
      {
        std::vector<std::string>::const_iterator
          pos = find( Tab_.begin(), Tab_.end(), s ); // O(n)

        if( pos == Tab_.end() ) //No entry with string s.
        {
          Tab_.push_back(s);
          NewVariable::getObject().set( s, ++num_ ); //Define action.
          Subject::inform_Observers(&NewVariable::getObject());//Inform observers about this action.
          return num_;
        }
        else return (pos - Tab_.begin() + 1);
      }
      else return 0;
    }
    void change_identifier(const std::string& old_id, const std::string& new_id)
    {
      if( old_id != "" && new_id != "" )
      {
        std::vector<std::string>::iterator
          pos = find( Tab_.begin(), Tab_.end(), old_id ); // O(n)

        if( pos != Tab_.end() ) //Entry with string 'old_id' found.
        {
          (*pos) = new_id;                                                    //Change identifier.
          IdentifierChanged::getObject().set( old_id, new_id, (pos - Tab_.begin() + 1) ); //Define action.
          Subject::inform_Observers(&IdentifierChanged::getObject());      //Inform observers about this action.
        }
      }
    }
    unsigned int look_up(const std::string& s) const
    {
      if( s != "" )
      {
        std::vector<std::string>::const_iterator
          pos = find( Tab_.begin(), Tab_.end(), s ); // O(n)

        if( pos == Tab_.end() ) //no entry with string s
          return 0;
        else return (pos - Tab_.begin() + 1);
      }
      else return 0;
    }
    const std::string& look_up(unsigned int i) const
    {
      static std::string empty_string("");

      if( 0 < i && i <= num_ ) return Tab_[ i-1 ];
      else                     return empty_string;
    }

    unsigned int size() const { return num_; }

  private:

    unsigned int num_;
    std::vector<std::string> Tab_;
  };
}

#endif

/*

  End of File: vartable.h

*/
