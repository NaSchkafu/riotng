/*

  File: actions.h, 2004/08/18

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

#ifndef ACTIONS_H_INCLUDED
#define ACTIONS_H_INCLUDED

#include "observer.h"
#include <string>
namespace riot {

/*
  Some actions that may occur in conjunction with the table of variables.
*/

  class IdentifierChanged : public Action
  {
  private:

    /*
      To make it a Singleton.
    */

    IdentifierChanged() {}
    IdentifierChanged(const IdentifierChanged&);
    ~IdentifierChanged() {}

    IdentifierChanged& operator = (const IdentifierChanged&);

  public:

    /*
      Data.
    */

    std::string old_identifier_, new_identifier_;
    unsigned int code_;

    /*
      Operations.
    */

    static IdentifierChanged& getObject()
    {
      static IdentifierChanged singleton;
      return singleton;
    }

    void set(const std::string& old_id, const std::string& new_id, unsigned int c)
    {
      old_identifier_ = old_id;
      new_identifier_ = new_id;
      code_ = c;
    }
  };

  class NewVariable : public Action
  {
  private:

    /*
      To make it a Singleton.
    */

    NewVariable() {}
    NewVariable(const NewVariable&);
    ~NewVariable() {}

    NewVariable& operator = (const NewVariable&);

  public:

    /*
      Data.
    */

    std::string  identifier_;
    unsigned int code_;

    /*
      Operations.
    */

    static NewVariable& getObject()
    {
      static NewVariable singleton;
      return singleton;
    }

    void set(const std::string& id, unsigned int c)
    {
      identifier_ = id;
      code_       = c;
    }
  };
}

#endif

/*

  End of File: actions.h

*/
