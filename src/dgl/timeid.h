/*

 File: timeid.h, 2004/08/19

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

#ifndef TIMEID_H_INCLUDED
#define TIMEID_H_INCLUDED

#include "vartable.h"
namespace riot{
class TimeIdentifier
{
 private:

  /*
    To make it a Singleton.
  */
  
  TimeIdentifier() : identifier_(""), code_(0) {}
  TimeIdentifier(const TimeIdentifier&);

  TimeIdentifier& operator = (const TimeIdentifier&);
  ~TimeIdentifier() {}

 public:

  static TimeIdentifier& getObject()
    {
      static TimeIdentifier singleton;
      return singleton;
    }

  unsigned int setID(const std::string& s) 
    { 
      identifier_ = s;
      code_       = Variables::getTable().insert(s);

      return code_;
    }

  unsigned int getID() const { return code_; }
  unsigned int changeID(const std::string& s) { return setID(s); } //In future applications maybe at this point other
                                                                   //objects can be informed about the change.
 private:

  /*
    Data.
  */

  std::string  identifier_;
  unsigned int       code_;
};
}
#endif

/*

End of File: timeid.h

*/
