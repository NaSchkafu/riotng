/*

  File: observer.h, 2004/09/01

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

#ifndef OBSERVER_H_INCLUDED
#define OBSERVER_H_INCLUDED

#include <vector>
#include <algorithm>
namespace riot
{

  class Action
  {
  public:

    virtual ~Action() {}

    virtual unsigned int id() const { return 0; }
  };

  class Subject;

  class Observer
  {
  public:

    virtual ~Observer() {}

    virtual void update(const Subject*,const Action*) = 0;
  };

  class Subject
  {
  public:

    void register_Observer(Observer *v)
    {
      std::vector<Observer*>::const_iterator
        pos = std::find( tab_.begin(), tab_.end(), v );

      if( pos == tab_.end() ) tab_.push_back( v );

      Action confirmation;
      v->update(this, &confirmation);
    }
    void remove_Observer(Observer *v)
    {
      std::vector<Observer*>::iterator
        pos = std::find( tab_.begin(), tab_.end(), v );

      if( pos != tab_.end() ) tab_.erase( pos );
    }
    void inform_Observers(Action* a) const
    {
      for(unsigned int i = 0; i < tab_.size(); i++) (tab_[i])->update(this, a);
    }

  protected:

  Subject() : tab_(0) {}

  private:

    std::vector<Observer*> tab_;
  };
}

#endif

/*

  End of File: observer.h

*/
