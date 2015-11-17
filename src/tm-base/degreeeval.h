/*

 File: degreeeval.h, 2004/08/26

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

#ifndef DEGREEEVAL_H_INCLUDED
#define DEGREEEVAL_H_INCLUDED

#include "observer.h"
#include "actions.h"
#include "vartable.h"
#include "monom.h"

namespace riot 
{
  

/*
  Monomgrad-Tester sind Objekte, die den Grad eines Monoms berechnen und pruefen 
  ob dieser kleiner dem angegebenen maximalen Grad ist.
  Objekte dieses Typs stehen mit der Variablentabelle in Verbindung, da unter 
  Umstaenden bei der Berechung bzw. Testung des Grades die ein oder andere Variable
  ausgezeichnet ist (d.h. eine Sonderbehandlung erfährt).

  Implementation with Barton and Nackman Trick is not possible!
*/

class DegreeBase //Gibt die Schnittstelle vor.
{
 public:

  virtual ~DegreeBase() {}

  //Delegate to leaf.
  virtual bool check(const Monom*, unsigned int) const = 0;
};

/*
  Objekte als Singleton.
*/

class Total_Degree : public DegreeBase //Behandelt alle Variablen gleich und berechnet den Gesamtgrad.
{
 public:

  static Total_Degree& getObject()
    {
      static Total_Degree singleton;
      return singleton;
    }

  void update(const Subject*, const Action*) {}

  bool check(const Monom* m, unsigned int max_deg) const
    {
      return m->expnt_sum() <= max_deg;
    }

 private:

  /*
    To make it a Singleton.
  */

  Total_Degree() {}
  Total_Degree(const Total_Degree&);

  ~Total_Degree() {}

  Total_Degree& operator = (const Total_Degree&);
};

}


#endif

/*

End of File: degreeeval.h

*/
