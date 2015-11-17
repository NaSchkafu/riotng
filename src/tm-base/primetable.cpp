/*

  File: primetable.h, 2004/08/18

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

#include "primetable.h"
#include <algorithm>
namespace riot
{

  PrimeNumbers::PrimeNumbers(unsigned long n)
  {
    /*
      Die Primzahlen werden mit dem Algorithmus "Sieb des Eratosthenes" berechnet.
    */

    // Lege ein Array mit n+1 Werten an und setzte alle Werte auf "true"
    // (wir nehmen zunaechst an, dass alle Zahlen Primzahlen sind)

    bool *sieb = new bool[n + 1];

    sieb[0] = false; sieb[1] = false; // 0 und 1 sind keine Primzahlen
    for(unsigned long int i = 2; i < (n + 1) ; i++) sieb[i] = true;

    // Für jede Zahl von 2 bis n, streiche die Vielfachen dieser Zahl

    for(unsigned long int faktor = 2; faktor < (n + 1); faktor++)
    {
      // Wenn der aktuelle faktor bereits ausgesiebt wurde,
      // sind seine Vielfachen ebenfalls bereits
      // gesiebt. Andernfalls sieben wir seine Vielfachen jetzt.

      if(sieb[faktor])
      {
        unsigned int aktuelles_element = faktor;

        // Erhöhe den Wert von aktuelles_element um faktor und siebe
        // die Zahl aus, weil sie ein Vielfaches von faktor ist.

        while((aktuelles_element += faktor) < (n + 1)) sieb[aktuelles_element] = false;
      }
    }

    //Eintragen der Primzahlen in die Tabelle

    for(unsigned long int i = 0; i < (n + 1); i++)
      if(sieb[i]) table_.push_back(i);

    delete[] sieb;

    locked_prime_ = 0;
  }

  unsigned long PrimeNumbers::find_SmallestPrimeGreaterThan(const double& k) const
  {
    std::vector<unsigned long>::const_iterator
      pos = std::upper_bound(table_.begin(),table_.end(),k); // O(log(n))

    if( pos == table_.end() ) return *(--pos);
    else
    {
      if( *pos == locked_prime_ ) return *(++pos); //Return the next greater prime.
      else                        return *pos;
    }
  }

  unsigned long PrimeNumbers::find_GreatestPrimeSmallerThan(const double& k) const
  {
    std::vector<unsigned long>::const_iterator
      pos = std::lower_bound(table_.begin(),table_.end(),k); // O(log(n))

    if( pos == table_.end() ) return 1;
    else
    {
      if( *pos == locked_prime_ ) return *(--pos); //Return the next lower prime.
      else                        return *pos;
    }
  }

/*

  End of File: primetable.cpp

*/
}
