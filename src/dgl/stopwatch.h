/*

 File: stopwatch.h, 2004/08/23

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

#ifndef _STOPWATCH_H_INCLUDED
#define _STOPWATCH_H_INCLUDED

#include <sys/times.h>
#include <unistd.h>

class Stopwatch
{
  struct tms start_cpu;
  struct tms end_cpu;

  double accu;

 public:
  
  Stopwatch() : accu(0.0) {}
  ~Stopwatch() {}

  void start() { times(&start_cpu); }
  void stop () 
    { 
      times(&end_cpu  ); 
      accu += elapsedTime();
    }
  double elapsedTime() const
    {
      return ((double)(end_cpu.tms_utime - start_cpu.tms_utime))/sysconf(_SC_CLK_TCK);
    }
  double accumulatedTime() const
    {
      return accu;
    }
  void clearAccu() { accu = 0.0; }
};

#endif


/*

  End of File: stopwatch.h

*/
