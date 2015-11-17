/*

 File: output.h, 2004/08/24

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

#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

#include <iostream>
#include <cstdio>
#include <string>

typedef FILE* FILE_ptr;

class Files
{
 private:

  unsigned int  size_;
  FILE_ptr     *file_;

 public:

  Files() : size_(0), file_((FILE_ptr*)0) {}
  Files(unsigned int n, const std::string s) : size_(n)
    {
      file_ = new FILE_ptr[n];
      char *c = new char[20];
      
      for(unsigned int i = 0; i < n; i++) 
	{
	  sprintf(c,"_%02u%s",i,".txt");
	  if( not (file_[i] = fopen( (s+c).c_str(), "wt" )) )
	    {
	      std::cout << "*** Can't create output files. ***" << std::endl;
	      std::exit(1);
	    }
	}
    }
  ~Files()
    {
      for(unsigned int i = 0; i < size_; i++) 
	{
	  fclose( file_[i] );
	}
      delete file_;
    }

  /*
    Access.
  */

  FILE_ptr operator [](unsigned int i) const { return file_[i]; }
};

#endif

/*

  End of File: output.h

*/
