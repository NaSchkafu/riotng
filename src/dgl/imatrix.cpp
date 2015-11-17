/*

 File: imatrix.cpp, 2005/01/28

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

#include "imatrix.h"
#include "linsys.h"

/*
  The following function computes an enclosure from the inverse
  of a given real matrix using the function 'LinSolve'.
*/

void imatinv(const Matrix& A, IMatrix& B, int& err_code)
{
  /*
    The size of 'A' and 'B' is not checked
  */

  int n = A.size(1);

  IVector col(n);
  Vector  rhs(n);
  
  for(int i = 0; i < n; i++)
    {
      rhs[i] = 1.0;
      err_code = 0;
      
      LinSolve( A, rhs, col, err_code );
      
      if( err_code != 0 )
	{
	  cout << LinSolveErrMsg(err_code) << endl;
	  return;
	}
      
      //Set column.
      for(int j = 0; j < n; j++) B(j,i) = col[j];
      
      rhs[i] = 0.0;
    }
}

/*

End of File: imatrix.cpp

*/
