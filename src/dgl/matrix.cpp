
/*

 File: matrix.cpp, 2006/03/10

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

#include "matrix.h"
#include <cmath>

const int
  NoError   = 0,     // No error occurred
  NotSquare = 1,     // Matrix to be inverted is not square
  Singular  = 2;     // Matrix to be inverted is probably singular

/*
  Computes the determinant if an (n x n)-Matrix.
*/

double det(const Matrix& A, int& Err)
{
  const double
    Tiny = 1E-200; // A determinant less than 'Tiny' is handled like zero

  double  Max, Temp;                       // Help variables
  int     n = A.size(1);                   // Length of the rows of 'A'
  int     m = A.size(2);                   // Length of the columns of 'A'
  int     i, j, k, l, kk;                  // For loops

  Err = NoError;                           // Initialization

  if (n != m)                              // Error: 'A' is not square
    { 
      Err = NotSquare; 
      return 0.0; 
    }

  if (n == 1)
    {
      return A(0,0);
    }

  if (n == 2) // Special case: (2,2)-matrix
    {                            
      double d = A(0,0)*A(1,1)-A(1,0)*A(0,1); // Compute the determinant of 'A'

      if (std::abs(d) < Tiny)                 // Error: 'A' is probably singular
	{ 
	  Err = Singular; 
	  return 0.0; 
	}
      else return d;
    }

  // Start usual case: Dimension of 'A' > 2
  //---------------------------------------
  Matrix R = A;                // Use a copy of 'A'
  double v[n], x[n];           // Help vectors
  int p = 0;

  // Start LU factorization
  //-----------------------
  i = 1;
  while ( (Err == NoError) && (i <= n) ) {
    // Compute the numerators of those elements of the i-th column
    // of L which are not updated so far and store them in 'v'.
    //------------------------------------------------------------
    for (k = i; k <= n; k++) {   
      v[k-1] = R(k-1,i-1);
      for (j = 1; j < i; j++) v[k-1] -= R(k-1,j-1)*R(j-1,i-1);
    }

    // Look for the column pivot
    //--------------------------
    j = i;  Max = std::abs(v[i-1]);
    for (k = i+1; k <= n; k++)
      if ( (Temp = std::abs(v[k-1])) > Max ) { j = k; Max = Temp; }

    // Swap rows of 'R' and 'v', store number of changes in 'p'
    //-------------------------------------------------------
    if (j != i) {
      for(int ii = 1; ii <= m; ii++)
	{
	  x[ii-1] = R(i-1,ii-1); R(i-1,ii-1) = R(j-1,ii-1);  R(j-1,ii-1) = x[ii-1];
	}
      ++p;
      Temp = v[i-1];  v[i-1] = v[j-1];  v[j-1] = Temp;
    }

    if (Max < Tiny) // Pivot element < 'Tiny', inversion
      {             // failed, matrix 'A' assumed to be singular
	Err = Singular; 
	return 0.0; 
      }     

    Temp = v[i-1];
    R(i-1,i-1) = Temp;    // Update the diagonal element of U

    // Update U's row and L's column elements
    //---------------------------------------
    for (k = i+1; k <= n; k++) {
      for (j = 1; j < i; j++) R(i-1,k-1) -= R(i-1,j-1)*R(j-1,k-1);
      R(k-1,i-1) = v[k-1] / Temp;
    }
    i++;
  } // while

  // Now 'R' is overwritten with the subdiagonal elements of L in its
  // lower left triangle and with the elements of U in its diagonal and
  // its upper right triangle. The determinant now is the product of the
  // diagonal elements of 'R'.
  //-------------------------------------------------------------------
  double d = R(0,0);
  for(i = 1; i < n; i++) d *= R(i,i);

  if( p & 1 ) return -d;
  else        return  d;
}

Matrix abs(const Matrix& A)
{
  unsigned int r = A.size(1), c = A.size(2);

  Matrix B(r,c);

  for(unsigned int i = 0; i < r; i++)
    for(unsigned int j = 0; j < c; j++)
      if( A(i,j) < 0 ) B(i,j) = -A(i,j);
      else             B(i,j) =  A(i,j);

  return B;
}

double Matrix::max() const
{
  if( row_*col_ > 0 )
    {
      double m = _M_[0];
      for(unsigned int i = 1; i < row_*col_; i++) m = (m < _M_[i])?_M_[i]:m;
      return m;
    }
  else return 0.0;
}

double Matrix::min() const
{
  if( row_*col_ > 0 )
    {
      double m = _M_[0];
      for(unsigned int i = 1; i < row_*col_; i++) m = (m > _M_[i])?_M_[i]:m;
      return m;
    }
  else return 0.0;
}

/*

End of File: matrix.cpp

*/
