//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz
//
// This program/module is free software for non-commercial use. For details
// on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
// This program/module is distributed WITHOUT ANY WARRANTY. For details,
// see the "Disclaimer / Legal Matters" of the book (page iv).
//
//============================================================================
//----------------------------------------------------------------------------
// File: matinv (implementation)
// Purpose: Computation of an approximate inverse of a real square matrix.
// Method: LU decomposition applying Crout's algorithm.
// Global functions:
//   MatInv()      : matrix inversion
//   MatInvErrMsg(): to get an error message text
//----------------------------------------------------------------------------
// Portet to Filib++ on 03/07/16
//----------------------------------------------------------------------------
#include <string>       // String handling
#include <cmath>
#include <cstdio>
#include <cstring>
#include "matinv.h"

using namespace std;

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError      = 0,     // No error occurred
  NotSquare    = 1,     // Matrix to be inverted is not square
  Singular     = 2;     // Matrix to be inverted is probably singular

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* MatInvErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    switch (Err) {
    case NotSquare:
      sprintf(Hlp,"Matrix to be inverted is not square");           break;
    case Singular:
      sprintf(Hlp,"Inversion failed, matrix is probably singular"); break;
    default:
      strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // MatInvErrMsg

//----------------------------------------------------------------------------
// Purpose: The function 'MatInv()' may be used for the computation of an
//    approximate inverse of a real square matrix.
// Parameters:
//    In : 'A'  : matrix to be inverted
//    Out: 'R'  : approximate inverse
//         'Err': error code
// Description:
//    Inversion of a regular matrix A stored in 'A' using LU decomposition.
//    For LU decomposition, formally a permutation matrix P is determined so
//    that the product P*A may be decomposed into a lower-triangular matrix L
//    and an upper-triangular matrix U with P*A = L*U. The diagonal elements
//    of L are set to 1. Using Crout's algorithm, the elements of both matri-
//    ces L and U are stored by temporary overwriting a copy of the input
//    matrix 'A'. The permutation matrix P is not explicitly generated. The
//    indices of row interchanges are stored in the index vector 'p' instead.
//    The i-th element of P*b may be accessed indirectly using the p[i]-th
//    entry of 'b'. The k-th column of the inverse R of P*A is computed by
//    forward/backward substitution with the k-th unit vector e_k as the
//    right-hand side of the system: U*y = P*e_k ==> y, L*x = y ==> x. For
//    error codes, see above.
//----------------------------------------------------------------------------
void MatInv (Matrix A, Matrix& R, int& Err)
{
  const double
    Tiny = 1E-200;          // A divisor less than 'Tiny' is handled like zero

  double  Max, Temp;                       // Help variables
  int     p1 = A.Lb(1), q1 = A.Lb(2);      // Lower bounds of 'A'.
  int     pn = A.Ub(1), qm = A.Ub(2);      // Upper bounds of 'A'.
  int     n = pn-p1+1;                     // Length of the rows of 'A'
  int     m = qm-q1+1;                     // Length of the columns of 'A'
  int     i, j, k, l, kk;                  // For loops

  Err = NoError;                           // Initialization

  if (n != m)                              // Error: 'A' is not square
    { Err = NotSquare; return; }

  // Normalize index range of 'A' to standard range (1..n,1..n)
  //-----------------------------------------------------------
  R = A;                                   // Resize 'R' to same shape as 'A'

  if (n == 1)
    {
      R(p1,q1) = 1.0 / A(0,0);
      return;
    }

  if (n == 2) 
    {                            // Special case: (2,2)-matrix
      double det = A(0,0)*A(1,1)-A(1,0)*A(0,1); // Compute the determinant of 'A'

      if (std::abs(det) < Tiny)                  // Error: 'A' is probably singular
	{ Err = Singular; return; }

      R(p1,q1) =  A(1,1) / det;  R(p1,qm) = -A(0,1) / det;
      R(pn,q1) = -A(1,0) / det;  R(pn,qm) =  A(0,0) / det;
      return;
  }

  // Start usual case: Dimension of 'A' > 2
  //---------------------------------------
  double v[n], x[n];           // Help vectors
  int*   p = new int [n+1];    // Dynamic array of type integer
                               // Note: p[0] not used !

  for (i = 1; i <= n; i++)  p[i] = i; // Initializing permutation vector

  // Start LU factorization
  //-----------------------
  i = 1;
  while ( (Err == NoError) && (i <= n) ) {
    // Compute the numerators of those elements of the i-th column
    // of L which are not updated so far and store them in 'v'.
    //------------------------------------------------------------
    for (k = i; k <= n; k++) {   
      v[k-1] = A(k-1,i-1);
      for (j = 1; j < i; j++) v[k-1] -= A(k-1,j-1)*A(j-1,i-1);
    }

    // Look for the column pivot
    //--------------------------
    j = i;  Max = std::abs(v[i-1]);
    for (k = i+1; k <= n; k++)
      if ( (Temp = std::abs(v[k-1])) > Max ) { j = k; Max = Temp; }

    // Swap rows of 'A' and 'v', store the permutation in 'p'
    //-------------------------------------------------------
    if (j != i) {
      for(int ii = 1; ii <= m; ii++)
	{
	  x[ii-1] = A(i-1,ii-1); A(i-1,ii-1) = A(j-1,ii-1);  A(j-1,ii-1) = x[ii-1];
	}
      p[i] ^= p[j] ^= p[i] ^= p[j];
      Temp = v[i-1];  v[i-1] = v[j-1];  v[j-1] = Temp;
    }

    if (Max < Tiny)                   // Pivot element < 'Tiny', inversion
      { Err = Singular; return; }     // failed, matrix 'A' assumed to be
                                      // singular
    Temp = v[i-1];
    A(i-1,i-1) = Temp;    // Update the diagonal element of U

    // Update U's row and L's column elements
    //---------------------------------------
    for (k = i+1; k <= n; k++) {
      for (j = 1; j < i; j++) A(i-1,k-1) -= A(i-1,j-1)*A(j-1,k-1);
      A(k-1,i-1) = v[k-1] / Temp;
    }
    i++;
  } // while

  // Now 'A' is overwritten with the subdiagonal elements of L in its
  // lower left triangle and with the elements of U in its diagonal and
  // its upper right triangle. The elements of the inverse matrix are
  // computed column by column using forward/backward substitution.
  //-------------------------------------------------------------------
  for (k = 1; k <= n; k++) {
    // Forward substitution: L*x = P*e_k, where e_k is the k-th unit
    // vector. Note: If P*e_k has m leading zeros, this results in
    // x_i = 0 for 1,..,l-1 and x_l = 1. Thus, forward substitution
    // starts at index l+1.
    //--------------------------------------------------------------
    l = 1;
    while (p[l] != k) { x[l-1] = 0.0; l++; }
    x[l-1] = 1.0;
    for (i = l+1; i <= n; i++) {
      x[i-1] = 0.0;
      for (j = l; j < i; j++) x[i-1] -= A(i-1,j-1)*x[j-1];
    }

    // Backward substitution: U * x = x, where the right-hand side is
    // the result of the forward substitution. It will be overwritten
    // by the solution of the system, the k-th column of the inverse
    // matrix.
    //---------------------------------------------------------------
    kk = q1+k-1;                 // Index offset for the k-th column of R

    for (i = n; i >= 1; i--) {
      for (j = i+1; j <= n; j++) x[i-1] -= A(i-1,j-1)*x[j-1];
      x[i-1] /= A(i-1,i-1);
      R(p1+i-1,kk) = x[i-1];      // Remember index offset !
    }
  } // for (k = 1,...

  delete [] p;   // Free dynamically allocated memory
} // MatInv
