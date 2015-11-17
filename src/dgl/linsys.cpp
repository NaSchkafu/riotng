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
// File: linsys (implementation)
// Purpose: Computation of a verified solution of a square linear system of
//    equations A*x = b with full real matrix A and real right-hand side b.
// Method: Transformation of A*x = b to fixed-point form and applying an
//    interval residual iteration.
// Global functions:
//    LinSolve()      : to get a verified enclosure of the solution (two
//                      versions)
//    LinSolveErrMsg(): to get an error message text
//----------------------------------------------------------------------------
// Portet to Filib++ on 03/07/16
//----------------------------------------------------------------------------
#include <string.h>         // String handling
#include <cmath>
#include <cstdio>
#include "imatrix.h"        // Interval matrix
#include "matinv.h"         // Matrix inversion
#include "linsys.h"

using namespace std;

int sign(const double& x)
{
  if( x > 0.0 ) return 1;
  else if( x < 0.0 ) return -1;
  else return 0;
}

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError      = 0,   // No error occurred
  NotSquare    = 1,   // System to be solved is not square
  DimensionErr = 2,   // Dimensions of A and b are not compatible
  InvFailed    = 3,   // System is probably singular
  VerivFailed  = 4;   // Verification failed, system is probably
                      // ill-conditioned

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* LinSolveErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    switch (Err) {
      case NotSquare:
        strcpy(Hlp,"System to be solved is not square"); break;
      case DimensionErr:
        strcpy(Hlp,"Dimensions of A and b are not compatible"); break;
      case InvFailed:
        strcpy(Hlp,"System is probably singular"); break;
      case VerivFailed:
        strcpy(Hlp,"Verification failed, system is probably ill-conditioned");
        break;
      default:
        strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // LinSolveErrMsg

//----------------------------------------------------------------------------
// Computes an upper bound for the maximum norm of a real matrix 'M'.
//----------------------------------------------------------------------------
static double MaxNorm(const Matrix& M)
{
  int          i, j;
  double       Max, Tmp;

  Max = 0.0;
  for(i = M.Lb(1); i <= M.Ub(1); i++) {
    Tmp = 0.0;
    for(j = M.Lb(2); j <= M.Ub(2); j++) 
#ifdef FILIB_VERSION
      Tmp = filib::fp_traits<double>::upward_plus( Tmp, std::abs(M(i,j)), true );
#endif
#ifdef CXSC_VERSION
      Tmp = cxsc::_double( cxsc::addup( Tmp, std::abs(M(i,j)) ) );
#endif

    if (Tmp > Max) Max = Tmp;
  }
  return Max;
} // MaxNorm

//----------------------------------------------------------------------------
// The vectors x and y are successive approximations for the solution of a
// linear system of equations computed by iterative refinement. If a compo-
// nent of y is diminished by more than 'Factor', it is a good candidate for
// a zero entry. Thus, it is set to zero.
//----------------------------------------------------------------------------
static void CheckForZeros (double* x, double* y, int n)
{
  const double Factor = 1E+5;
  int          i;

  for(i = 0; i <= n-1; i++)
    if( std::abs(y[i])*Factor < std::abs(x[i]) ) y[i] = 0.0;
}

//----------------------------------------------------------------------------
// The vectors x and y are successive iterates. The function returns true if
// the relative error of all components x_i and y_i is <= 10^(-12), i.e. y_i
// has about 12 correct decimals. If x_i or y_i vanishes, the relative error
// is implicitly set to zero.
//----------------------------------------------------------------------------
static int Accurate (double* x, double* y, int n)
{
  const double Delta = 1E-12;   // Relative error bound
  int          i, ok;
  double       abs_yi;

  i = 0;
  do {
    if(sign(x[i])*sign(y[i]) <= 0)       // Relative error set to 0
      ok = 1;
    else {
     abs_yi = std::abs(y[i]);            // Relative error > Delta ?
     ok = (std::abs(y[i] - x[i]) <= Delta * abs_yi );
    }
    i++;
  } while (ok && (i < n));
  return ok;
} // Accurate

void Blow(const IVector& x, IVector& y, int n, double eps)    // Epsilon deflation
{                                                         //------------------
  for(int i = 0; i < n; i++) y[i] = blow(x[i],eps);
}

int in (const IVector& x, const IVector& y, int n)   // Inner inclusion for two ivectors
{                                          //---------------------------------
  int   i = 0, ok = 1;

  while(ok && i < n) { ok = interior(x[i],y[i]); i++; }
  return ok;
}

//----------------------------------------------------------------------------
// This function 'VerificationStep()' performs the iteration
// [y] = Blow([x],Eps), [x] = [z] + [C]*[y] for k = 1,2,... until the new
// iterate [x] lies in the interior of [y] or until the maximum number of
// iterations is exceeded. The flag 'IsVerified' is set if an inclusion in
// the interior could be established.
//----------------------------------------------------------------------------
static void VerificationStep ( IVector& xx, const IVector& zz, int n,
                               IMatrix& C,  int& IsVerified )
{
  const int    MaxIter = 3;          // Maximal number of iteration steps
  const double Epsilon = 1000.0;     // Factor for the epsilon inflation

  int     p;
  IVector yy(n);

  for(int i = 0; i < n; i++) xx[i] = zz[i]; // Initialize:  [x] = [z]
  p = 0;
  do {
    Blow(xx,yy,n,Epsilon);             // Epsilon inflation
    for(int i = 0; i < n; i++)         // New iterate: [x] = [z]+[C]*[y]
      {
	xx[i] = zz[i];
	for(int j = 0; j < n; j++) xx[i] += C(i,j)*yy[j];
      }
    IsVerified = in(xx,yy,n);          // Inclusion in the interior?
    p++;
  } while (!IsVerified && (p < MaxIter));
}

int Zero (const IVector& x, int n)                               // Check for zero vector
{                            
  int i, ok;
  for (i = 0, ok = 1; i < n && ok; i++) ok = (x[i] == 0.0);
  return ok;
}

//----------------------------------------------------------------------------
// Purpose: The function 'LinSolveMain()' computes a verified solution of a
//    square linear system of equations A*x=b.
// Parameters:
//    In : 'A'          : matrix of the system, passed as reference
//                        parameter to save time for copying it
//         'b'          : right-hand side of the system, passed as
//                        reference parameter to save time for copying it
//         'ComputeCond': flag signalling whether a condition number
//                        estimate is to be computed
//    Out: 'xx'         : enclosure of the unique solution, resized to
//                        standard index range with lower index bound 1
//         'Cond'       : condition number estimate
//         'Err'        : error code
// Description: An approximate inverse 'R' of 'A' is computed by calling
//   function 'MatInv()'. Then an approximate solution 'x' is computed
//   applying a conventional real residual iteration. For the final verifica-
//   tion, an interval residual iteration is performed. An enclosure of the
//   unique solution is returned in the interval vector 'xx'. The function
//   also returns a condition number estimate 'Cond' if the flag 'ComputeCond'
//   is set. 'Cond' is initialized by -1. A negative value for 'Cond' signals
//   that an estimate could not be computed.
//----------------------------------------------------------------------------
static void LinSolveMain (const Matrix&   A,
                          const Vector&   b,
                          IVector&       xx,
			  double&      Cond,
                          int   ComputeCond,
                          int&          Err )
{
  const int
    MaxResCorr = 10;    // Maximal number of real residual corrections

  int     n = A.size(1);                  // Length of the columns of 'A'
  int     m = A.size(2);                  // Length of the rows of 'A'
  int     IsVerified, i, j, k;

  Matrix R;                             // To store the inverse of 'A'

  Cond = -1.0;                          // Initial value for the condition

  if (n != m)                           // Error: 'A' is not square
    { Err = NotSquare; return; }

  if (n != b.size())                    // Error: Dimensions of 'A' and 'b'
    { Err = DimensionErr; return; }     //        are not compatible

  MatInv(A,R,Err);
  if (Err != NoError)                   // Error: Inversion failed
    { Err = InvFailed; return; }

  // Start algorithm
  //----------------
  double  x[n], y[n], d[n];       // Allocate dynamic arrays
  IVector yy(n), zz(n);
  IMatrix C(n,n);

  if(ComputeCond)                         // Compute condition number
#ifdef FILIB_VERSION
    Cond = filib::fp_traits<double>::upward_multiplies( MaxNorm(A), MaxNorm(R), true );
#endif
#ifdef CXSC_VERSION
    Cond = cxsc::_double( cxsc::multup( MaxNorm(A), MaxNorm(R) ) );
#endif

  // Normalize index range of 'A', 'R' and 'b' to standard range 1..n
  // and resize 'xx' to length 'n'. Save lower bounds of 'b' and 'A'
  // to restore them before leaving 'LinSolveMain()'.
  //-----------------------------------------------------------------

  for(int ii = 0; ii < n; ii++) // Real residual iteration
    {
      x[ii] = 0.0;
      for(int jj = 0; jj < n; jj++) x[ii] += R(ii,jj)*b[jj];
    }
  
  k = 0;
  do {
    for(int ii = 0; ii < n; ii++) y[ii] = x[ii];
    for(i = 1; i <= n; i++) {   // Compute: d = #*(b-A*y)
      d[i-1] = b[i-1];          //-----------------------
      for(int jj = 0; jj < n; jj++) d[i-1] -= A(i-1,jj)*y[jj];
    }
    for(i = 1; i <= n; i++) {   // Compute: x = #*(y+R*d)
      x[i-1] = y[i-1];          //-----------------------
      for(int jj = 0; jj < n; jj++) x[i-1] += R(i-1,jj)*d[jj];
    }
    CheckForZeros(y,x,n);
    k++;
  } while (!Accurate(y,x,n) && (k < MaxResCorr));

  // Prepare verification step, i.e. compute enclosures [C]
  // and [z] of C = (I - R*A) and z = R*(b - A*x).
  //-------------------------------------------------------
  for(i = 1; i <= n; i++)
    for(j = 1; j <= n; j++) {            // Compute [C] = ##(I-R*A)
      C(i-1,j-1) = (i == j) ? 1.0 : 0.0;
      for(int ikj = 0; ikj < n; ikj++) C(i-1,j-1) -= Interval(R(i-1,ikj)) * A(ikj,j-1);
    }
  for(i = 1; i <= n; i++) {              // Compute: d = #*(b-A*x)
    d[i-1] = b[i-1];
    for(int jj = 0; jj < n; jj++) d[i-1] -= A(i-1,jj)*x[jj];
  }
  for(i = 1; i <= n; i++) {              // Compute: [y] = ##(b-A*x-d)
    yy[i-1] = b[i-1];
    for(int jj = 0; jj < n; jj++) yy[i-1] -= Interval(A(i-1,jj)) * x[jj];
    yy[i-1] -= d[i-1];
  }

  // Now b-A*x lies in d+[y]. Thus R*(b-A*x) lies in [z] = ##(R*d+R*[y])
  //--------------------------------------------------------------------
  for (i = 1; i <= n; i++) {
    zz[i-1] = 0.0;
    for(int jj = 0; jj < n; jj++) zz[i-1] += Interval(R(i-1,jj)) * d[jj];
    for(int jj = 0; jj < n; jj++) zz[i-1] += R(i-1,jj) * yy[jj];
  }

  // If R*(b-A*x) = 0 then x is an exact solution else try to find a
  // verified enclosure [x] for the absolute error of x.
  //----------------------------------------------------------------
  if(Zero(zz,n))
    for(int ii = 0; ii < n; ii++) xx[i] = x[i];
  else {
    VerificationStep(xx,zz,n,C,IsVerified);   // Attempt to compute [x]
    if ( !IsVerified )
      Err = VerivFailed;
    else
      for(int ii = 0; ii < n; ii++) xx[ii] += x[ii];  // The exact solution lies x+[x]
  }
} //LinSolveMain

//----------------------------------------------------------------------------
// Purpose: The function 'LinSolve()' computes a verified solution of a
//    square linear system of equations A*x=b without returning a condition
//    number estimate.
// Parameters:
//    In : 'A'  : matrix of the system, passed as reference
//                parameter to save time for copying it
//         'b'  : right-hand side of the system, passed as
//                reference parameter to save time for copying it
//    Out: 'xx' : enclosure of the unique solution
//         'Err': error code
// Description: Calls 'LinSolveMain()' for solving the linear system with the
//    flag 'ComputeCond' not set.
//----------------------------------------------------------------------------
void LinSolve(const Matrix& A, const Vector& b, IVector& xx, int& Err )
{
  double DummyCond;         // Dummy parameter for call of 'LinSolveMain()'

  LinSolveMain(A,b,xx,DummyCond,0,Err);
}

//----------------------------------------------------------------------------
// Purpose: The function 'LinSolve()' computes a verified solution of a
//    square linear system of equations A*x=b and returns a condition
//    number estimate.
// Parameters:
//    In : 'A'   : matrix of the system, passed as reference
//                 parameter to save time for copying it
//         'b'   : right-hand side of the system, passed as
//                 reference parameter to save time for copying it
//    Out: 'xx'  : enclosure of the unique solution
//         'Cond': condition number estimate
//         'Err' : error code
// Description: Calls 'LinSolveMain()' for solving the linear system with the
//    flag 'ComputeCond' set.
//----------------------------------------------------------------------------
void LinSolve (const Matrix& A, const Vector& b, IVector& xx, double& Cond, int& Err)
{
  LinSolveMain(A,b,xx,Cond,1,Err);
}
