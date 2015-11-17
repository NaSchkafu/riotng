/*

 File: linsys.h, 2003/09/25

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

#ifndef LINSYS_H_INCLUDED
#define LINSYS_H_INCLUDED

#include "matrix.h"
#include "dvector.h"
#include "ivector.h"
#include "adaptintval.h"

using Adapt::Interval;
using namespace std;

extern char* LinSolveErrMsg ( int );
extern void  LinSolve (const Matrix&, const Vector&, IVector&, int& );
extern void  LinSolve (const Matrix&, const Vector&, IVector&, double&, int& );

#endif

/*

End of File: linsys.h

*/
