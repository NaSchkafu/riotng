/*
  
 File: dglalg.h, 2005/02/21

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

#ifndef DGLALG_H_INCLUDED
#define DGLALG_H_INCLUDED

#include "data.h"
#include "vvtaylormodel.h"
#include "ivector.h"
#include "output.h"

namespace riot{

class DGLSolver
{
 public:

  static DGLSolver& Algorithm()
    {
      static DGLSolver singleton;
      return singleton;
    }

  void prepareAlgorithm(Data*);
  void run(const Files*);

 private:

  DGLSolver() : AlgorithmPrepared(false) {}
  DGLSolver(const DGLSolver&);

  DGLSolver& operator = (const DGLSolver&);
  ~DGLSolver() {}

  //Functions.
  bool terminate(const IVector&,const IVector&); //Termination criteria in iterative refinement of the solution.
  void inflate(Interval&); //Epsilon Inflation.

  //Data.
  bool       AlgorithmPrepared;
  Data              *Data_Ptr_;
  unsigned int       TimeCode_;

  VVTaylorModel    InitValsTM_;
  VVTaylorModel InitValsDGLTM_;
  VVTaylorModel      Solution_;
  VVTaylorModel        currTM_;
  TaylorModel          TimeTM_;

  IVector              currIV_;

  double               h_curr_;
  double      t_curr_, t_next_;
  Interval              t_end_;
};
}
#endif

/*

  End of File: dglalg.h

*/
