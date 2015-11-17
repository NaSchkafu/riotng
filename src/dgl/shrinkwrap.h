/*

 File: shrinkwrap.h, 2005/02/21

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

#ifndef SHRINKWRAP_H_INCLUDED
#define SHRINKWRAP_H_INCLUDED

#include "vvtaylormodel.h"
#include "monom.h"
#include "data.h"
#include "output.h"
namespace riot{
/*

Implementation of the Shrink Wrapping algorithm.

*/

class ShrinkWrapping
{
 public:

  static ShrinkWrapping& Algorithm()
    {
      static ShrinkWrapping singleton;
      return singleton;
    }

  void prepareAlgorithm(const Data*);
  void run(VVTaylorModel*,const Files*);

 private:

  ShrinkWrapping() : MonomList_(0), AlgorithmPrepared(false) {}
  ShrinkWrapping(const ShrinkWrapping&);

  ShrinkWrapping& operator = (const ShrinkWrapping&);
  ~ShrinkWrapping() {}

  //Functions.
  void splitVVTM(VVTaylorModel*);
  void wrap     (VVTaylorModel*,int&);

  //Data.
  bool             AlgorithmPrepared;
  VVTaylorModel          InitValsTM_;
  Monom                  *MonomList_;
  unsigned int         *VarCodeList_;

  unsigned int            Dimension_;

  Vector                   Constant_;
  Matrix                 Linear_Mat_;
  Matrix             Linear_Mat_Inv_;
  IMatrix        Linear_Mat_Inv_Inv_;
  VVTaylorModel           Nonlinear_;
  IVector                 Remainder_;

  double                     d,s,t,q;
};
}
#endif

/*

  End of File: shrinkwrap.h

*/
