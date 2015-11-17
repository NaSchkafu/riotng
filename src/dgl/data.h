/*

 File: data.h, 2005/02/12

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

#ifndef DATA_H_INCLUDED
#define DATA_H_INCLUDED

#include <vector>
#include <string>

#include "adaptintval.h"
#include "taylormodel.h"

#include "ivector.h"
#include "vvtaylormodel.h"
#include "degree.h"

namespace riot {
enum StepSizeControl {
  AUTO  = 0,
  CONST = 1
};

enum BounderType {
  LDB   = 0,
  NAIVE = 1,
  MVF   = 2 
};

enum DegreeType {
  TOTALDEGREE = 0
};

enum Switch {
  ON  = 0,
  OFF = 1
};

typedef TaylorModel (*RHS)(const VVTaylorModel&,const TaylorModel&);

struct Data
{
  std::string          Title_; //Name of the problem. Used for naming the output files.

  unsigned int     Dimension_; //Dimension of the differential equation.
  unsigned int         Order_; //The order of the Taylor models.

  BounderType        Bounder_; //Algorithm for calculating a range enclosure of the polynomials.
  DegreeType     Order_Check_; //The way of how to determine the order of terms.
  Switch     Shrink_Wrapping_; //Shrink wrapping on or off.
  
  double            Sparsity_; //The sparsity tolerance.
  double                 t_0_; //Starting point.
  StepSizeControl Step_cntrl_; //Kind of step size control.
  double             h_start_; //Initial step size.
  double               h_min_; //Minimal step size.
  unsigned int    Step_sizes_; //Number of given step sizes.
  double           *h_sequel_; //Sequel of step sizes.
  double     Local_error_tol_; //Error tolerance for step size control.
  Interval         *result_t_; //Array with grid points where a solution should be printed.
  unsigned int Result_points_; //Number of grid point where a solution should be printed.

  std::string *InitVals_Identifier_; //Contains the names of the initial values.
  IVector      InitVals_Value_;      //Contains the values of the initial values.

  std::string Time_Identifier_; //Name of time variable.

  RHS *Func_Ptr_RHS_;

  unsigned int out_prec_;

  //Constructor which sets all needed information. See file 'data.cpp'.
  Data();
};

}
#endif

/*

  End of File: data.h

*/


