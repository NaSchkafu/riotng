/*

 File: data.cpp, 2005/02/24

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

#include "data.h"

/*

This file contains all data for enclosing the solution of an initial
value problem. The data consists of:  

- The right hand side of the differential equation. Corresponding to the number
  of equations we have to define functions 'RHS_Equation_No_i'. The letter 'i' 
  stands for the number of the equation. The arguments of that functions are a
  vector 'x' with Taylor models and a Taylor model 't' which represents the 
  time variable. The index of the vector starts at 0 and 'x[i]' stands for the 
  unknown in the i-th equation. The expression of each function could consist 
  of the following operators and intrinsic functions: 
      
      - the unary operator {-};
      - the binary operators {+,-,*,/};
      - the intrinsic functions:
          
            sqr    ( TaylorModel      ),
	    sqrt   ( TaylorModel      ),
	    invsqrt( TaylorModel      ), (this is 1 / sqrt(...)),
	    exp    ( TaylorModel      ),
	    sin    ( TaylorModel      ),
	    cos    ( TaylorModel      ),
	    power  ( TaylorModel, int ).
*/
namespace riot {
TaylorModel& TemporaryVariable()
{
  static TaylorModel tmp;
  return tmp;
}

TaylorModel RHS_Equation_No_0 (const VVTaylorModel& y, const TaylorModel& t)
{
  return y[2];
}
TaylorModel RHS_Equation_No_1 (const VVTaylorModel& y, const TaylorModel& t)
{
  return y[3];
}
TaylorModel RHS_Equation_No_2 (const VVTaylorModel& y, const TaylorModel& t)
{
  TemporaryVariable()=-0.1000000000e-1*(200.*sin(y[1])*y[2]*y[3]+100.*sin(y[1])*sqr(y[3])-1962.*sin(y[0])+100.*cos(y[1])*sin(y[1])*sqr(y[2])+981.*sin(y[0])*sqr(cos(y[1]))+981.*cos(y[1])*cos(y[0])*sin(y[1])+100.*sin(y[1])*sqr(y[2]))/(-2.+sqr(cos(y[1])));
  return TemporaryVariable();
}
TaylorModel RHS_Equation_No_3 (const VVTaylorModel& y, const TaylorModel& t)
{
  TemporaryVariable()=0.1000000000e-1*(200.*cos(y[1])*sin(y[1])*y[2]*y[3]+100.*cos(y[1])*sin(y[1])*sqr(y[3])+981.*sin(y[0])*sqr(cos(y[1]))+981.*cos(y[1])*cos(y[0])*sin(y[1])+200.*sin(y[1])*y[2]*y[3]+100.*sin(y[1])*sqr(y[3])-1962.*sin(y[0])+1962.*cos(y[0])*sin(y[1])+300.*sin(y[1])*sqr(y[2])+200.*cos(y[1])*sin(y[1])*sqr(y[2]))/(-2.+sqr(cos(y[1])));
  return TemporaryVariable();
}
/*TaylorModel RHS_Equation_No_4 (const VVTaylorModel& x, const TaylorModel& t)
{
  return TemporaryVariable() * x[1];
}
TaylorModel RHS_Equation_No_5 (const VVTaylorModel& x, const TaylorModel& t)
{
  return TemporaryVariable() * x[2];
}*/


Data::Data()
{

                                                                               /*

- the number of equations. Must be equal to the number of the defined functions
  above.

                                                                               */
  Dimension_ =  4;

                                                                               /*

- the order of the Taylor models.

                                                                               */
  Order_     =  10;

                                                                               /*

- the way of how to determine the order. Possible are:

      - TOTALDEGREE : sum of the exponents of all variables in a term.

                                                                               */

  Order_Check_ = TOTALDEGREE;

                                                                               /*

- the algorithm for calculating the range enclosures of polynomials. You can 
  choose between:
             
      - LDB   : for the Linear Dominated Bounder;
      - NAIVE : for the naive interval evaluation;
      - MVF   : for the mean value form.

                                                                               */

  Bounder_ = LDB;

                                                                               /*

- the information if you want to use Shrink wrapping (choose 'ON') 
  or if you don't want to use it (choose 'OFF').

                                                                               */

  Shrink_Wrapping_ = ON;


                                                                               /*

- the bound for generating sparsity.

                                                                               */

  Sparsity_ = 1e-15;

                                                                               /*

- the step size control:

     - AUTO  : for automatic step size control depending on the local 
               error but with minimal step size.
     - CONST : for constant step size. You can set one step size or 
               a sequel of step sizes. In the first case the end point 
	       of integration will be taken into account, in the second
	       case the end position will be computet as sum of the 
	       starting point with the all step sizes added.

                                                                               */
  
  Step_cntrl_ = AUTO;

  //
  // In case of automatic step size control set the following parameters.
  //

  h_start_ = 1e-3; 
  h_min_   = 1e-4; 

  //
  // In case of a constant step size give either one step size or a 
  // sequel of step sizes. In the second case separate the step
  // sizes with commas but set no comma behind the last step size.
  //

  /*Step_sizes_ = 1;              
  double local_h_sequel[] = //Only for a simpler setting of step sizes 
    {                       //(e.g. with copy and paste).
      0.1E-3
    };

  h_sequel_ = new double[Step_sizes_];
  for(unsigned int i = 0; i < Step_sizes_; i++) h_sequel_[i] = local_h_sequel[i];*/
  
                                                                               /*

- a bound for the local error in each step size. If you have chosen the 
  automatic step size control and the done error is greater
  than the given bound, then the integration step will be repeated.

                                                                               */

  Local_error_tol_ = 1e-11;

                                                                               /*

- the initial point of the initial value problem.

                                                                               */
  t_0_ = 0;

                                                                               /*

- the points where to calculate the solution enclosure. You can specify
  more than one point, that means that you can add intermediate points at 
  which the result will be calculated and printed.
  If the grid point is not machine representable use the 
  form ' "[a,a]" >> result_t_[i] '! Otherwise you can use 
  ' result_t_[i] = Interval(a,a) '. Set the points in ascending order!

                                                                               */

  Result_points_ = 1;//How many points?
  result_t_ = new Interval[ Result_points_ ]; 
  "[1,1]" >> result_t_[0];



/*

- a letter for the time variable.

                                                                               */

  Time_Identifier_ = "t";

                                                                               /*

- and last but not least letters for and values of the initial values.

                                                                               */
  

  InitVals_Identifier_ = new std::string[ Dimension_ ];
  InitVals_Value_      = IVector(Dimension_);

  InitVals_Identifier_[0]="x0";
  InitVals_Identifier_[1]="x1";
  InitVals_Identifier_[2]="x2";
  InitVals_Identifier_[3]="x3";

  "[2.3561925,2.3561925]">>InitVals_Value_[0];    
  "[-1.7265335380000000,-1.7265335380000000]">>InitVals_Value_[1];   
  "[0.43,0.43]">>InitVals_Value_[2]; 
  "[0.670,0.670]">>InitVals_Value_[3];

                                                                               /*

At the end we must set the pointers to the functions definied above and give 
the problem/data a name for naming the output files. 

                                                                               */

  Func_Ptr_RHS_ = new RHS[ Dimension_ ];
  Func_Ptr_RHS_[0]=&RHS_Equation_No_0;
  Func_Ptr_RHS_[1]=&RHS_Equation_No_1;
  Func_Ptr_RHS_[2]=&RHS_Equation_No_2;
  Func_Ptr_RHS_[3]=&RHS_Equation_No_3;
  

  Title_ = "MDP";

  out_prec_ = 10; //Set the precision of the taylorcoefficients for
                  //output in file 'Title'_03.txt
}

/*

  End of File: data.cpp

*/
}
