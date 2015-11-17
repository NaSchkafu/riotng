/*

 File: shrinkwrap.cpp, 2005/02/21

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

#include "shrinkwrap.h"
#include "matinv.h"

#include <string>

#include "adaptintval.h"
namespace riot{
void ShrinkWrapping::prepareAlgorithm(const Data *D)
{
  delete [] VarCodeList_;
  delete [] MonomList_;

  Dimension_ = D->Dimension_;

  InitValsTM_  = VVTaylorModel(Dimension_);
  MonomList_   = new Monom[Dimension_];
  VarCodeList_ = new unsigned int[Dimension_];

  for(unsigned int i = 0; i < Dimension_; i++)
    {
      InitValsTM_ [i] = TaylorModel( D->InitVals_Identifier_[i], 0.0, Interval(-1,1) );
      MonomList_  [i] = Monom( D->InitVals_Identifier_[i] );
      VarCodeList_[i] = Variables::getTable().look_up( D->InitVals_Identifier_[i] );
    }

  //Set the size of data.
  Constant_            = Vector(Dimension_);
  Linear_Mat_          = Matrix(Dimension_,Dimension_);
  Linear_Mat_Inv_      = Matrix(Dimension_,Dimension_);
  Linear_Mat_Inv_Inv_  = IMatrix(Dimension_,Dimension_);
  Nonlinear_           = VVTaylorModel(Dimension_);
  Remainder_           = IVector(Dimension_);

  AlgorithmPrepared = true;
}

/*

The following functions controls the Shrink wrapping process.

*/

void ShrinkWrapping::run(VVTaylorModel *T,const Files *file)
{
  /*
    The algorithm of this function is from the COSY-VI 
    package (3/22/2004), file 'tm.fox', see lines 806 to 856.
  */

  if( AlgorithmPrepared )
    {
      unsigned int kind = 0; //Will be used for output purposes.
      
      /*
	Split the Taylor model in constant, linear and remainder.
	'Nonlinear_' contains the modified Taylor model.
      */
      
      splitVVTM( T );
      
      /*
	Try shrink wrapping.
      */
      
      int err_code = 0;
      wrap(T,err_code);
      
      /*
	Calculate the condition number of the linear part.
      */

      double ConditionNumber = 0.0;

      if( err_code != 2 ) ConditionNumber = max(abs(Linear_Mat_))*max(abs(Linear_Mat_Inv_));
      
      /*
	Look if shrink wrapping was possible and if not
	calculate either a box or a parallelepiped enclosure.
      */
      
      double areaRemainder = 0.0, areaParallelepiped = 0.0, areaBox = 0.0;

      if( err_code != 0 ) 
	{
	  /*
	    First we check if we really have to use a box or 
	    parallelepiped enclosure or if we can use the 
	    given Taylor model for the next integration 
	    step without doing Shrink wrapping.

	    We calculate the area of the remainder box.
	  */

	  double fac   = 0.01; //Needed for the condition later (see below).
	  double area  = 2.0;  //Lenght of interval [-1,1].
	  
	  areaRemainder = ((*T)[0].interval_part()).diam();
	  for(int i = 1; i < Dimension_; i++)
	    {
	      fac  *= 0.01;
	      area *= 2.0;
	      areaRemainder *= ((*T)[i].interval_part()).diam();
	    }
	  
	  /*
	    Next we calculate the area of the parallelepiped given trough
	    the linear part of 'T'.
	  */

	  int err_code_tmp = 0;
	  areaParallelepiped = det(Linear_Mat_,err_code_tmp) * area;

	  /*
	    Now we check if the area of the remainder box is in relation 
	    to the parallelepiped area too big.
	  */

	  if( err_code_tmp == 0 and areaRemainder < fac * areaParallelepiped ) //and ConditionNumber < 1000 )
	    {
	      kind = 1; // '1' for no wrapping.
	    }
	  else
	    /*
	      An error occurded while calculating the parallelepiped area or the
	      remainder term is relative to the parallelepiped area much bigger in size.
	      
	      So we calculate a parallelepiped and a box enclosure and take
	      the one with the smaller area.
	    */
	    {
	      /*
		Calculate the box and the corresponding area.
	      */

	      IVector Box(Dimension_);
	      areaBox = 1.0; //Reset.
	      for(int i=0; i < Dimension_; i++)
		{ 
		  Box[i] = (*T)[i].eval();
		  areaBox *=  Box[i].diam();
		}
	      
	      /*
		Calculate the parallelepiped enclosure and the corresponding area if
		possible.
	      */

	      Vector Q(Dimension_);
	      Vector D(Dimension_);

	      if( err_code != 2 ) //and ConditionNumber < 1000 )
		/*
		  Matrix inversion didn't fail. The calculation of the parallelepiped
		  is possible.
		*/
		{
		  /*
		    Shrink wrapping was not possible because of the influence of the 
		    nonlinear part. So now we calculate an enclosure of the nonlinear 
		    part of 'T', then the shrink factor 'q' and after that the 
		    enclosing parallelepiped described by 'q * Linear_Mat_ * x'.
		    
		    'Remainder_' contains 'Linear_Mat_Inv_ * Remainder_' and 'Nonlinear_'
		    contains 'Linear_Mat_Inv_ * Nonlinear_'.
		  */

		  Remainder_ += Nonlinear_.eval();
		  D = sup(abs(Remainder_));
		  
		  for(int i=0; i < Dimension_; i++)
		    { 
		      Q[i] = 1.0 + D[i];
		      areaParallelepiped *= Q[i]; //Because: det(q*A) = q^n * det(A).
		    }

		  d = max(D);
		  q = max(Q);
		}
	      else areaParallelepiped = areaBox + 1.0; //Take the box enclosure (see next if-statement).

	      if( areaBox < areaParallelepiped )
		/*
		  Calculate the box representation with Taylor models each defined on [-1,1].
		*/
		{
		  kind = 2; // '2' for box enclosure.

		  for(int i=0; i < Dimension_; i++)
		    { 
		      double mid_point = Box[i].mid();
		      
		      Box[i] -= mid_point;
		      
		      double scale = Box[i].sup();
		      
		      (*T)[i] = mid_point + scale * InitValsTM_[i];
		    }
		}
	      else
		/*
		  Calculate the parallelepiped representation.
		*/
		{
		  kind = 3; // '3' for parallelepiped representation.

		  /*
		    First compute an enclosure of the inverse of 'Linear_Inv_'.
		  */
		  
		  imatinv( Linear_Mat_Inv_, Linear_Mat_Inv_Inv_, err_code_tmp );
		  
		  /*
		    Transform back.
		  */
		  
		  (*T) = (Linear_Mat_Inv_Inv_ * InitValsTM_);
		  
		  for(unsigned int i = 0; i < Dimension_; i++)
		    {
		      std::cout << "Q["<<i<<"] = " << Q[i] << std::endl;
		      
		      (*T)[i] *= Q[i];
		      (*T)[i] += Constant_[i];
		    }
		}
	    }
	}
			    
      /*
	Write the calculated values in the output file.
      */
      
      fprintf( 
	      (*file)[1],
	      
	      "%10.9E %10.9E %10.9E %10.9E %10s %10s",
	      
	      s,t,d,q-1.0, (err_code != 0)? "yes":"no", (kind < 2)? "--  ":((kind == 2)? "Box":"Epiped")
	      );
      fprintf( (*file)[2],"%15.9E %15.9E %15.9E",ConditionNumber,areaParallelepiped,areaBox );
      fflush ( (*file)[1] );
      fflush ( (*file)[2] );
    }
  else
    std::cout << "ShrinkWrapping::run(): *** Call prepareAlgorithm first. ***" << std::endl;
}

/*

Der Algorithmus des Shrink Wrapping. Verwendet die Toolbox-Funktionen zur Matrix Invertierung.

*/

void ShrinkWrapping::wrap(VVTaylorModel *T, int& err_code)
{
  /*
    First we have to compute the approximate inverse of 'Linear_'.
  */

  MatInv( Linear_Mat_, Linear_Mat_Inv_, err_code ); //Call to Toolbox function.
      
  if( err_code == 2 )
    {
      std::cout << "Matrix may be singular." << std::endl;
      return; //Quit shrink wrapping.
    }

  /*
    After that we compute a Taylor model for the nonlinear part of the transformed 'T'.
    
    T(x_1,...,x_n) = Linear_Mat_ * (x_1,...,x_n) + Nonlinear_(x_1,...,x_n)        (*)
    
    so in exact arithmetic we get
    
    Linear_Mat_Inv_ * T(x_1,...,x_n) = (x_1,...,x_n) + Linear_Mat_Inv_ * Nonlinear_(x_1,...,x_n).
    
    We compute the nonlinear part of the transformed 'T' as follows
    
    Linear_Mat_Inv_ * Nonlinear_(x_1,...,x_n) = Linear_Mat_Inv_ * T(x_1,...,x_n) - (x_1,...,x_n)
    
    to catch also linear corrections due to the error in inversion of the linear part.

    NOTE: 'Nonlinear_' contains 'T' as shown above (see (*)).
  */

  Nonlinear_ = Linear_Mat_Inv_ * Nonlinear_ - InitValsTM_;

  /*
    Transformation of remainder.
  */
	  
  Remainder_ = Linear_Mat_Inv_ * Remainder_;
      
  /*
    In the next step we compute an upper bound 's' of the absolute value of the nonlinear part, e.g.
    
    s >= | Nonlinear_[i](x_1,...,x_n) |    for all i=1,..,n, and all (x_1,...,x_n),
    
    and an upper bound 't' for the absolute value of all partial derivations of order one, e.g.
    
    t >= | D Nonlinear_[i](x_1,...,x_n) / D x_j |   for all i,j=1,..,n, and all (x_1,...,x_n).
  */
	  
  s = max(sup(abs(Nonlinear_.eval())));
  t = 0.0;
  for(int i = 0; i < Dimension_; i++)
    {
      VVTaylorModel Derivative = derivate( Nonlinear_, VarCodeList_[i] );
      
      double tmp = max(sup(abs(Derivative.eval())));
      
      t = ( t < tmp ) ? tmp : t;
    }
  Vector D = sup(abs(Remainder_));
  d = max(D);
      
  /*
    Look if the map is shrinkable and if so, calculate the shrink factor q.
  */

  Vector Q(Dimension_);
  q = 2.0;
  if( s < 1 and (Dimension_-1)*t < 1 ) 
    { 
      for(unsigned int i = 0; i < Dimension_; i++)
	{
	  Q[i] = sup( 1 + D[i] * (1 + (Dimension_-1)*Interval(t)) / ( (1-Interval(s)) * (1 - (Dimension_-1)*Interval(t)) ) );
	}

      q = max(Q);
    }
  if( q >= 1.01 ) 
    {
      err_code = 1;
      return; //Quit shrink wrapping.
    }

  /*
    Transform back.
    First compute an enclosure of the inverse of 'Linear_Inv_'.
  */
  
  imatinv( Linear_Mat_Inv_, Linear_Mat_Inv_Inv_, err_code );

  /*
    Now we can transform back and add the constant.
  */
  
  (*T) = (Linear_Mat_Inv_Inv_ * (InitValsTM_ + Nonlinear_));

  for(unsigned int i = 0; i < Dimension_; i++)
    {
      (*T)[i] *= Q[i];
      (*T)[i] += Constant_[i];
    }
}

/*

The following function split the given vector valued Taylor model in 
the constant, the linear and the remainder part.

*/

void ShrinkWrapping::splitVVTM(VVTaylorModel* T)
{
  Nonlinear_ = *T;

  for(int i = 0; i < Dimension_; i++)
    {
      Constant_ [i] = Nonlinear_[i].read_coefficient_and_remove( Monom() );

      for(int j = 0; j < Dimension_; j++) 
	Linear_Mat_(i,j) = Nonlinear_[i].coefficient_of( MonomList_[j] );
      
      Remainder_[i] = Nonlinear_[i].interval_part();

      Nonlinear_[i].replace_interval_part_with( Adapt::ZERO_INTERVAL() );
    }
}

/*

  End of File: shrinkwrap.cpp

*/

}
