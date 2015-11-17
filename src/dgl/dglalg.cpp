/*

 File: dglalg.cpp, 2006/03/08

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

#define OUTPUT

#include "dglalg.h"

#include "polyeval.h" //From TM-Base package.

#include "timeid.h"
#include "stopwatch.h"
#include "shrinkwrap.h"

#include <stdlib.h>
namespace riot{
void DGLSolver::prepareAlgorithm(Data *D)
{
  Data_Ptr_ = D;

  /*
    Set output precision.
  */
  
  //
  // For output of 'double'.
  //
  std::cout.precision(15);
  std::cout.setf( std::ios::scientific );
#ifdef FILIB_VERSION
  //
  // For output of 'Interval'.
  //
  filib::interval<double>::precision(15);
#endif

  /*
    Set order of Taylor model.
  */

  if( D->Order_ > 0 )
    {
      TaylorModel::set_order( D->Order_ );
    }
  else
    {
      std::cout << "ERROR: *** Order_ must be greater than zero *** " << std::endl;
      std::exit(1); //Quit.
    }

  /*
    Set method for checking the order of the polynomial, i.e. 
    for checking the structure of the monomials.
  */

  switch( D->Order_Check_ )
    {
    case TOTALDEGREE: 
      {
	TaylorModel::set_degree_check(Total_Degree::getObject());
	break;
      }
    default:
      {
	std::cout << "ERROR: *** The chosen way of checking the polynomial order is not available. ***" << std::endl;
	std::exit(1); //Quit.
      }
    }

  /*
    Set polynomial range evaluation.
  */

  switch( D->Bounder_ )
    {
    case LDB: 
      {
	TaylorModel::set_polynomial_range_evaluation(LDB::getObject()); 
	break;
      }
    case NAIVE:
      {
	TaylorModel::set_polynomial_range_evaluation(IntervalEvaluation::getObject()); 
	break;
      }
    case MVF:
      {
	TaylorModel::set_polynomial_range_evaluation(MeanValueForm::getObject());
	break;
      }
    default:
      {
	std::cout << "ERROR: *** The chosen algorithm for polynomial range evaluation is not available. ***" << std::endl;
	std::exit(1); //Quit.
      }
    }

  /*
    Set sparsity tolerance.
  */
  
  if( 0 < D->Sparsity_ && D->Sparsity_ < 1E-5 )
    {
      TaylorModel::set_sparsity_tol( D->Sparsity_ );
    }
  else
    {
      std::cout << "ERROR: *** Sparsity_ must be positive and lower than 1E-5. *** " << std::endl;
      std::exit(1); //Quit.
    }

  /*
    Check step size control parameters and set 't_end_'.
  */

  switch( D->Step_cntrl_ )
    {
    case AUTO:
      {
	h_curr_ = D->h_start_;

	if( h_curr_ <= 0.0 or D->h_min_ < 0.0 )
	  {
	    std::cout << "ERROR: *** Wrong value in h_start_ or h_min_. *** " << std::endl;
	    std::exit(1); //Quit.
	  }
	else if( D->Result_points_ <= 0 )
	  {
	    std::cout << "ERROR: *** Number of result points has wrong value. *** " << std::endl;
	    std::exit(1); //Quit.
	  }
	else
	  {
	    /*
	      Chose last point as end point.
	    */

	    t_end_ = D->result_t_[ D->Result_points_ - 1 ];
	  }
	break;
      }
    case CONST:
      {
	if( D->Step_sizes_ <= 0 )
	  {
	    std::cout << "ERROR: *** Wrong number of step sizes. *** " << std::endl;
	    std::exit(1); //Quit.
	  }
	else if( D->Step_sizes_ == 1 )
	  {
	    if( D->Result_points_ <= 0 )
	      {
		std::cout << "ERROR: *** Number of result points has wrong value. *** " << std::endl;
		std::exit(1); //Quit.
	      }
	    else
	      {
		/*
		  Chose last point as end point.
		*/

		t_end_ = D->result_t_[ D->Result_points_ - 1 ];
	      }
	  }
	else
	  {
	    /*
	      Chose sum of given step sizes as end point.
	    */

	    double sum = D->h_sequel_[0];
	    for(unsigned int i = 1; i < D->Step_sizes_; i++) sum += D->h_sequel_[i];

	    t_end_ = D->t_0_ + sum;
	  }

	h_curr_ = D->h_sequel_[0];
	break;
      }
    default:
      {
	std::cout << "ERROR: *** The chosen step size control is not available. ***" << std::endl;
	std::exit(1); //Quit.
      }
    }

  /*
    Check and set the algorithm specific data.
  */

  if( D->t_0_ >= sup(t_end_) )
    {
      std::cout << "ERROR: *** Wrong value in t_0. *** " << std::endl;
      std::exit(1); //Quit.
    }
  else
    {
      t_curr_ = D->t_0_;
      t_next_ = D->t_0_;

      InitValsTM_    = VVTaylorModel(D->Dimension_);
      InitValsDGLTM_ = VVTaylorModel(D->Dimension_);
      currTM_        = VVTaylorModel(D->Dimension_);
      Solution_      = VVTaylorModel(D->Dimension_);
      
      for(int i = 0; i < D->Dimension_; i++) 
	{
	  InitValsTM_[i] = TaylorModel( D->InitVals_Identifier_[i], 0.0, Interval(-1,1) );
	}

      //
      // Prepare the shrink wrapping algorithm.
      //
      ShrinkWrapping::Algorithm().prepareAlgorithm( D );

      //
      // Set the time identifier.
      //
      TimeCode_ = TimeIdentifier::getObject().setID(D->Time_Identifier_);

      currIV_ = IVector(D->Dimension_);
    }

  /*
    Create a copy of the data file.
  */

  system( ("cp data.cpp " + Data_Ptr_->Title_ + "_Data.txt").c_str() );
  
  AlgorithmPrepared = true;
}

void DGLSolver::run(const Files *file)
{
  if( AlgorithmPrepared )
    {
      Stopwatch watch_0, watch_1, watch_2, watch_3, watch_4, watch_5, watch_6;

      /*---------------------Prepare output files: write the headers.-------------------------*/
      
      fprintf( (*file)[0],"%-8s %15s %15s %6s %8s %6s %15s %11s %10s %10s %10s %7s %7s \n",
	       "Step-Nr.","t_curr_","h_curr_     ","iter_1","iter_sum","iter_2","local_err    ","h decreased",
	       "Time Step","Time TM","Time I0","Time I*","Time J*"
	       );
      fflush ( (*file)[0] );
      fprintf( (*file)[1],"%-8s %15s %15s %15s %15s %10s %10s %10s \n",
	       "Step-Nr.","s        ","t        ","d        ","q        ","1# failed", "Encl.", "Time "
	       );
      fflush ( (*file)[1] );
      fprintf( (*file)[2],"%-8s %15s %15s %15s \n",
	       "Step-Nr.","cond     ","ParArea  ","BoxArea   ");
      fflush ( (*file)[2] );

      /*
	Write init values for plotting.
      */

      fprintf( (*file)[3], "%s %+16.15E ", "*", t_curr_ );
      for(int i = 0; i < Data_Ptr_->Dimension_; i++) 
	{
	  fprintf( (*file)[3], "%+16.15E %+16.15E ", 
		   Data_Ptr_->InitVals_Value_[i].inf(), 
		   Data_Ptr_->InitVals_Value_[i].sup() );
	}
      fprintf( (*file)[3], "\n" );
      fflush ( (*file)[3] );

      /*--------------------------------------------------------------------------------------*/

      const IVector ZERO_IVECTOR(Data_Ptr_->Dimension_);
      IVector Enclosure(Data_Ptr_->Dimension_);
      unsigned int step_cntr = 0;
      unsigned int t_cntr    = 0;

      watch_0.start();

      /*
	The variable 't_curr_' denotes the current grid point and the
	variable 't_next_' the following grid point. The values are equal 
	when entering the loop.
	The variable 't_end_' is of type interval and contains the possibly 
	not machine representable right end of the domain for integration.
	'Data_Ptr_' is a pointer to the data structur.
      */
      //CM1>
      bool Repeat_Step = false;
      bool TransformInitVals = true;

      while( t_curr_ < sup(t_end_) )
	{
	  ++step_cntr;    //Counts steps of integration.

	  //<CM1
	  watch_1.start();//Measures time of the current integration step.
#ifdef OUTPUT
	  std::cout << step_cntr << " ------------------" << std::endl;
#endif	  
	  /*
	    The initial values will be represented with Taylor models each defined on 
	    [-1,1]. Is the scaling factor 'scale' lower than the user defined
	    sparsity tolerance 'Sparsity_' the Taylor model stays constant, because
	    internal the variable term will be swept into the remainder term.
	    The variable 'InitValsTM_' is a vector with Taylor models each of the form
	    [ \alpha_i , [0,0] ], \alpha_i \in [-1,1]. 
	    'InitValsDGLTM_' constains the Taylor models representing 
	    the initial value set in the current time step.
	    The variable 'Dimension_' keeps the size of the ODE-System.
	  */
	  //CM2>
          if( TransformInitVals and not Repeat_Step )
            {
              for(int i=0; i < Data_Ptr_->Dimension_; i++)
                { 
		  double mid_point = Data_Ptr_->InitVals_Value_[i].mid();
                  Interval domain = Data_Ptr_->InitVals_Value_[i] - mid_point;
                  double scale = domain.sup();
		  
                  TransformInitVals &= (scale < Data_Ptr_->Sparsity_);
                  InitValsDGLTM_[i] = mid_point + scale * InitValsTM_[i];
                }
	    }
	  //<CM2//CM3>
	  /*
	    Set 't_next_' and check if the step was too large.
	  */

	  t_next_ = t_curr_ + h_curr_;
	  
	  if( t_next_ >= inf(t_end_) ) 
	    {
	      t_next_ = sup(t_end_);

	      /*
		Calculate new step size.
	      */

	      h_curr_ = t_next_ - t_curr_;
	    }
	  //<CM3//CM4>
	  /*
	    Next we generate a Taylor model for the time variable:
	    TimeTM_=[ t_curr_+(t-t_curr_) , [0,0] ], 
	    t \in [t_curr_,t_next_].
	  */

	  Interval time_interval(t_curr_,t_next_);

	  if( Repeat_Step )
	    {
	      TimeTM_.change_var_data(TimeCode_, t_curr_, time_interval);
	    }
	  else
	    {
 	      TimeTM_ = t_curr_ + 
 		TaylorModel( 
			    Data_Ptr_->Time_Identifier_, 
			    t_curr_, 
			    time_interval 
			   );

	      /*
		Initialize 'Solution_' and store information about the time 
		variable (the subdomain of integration) in it.
		This will be needed for integration.
	      */
	      
	      Solution_ = InitValsDGLTM_ + 0.0 * TimeTM_;
	    }
	  //<CM4
	  /*
	    Calculate the Taylor polynomial of the Solution with 
	    Picard Iteration.
	    
	    Start with polynomial order 1 and increase it by one in each step
	    till the user defined polynomial order is reached. This saves
	    unnecessary calculations and storage of terms.
	  */
	  //CM5>
	  if( not Repeat_Step )
	    /*
	      Calculate the Taylor polynomial.
	    */
	    {
	      //<CM5

	      watch_2.start();
	      
	      //CM5>
	      TaylorModel::set_mode( NO_REMAINDER );
	      
	      for(int k = 1; k <= Data_Ptr_->Order_; k++)
		{
		  TaylorModel::set_order( k );

		  /*
		    Evaluate the Picard operator.
		  */

		  for(int i = 0; i < Data_Ptr_->Dimension_; i++)
		    {
		      currTM_[i] = InitValsDGLTM_[i] +
			integrate( 
			       Data_Ptr_->Func_Ptr_RHS_[i]( Solution_, TimeTM_ ), 
			       TimeCode_ 
			      );
		    }

		  Solution_ = currTM_;
		}
	      //<CM5

	      watch_2.stop();

	      //CM5>
	    }
	  else
	    /*
	      The integration step will be repeated. We can reuse 
	      the already calculated Taylor polynomial.
	    */
	    {
	      Solution_ = currTM_;
	    }

	  TaylorModel::set_mode( NORMAL_MODE );
	  //<CM5//CM6>
	  int  iter_1, iter_1_sum = 0, iter_1_max = 3;
	  bool Inclusion = true;
	  bool Decreased = false;
	  TaylorModel help;

	  do {
	    //<CM6//CM7>
	    /*
	      Calculate the inner approximation.
	    */
	    //<CM7

	    watch_3.start();

	    //CM7>
	    // Create P+[0,0].
	    Solution_.replace_interval_part_with( ZERO_IVECTOR );
	    TaylorModel::set_mode( REMAINDER_REC ); // Set Record mode.

	    for(int i = 0; i < Data_Ptr_->Dimension_; i++)
	      {
		/*
		  Evaluate the i-th component of the
		  integral operator: K_i( P+[0,0] ) - P_i.
		  The calculated inner approximation will be stored
		  in the interval vector 'currIV_'.
		*/

		help  = InitValsDGLTM_[i];
		help += integrate( 
			       Data_Ptr_->Func_Ptr_RHS_[i]( Solution_, TimeTM_ ),
			       TimeCode_ 
			      );
		help.subtract_polynomial_part_of( Solution_[i] );

		currIV_[i] = help.eval();
	      }
	    //<CM7

	    watch_3.stop();

	    //CM8>
	    /*
	      Calculate enclosures of the solution.
	    */

	    //<CM8

	    watch_4.start();

	    //CM8>
	    TaylorModel::set_mode( REMAINDER_PLAY ); //Play mode.

	    iter_1 = 0;
	    Inclusion = false;
	    
	    while( not Inclusion and iter_1 < iter_1_max )
	      {
		++iter_1;

		/*
		  Replace the remainder part with the inflated
		  inner approximation I_k.
		*/

		for(int i = 0; i < Data_Ptr_->Dimension_; i++) 
		  {
		    inflate( currIV_[i] ); //Epsilon Inflation.

		    Solution_[i].replace_interval_part_with( currIV_[i] );
		  }

		/*
		  Evaluate the integral operator: K( P + I_k ) - P =: J_k
		*/

		Inclusion = true;

		for(int i = 0; i < Data_Ptr_->Dimension_; i++)
		  {
		    help  = InitValsDGLTM_[i];
		    help += integrate( 
				   Data_Ptr_->Func_Ptr_RHS_[i]( Solution_, TimeTM_ ), 
				   TimeCode_ 
				  );
		    help.subtract_polynomial_part_of( Solution_[i] );

		    currIV_[i] = help.eval();

		    /*
		      Check if there is inclusion.
		    */

		    Inclusion &= (currIV_[i]).subset( Solution_[i].interval_part() );

		    /*
		      If inclusion holds in the current component
		      use the calculated (better) enclosure in the following
		      calculations.
		    */ 

		    if( Inclusion ) 
		      Solution_[i].replace_interval_part_with( currIV_[i] );
		  }
	      }
	    //<CM8

	    watch_4.stop();

	    //CM9>
	    /*
	      If failure then decrease the step size or in 
	      case of constant step size exit the algorithm.
	    */

	    if( not Inclusion ) 
	      {
		/*
		  No inclusion has been reached, because the time interval
		  of the integration step was too big.
		  Decrease the step size in case of automatic step 
		  size control and exit in case of a constant step size.
		*/
		
		switch( Data_Ptr_->Step_cntrl_ )
		  {
		  case AUTO:
		    {
		      h_curr_ *= 0.7; // Use the same factor as in COSY.

		      if( h_curr_ < Data_Ptr_->h_min_ )
			{
			 std::cout << "DROPOUT: *** "
			       << "Step size undergoes minimal step size." 
			       << "Inclusion is not possible. *** " 
			       << std::endl;
			 std::exit(1); //Quit.
			}

		      Decreased = true;
		
		      t_next_ = t_curr_ + h_curr_;

		      /*
			Change the data of the time variable.
		      */
		      
		      TimeTM_.change_var_data(
					      TimeCode_, 
					      t_curr_, 
					      Interval(t_curr_,t_next_)
					     );
		      
		      break;
		    }
		  case CONST:
		    {
		     std::cout << "DROPOUT: ***"
			   << "With the given step size"
			   << "inclusion is not possible. *** " 
			   << std::endl;
		     std::exit(1); //Quit.

		     break;
		    }
		  }
	      }

	    iter_1_sum += iter_1;

	  } while ( not Inclusion );
	  //<CM9

	  watch_5.start();

	  //CM10>
	  /*
	    Improve the enclosure of the solution.
	  */

	  TaylorModel::set_mode( REMAINDER_PLAY ); //Play mode.

	  int iter_2 = 0;
	  Interval better_enc;

	  do {

	    /*
	      Evaluate the integral operator.
	    */

	    for(int i = 0; i < Data_Ptr_->Dimension_; i++)
	      {
		//Save remainder interval.
		currIV_[i] = Solution_[i].interval_part(); 

		help  = InitValsDGLTM_[i];
		help += integrate( 
			      Data_Ptr_->Func_Ptr_RHS_[i]( Solution_, TimeTM_ ),
			      TimeCode_ 
			     );
		help.subtract_polynomial_part_of( Solution_[i] );
		
		better_enc = help.eval();
		
		/*
		  Set new remainder interval.
		*/
		
		Solution_[i].replace_interval_part_with( better_enc );
	      }

	    ++iter_2;

	  } while ( not terminate( currIV_, Solution_.interval_part() ) );

	  TaylorModel::set_mode( NORMAL_MODE ); //Normal mode.
	  //<CM10

	  watch_5.stop();

	  //
	  // Print some results.
	  //
	  fprintf( (*file)[0],"%-8u %15.9E %15.9E %6u %8u %6u ",step_cntr,t_curr_,h_curr_,iter_1,iter_1_sum,iter_2 );
	  fflush ( (*file)[0] );
	  
	  //CM11>
	  /*
	    Calculate the solution set at 't_next_' and 
	    the local error increase.
	  */	 

	  Interval grid_point;
	  double local_err = 0.0;

	  if( t_next_ == sup(t_end_) ) grid_point = t_end_;
	  else                         grid_point = t_next_;

	  for(int i = 0; i < Data_Ptr_->Dimension_; i++) 
	    {
	      help = substitute( Solution_[i], TimeCode_, grid_point );

	      Interval tmp = abs( 
				 subtract_bounds(
						 help.interval_part(),
						 InitValsDGLTM_[i].interval_part()
						) 
				);

	      local_err = ( local_err > tmp.sup() ) ? local_err : tmp.sup();
	    }
	  //<CM11
#ifdef OUTPUT
	  //std::cout << "ip= " << t_next_ << " " << wert << std::endl;//------------------------------
	  std::cout << "local_err = " << local_err << std::endl;
#endif
	  //CM12>
	  /*
	    If automatic step size control was chosen then 
	    decide on the local error increase if the time 
	    step should be repeated and correct the step 
	    size due to the local error increase. In case 
	    of a constant step size the size of the local 
	    error increase will be ignored.
	  */

	  switch( Data_Ptr_->Step_cntrl_ )
	    {
	    case AUTO:
	      {
		Interval tmp = pow( 
		    Interval(Data_Ptr_->Local_error_tol_) / (10*local_err), 
		    Interval(1.0)/(TaylorModel::order()+1) 
		   );
	  
		double c = tmp.inf(); //Factor for step size correction.

		if( local_err > Data_Ptr_->Local_error_tol_ 
		    and h_curr_ > Data_Ptr_->h_min_ ) 
		  /*
		    Repeat the time step.
		  */
		  {
		    Repeat_Step = true;
		    --step_cntr;

		    //
		    // Reset 't_next_'.
		    //
		    t_next_ = t_curr_;
		  }
		else
		  /*
		    The local error increase is smaller than the user 
		    defined upper bound. The time step need not be 
		    repeated.
		  */
		  {
		    Repeat_Step = false;
		    t_curr_     = t_next_;

		    //
		    // Make sure the factor for 
		    // step size correction is not too big.
		    //
		    if( c > 1.1 ) c = 1.1;
		  }

		/*
		  Correct the step size and make sure the step size is
		  not too small.
		*/

		h_curr_ *= c;
		if( h_curr_ < Data_Ptr_->h_min_ ) 
		  {
		    h_curr_ = Data_Ptr_->h_min_;
		  }

		if( c < 1.0 ) Decreased = true; //For output only.

		break;
	      }
	    case CONST:
	      {
		t_curr_ = t_next_;

		if( Data_Ptr_->Step_sizes_ > 1 )
		  h_curr_ = Data_Ptr_->h_sequel_[ step_cntr ];
	      }
	    }
	  //<CM12

	  //
	  // Print some results.
	  //
	  fprintf( (*file)[0],"%15.9E %11s",local_err,(Decreased?"yes":"no") );
	  fflush ( (*file)[0] );

	  /*
	    If the time step need not to be repeated then calculate the 
	    solution set at 't_next_' and at the intermediate grid points 
	    between 't_curr_' and 't_next_', do shrink wrapping and set 
	    the new init value set.
	  */
	  //CM13>

	  if( not Repeat_Step )
	    {
	      /*
		Calculate the solution set.
	      */
	      //<CM13

	      //
	      //Intermediate grid points first.
	      //
	      while( (t_cntr < Data_Ptr_->Result_points_) and
		     (sup(Data_Ptr_->result_t_[ t_cntr ]) <= t_next_) )
		{
		  grid_point = Data_Ptr_->result_t_[ t_cntr++ ];

		  //
		  // Print header of solution output.
		  //
		  fprintf( (*file)[3], "%s %+16.15E ", "*", grid_point.mid() );
		  fprintf( (*file)[4], "Step no. %8i:   t = %10.6f --------------------\n", 
			   step_cntr,
			   grid_point.mid() );

		  //
		  // Now calculate solution set.
		  //
		  for(int i = 0; i < Data_Ptr_->Dimension_; i++) 
		    {
		      help = substitute( Solution_[i], TimeCode_, grid_point );
		  
		      //
		      // Print the solution set and an interval enclosure.
		      //
		      Interval enclosure = help.eval();
		      fprintf( (*file)[3], "%+16.15E %+16.15E ", enclosure.inf(), enclosure.sup() );
		  
		      std::ostringstream out;
		      out.precision( Data_Ptr_->out_prec_ );
		      out.setf( std::ios::scientific );
		      out << help;
		      fprintf( (*file)[4], "%s\n", (out.str()).c_str() );
		    }

		  //
		  // Print foot of solution output.
		  //
		  fprintf( (*file)[3], "\n" );
		  fflush ( (*file)[3] );
		  fprintf( (*file)[4], "\n" );
		  fflush ( (*file)[4] );
		}
	      //CM13>

	      if( t_next_ == sup(t_end_) ) grid_point = t_end_;
	      else                         grid_point = t_next_;
	      //<CM13
	      //
	      // Print header of solution output.
	      //
	      fprintf( (*file)[3], "%s %+16.15E ", "-", grid_point.mid() );
	      //CM13>

	      for(int i = 0; i < Data_Ptr_->Dimension_; i++) 
		{
		  Solution_[i] = substitute( Solution_[i], TimeCode_, grid_point );

		  /*
		    'Enclosure' is an interval vector which bounds 
		    the solution set, described by the calculated 
		    Taylor models.
		  */

		  Enclosure[i] = Solution_[i].eval();
		  //<CM13
		  //
		  // Print the results for plotting.
		  //
		  fprintf( (*file)[3], "%+16.15E %+16.15E ", Enclosure[i].inf(), Enclosure[i].sup() );
#ifdef OUTPUT
		  std::cout << "Enclosure[" << i << "] = " << Enclosure[i] 
			    << ",  diam = " << Enclosure[i].diam() 
			    << std::endl;
#endif
		  //CM13>
		}
	      //<CM13

	      fprintf( (*file)[3], "\n" );
	      fflush ( (*file)[3] );

	      //CM13>

	      /*
		Use shrink wrapping only at intermediate grid points.
		When shrink wrapping at the end point the 
		enclosure property of the Taylor models for 
		fixed initial values is lost.
	      */
	      
	      if( Data_Ptr_->Shrink_Wrapping_ == ON and t_next_ < sup(t_end_) )
		{
		  //<CM13

		  fprintf( (*file)[1],"%-8u ",step_cntr );
		  fprintf( (*file)[2],"%-8u ",step_cntr );

		  watch_6.start();

		  //CM13>
		  ShrinkWrapping::Algorithm().run( &Solution_, file );
		  //<CM13

		  watch_6.stop();
		  
		  fprintf( (*file)[1],"%10.2f \n",watch_6.accumulatedTime() );
		  fprintf( (*file)[2],"\n" );
		  fflush ( (*file)[1] );
		  fflush ( (*file)[2] );

		  //CM13>
		}

	      /*
		An enclosure of the solution has been calculated.
		If the init values of the current time step are points,
		then use the intervals in 'Enclosure' for the following
		time step. Otherwise use the calculated representation
		of the shrink wrapped solution set by Taylor models.
	      */
	      
	      if( TransformInitVals )
		{
		  Data_Ptr_->InitVals_Value_ = Enclosure;
		}
	      else
		{
		  InitValsDGLTM_ = Solution_;
		}
	    }
	  //<CM13

	  watch_1.stop();

	  /*
	    Print the calculation times.
	  */

	  fprintf( (*file)[0],"%10.2f %10.2f %10.2f %7.2f %7.2f \n",
		   watch_1.accumulatedTime(), //Total time for the time step.
		   watch_2.accumulatedTime(), //Total time for calculating the Taylor polynomial.
		   watch_3.accumulatedTime(), //Total time for inner approximation.
		   watch_4.accumulatedTime(), //Total time for solution enclosure.
		   watch_5.accumulatedTime()  //Total time for iterative refinement.
		   );
	  fflush ( (*file)[0] );
	      
	  watch_1.clearAccu();
	  watch_2.clearAccu();
	  watch_3.clearAccu();
	  watch_4.clearAccu();
	  watch_5.clearAccu();
	  watch_6.clearAccu();
	  //CM13>
	  
	} //End of while-loop (Time-Loop).
      //<CM13

      watch_0.stop();

      /*
	In case of a sequel of step sizes we have to print
	the solution at the calculated end point. Notice: In all
	other cases of step size control the last grid point from 
	the grid point array is used as end point and so the 
	corresponding solution is already printed (see above).
      */

      if( Data_Ptr_->Step_cntrl_ == CONST and Data_Ptr_->Step_sizes_ > 1 )
	{
	  //
	  // Print header of solution output.
	  //
	  fprintf( (*file)[3], "%s %+16.15E ", "*", t_end_.mid() );
	  fprintf( (*file)[4], "Step no. %8i:   t = %10.6f --------------------\n", 
		   step_cntr,
		   t_end_.mid() );

	  //
	  // Print the solution set and an interval enclosure.
	  //
	  for(int i = 0; i < Data_Ptr_->Dimension_; i++) 
	    {
	      fprintf( (*file)[3], "%+16.15E %+16.15E ", Enclosure[i].inf(), Enclosure[i].sup() );
	  
	      std::string s;
	      s << Solution_[i];
	      fprintf( (*file)[4], "%s\n", s.c_str() );
	    }
	  
	  //
	  // Print foot of solution output.
	  //
	  fprintf( (*file)[3], "\n" );
	  fflush ( (*file)[3] );
	  fprintf( (*file)[4], "\n" );
	  fflush ( (*file)[4] );
	}

      /*
	Print the total calculation time.
      */
      
      fprintf( (*file)[0],"\n%s %10.2f","Rechenzeit insgesamt:",watch_0.accumulatedTime());
      fflush ( (*file)[0] );

      AlgorithmPrepared = false;
    }
  else
    std::cout << "DGLSolver::run(): *** Call prepareAlgorithm first. ***" << std::endl;
}

/*

Private functions for internal use only.

*/

/*
  
Termination criteria in iterative refinement.

*/

bool DGLSolver::terminate(const IVector& prev, const IVector& curr)
{
  static double eps = 0.01; // 1%
  bool ok = true;

  for(unsigned int i = 0; i < prev.size(); i++) 
    {
      Interval 
	x = curr[i],
	y = prev[i];
      double d = x.diam();

      if( d != 0.0 )
	{
	  ok &= ( ((x.inf()-y.inf())/d) < eps );
	  ok &= ( ((y.sup()-x.sup())/d) < eps );
	}
    }

  return ok;
}

/*

Epsilon-Inflation.

*/

void DGLSolver::inflate(Interval& x)
{
  static double eps1 = 1e-2;
  static double eps2 = 1e-6;

  x = (1+eps1)*x - eps1*x + Interval(-eps2,eps2);
}


/*

  End of File: dglalg.cpp

*/

// 	      Interval tmp = Interval( abs( Solution_[i].interval_part() ).sup() ) //Das berechnet Berz.
// 		- Interval( abs( InitValsDGLTM_[i].interval_part() ).sup() );

//	      wert = (wert < (sup(abs(help.interval_part())))) ? (sup(abs(help.interval_part()))) : wert;

}
