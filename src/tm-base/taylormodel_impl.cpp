/*

  File: taylormodel_impl.cpp, 2005/02/10

  Copyright (C) Ingo Eble,    ingoeble@web.de

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

#include "taylormodel_impl.h"
#include "factorialtable.h"
#include "polyeval.h"
#include <cmath>
#include <cstdlib>

namespace riot
{


/*-------------------------------------------------------------------------------------+

  Before we define the TaylorModel-Implementation, we have to introduce some
  helping objects (function objects):

  The following objects decide if an element will be inserted into the hash table when
  the key already exists or in all other cases if the value is adequate.
  See the functions 'insert' and 'insert_if' in file 'hashtable.h' for details.

  ---------------------------------------------------------------------------------------*/

  class IntervalAdditionNoValueCheck : public BasicAction<Monom,Interval,IntervalAdditionNoValueCheck>
  {
  private:

    /*
      To make it a Singleton.
    */

    IntervalAdditionNoValueCheck() {}
    ~IntervalAdditionNoValueCheck() {}
    IntervalAdditionNoValueCheck(const IntervalAdditionNoValueCheck&);
    IntervalAdditionNoValueCheck& operator = (const IntervalAdditionNoValueCheck&);

  public:

    typedef HashTableData<Monom,Interval> Node_data;

    static IntervalAdditionNoValueCheck& getObject()
      {
        static IntervalAdditionNoValueCheck singleton;
        return singleton;
      }

    /*
      The first operator is called when the key of the element which should be
      inserted already exists, the second operator is called from the
      'insert_if' function when inserting a new element. The function decides
      whether to do so or not.
    */

    bool operator () (Node_data& x, const Node_data& y)  { x.value() += y.value(); return true; }
    bool operator () (const Monom&, const Interval&   )  {                         return true; }
  };

  class RealAdditionValueCheck : public BasicAction<Monom,double,RealAdditionValueCheck>
  {
  private:

    /*
      To make it a Singleton.
    */

    RealAdditionValueCheck() : tol_(0.0) {}
    ~RealAdditionValueCheck() {}
    RealAdditionValueCheck(const RealAdditionValueCheck&);
    RealAdditionValueCheck& operator = (const RealAdditionValueCheck&);

  public:

    typedef HashTableData<Monom,double> Node_data;

    static RealAdditionValueCheck& getObject()
      {
        static RealAdditionValueCheck singleton;
        return singleton;
      }

    /*
      Set the tolerance for the value check.
    */

    void setTol(const double& d) { tol_ = d; }

    /*
      The first operator is called when the key of the element which should be
      inserted already exists, the second operator is called from the
      'insert_if' function when inserting a new element. The function decides
      whether to do so or not.
    */

    bool operator () (Node_data& x, const Node_data& y)
      {
        x.value() += y.value();
        return x.value() <= -tol_ || tol_ <= x.value();
      }
    bool operator () (const Monom&, const double&  val)
      {
        return val <= -tol_ || tol_ <= val;
      }

  private:

    double tol_;
  };

  class IntervalAdditionValueCheckWithContainer : public BasicAction<Monom,double,IntervalAdditionValueCheckWithContainer>
  {
  private:

    /*
      To make it a Singleton.
    */

    IntervalAdditionValueCheckWithContainer() {}
    ~IntervalAdditionValueCheckWithContainer() {}
    IntervalAdditionValueCheckWithContainer(const IntervalAdditionValueCheckWithContainer&);
    IntervalAdditionValueCheckWithContainer& operator = (const IntervalAdditionValueCheckWithContainer&);

  public:

    typedef HashTableData<Monom,double>            Node_data;
    typedef HashTable<Monom,Interval,MonomHasher2> ipolynomial_type;

    static IntervalAdditionValueCheckWithContainer& getObject()
      {
        static IntervalAdditionValueCheckWithContainer singleton;
        return singleton;
      }

    /*
      Set the pointer to the container.
    */

    void setPtr(ipolynomial_type *ptr) { ptr_ = ptr; }

    /*
      Set the tolerance for the value check.
    */

    void setTol(const double&       d) { tol_ = d;  }

    /*
      The first operator is called when the key of the element which should be
      inserted already exists, the second operator is called from the
      'insert_if' function when inserting a new element. The function decides
      whether to do so or not.
    */

    bool operator () (Node_data& x, const Node_data& y)
      {
        /*

          x * Monom + y * Monom \in ([x]+[y]) * Monom = mid([x]+[y]) * Monom + (([x]+[y]) - mid([x]+[y])) * Monom

          ----------+---------   ----------------+-----------------
          |                            |
          into Polynomial            into Restpolynomial

        */

        //Compute [x]+[y]
        Interval ix = x.value();
        ix += y.value();

        //Compute mid([x]+[y])
        double ixmid = ix.mid();

        if( -tol_ < ixmid && ixmid < tol_ )
        {
          ptr_->insert( x.key(), ix, IntervalAdditionNoValueCheck::getObject() );//All into restpolynomial.
          return false; //Tell the calling function to delete 'x' in the polynomial.
        }
        else
        {
          x.value() = ixmid; //Set result. (Into Polynomial, see above.)
          ptr_->insert( x.key(), ix - ixmid, IntervalAdditionNoValueCheck::getObject() ); //Into restpolynomial.
          return true; //Tell the calling function to not delete 'x'.
        }
      }

    bool operator () (const Monom& m, const double& val)
      {
        bool insert = ( val <= -tol_ || tol_ <= val );
        if( ! insert ) ptr_->insert( m, val, IntervalAdditionNoValueCheck::getObject() ); //Into restpolynomial.
        return insert;
      }

  private:

    /*
      Data.
    */

    ipolynomial_type *ptr_;
    double            tol_;
  };

  class ValueCheckOnly : public BasicAction<Monom,double,ValueCheckOnly>
  {
  private:

    /*
      To make it a Singleton.
    */

    ValueCheckOnly() {}
    ~ValueCheckOnly() {}
    ValueCheckOnly(const ValueCheckOnly&);
    ValueCheckOnly& operator = (const ValueCheckOnly&);

  public:

    typedef HashTableData<Monom,double> Node_data;

    static ValueCheckOnly& getObject()
      {
        static ValueCheckOnly singleton;
        return singleton;
      }

    /*
      Set the tolerance for the value check.
    */

    void setTol(const double& d) { tol_ = d; }

    /*
      The first operator is called when the key of the element which should be
      inserted already exists, the second operator is called from the
      'insert_if' function when inserting a new element. The function decides
      whether to do so or not.
    */

    bool operator () (Node_data& x, const Node_data& y) { return true; }
    bool operator () (const Monom&, const double&  val) { return val <= -tol_ || tol_ <= val; }

  private:

    /*
      Data.
    */

    double tol_;
  };

  class ValueCheckOnlyWithContainer : public BasicAction<Monom,double,ValueCheckOnlyWithContainer>
  {
  private:

    /*
      To make it a Singleton.
    */

    ValueCheckOnlyWithContainer() {}
    ~ValueCheckOnlyWithContainer() {}
    ValueCheckOnlyWithContainer(const ValueCheckOnlyWithContainer&);
    ValueCheckOnlyWithContainer& operator = (const ValueCheckOnlyWithContainer&);

  public:

    typedef HashTableData<Monom,double> Node_data;
    typedef HashTable<Monom,Interval,MonomHasher2> ipolynomial_type;

    static ValueCheckOnlyWithContainer& getObject()
      {
        static ValueCheckOnlyWithContainer singleton;
        return singleton;
      }

    /*
      Set the pointer to the container.
    */

    void setPtr(ipolynomial_type *ptr) { ptr_ = ptr; }

    /*
      Set the tolerance for the value check.
    */

    void setTol(const double&      d)  { tol_ =   d; }

    /*
      The first operator is called when the key of the element which should be
      inserted already exists, the second operator is called from the
      'insert_if' function when inserting a new element. The function decides
      whether to do so or not.
    */

    bool operator () (Node_data& x, const Node_data& y) { return true; }
    bool operator () (const Monom& m,const double& val)
      {
        bool insert = ( val <= -tol_ || tol_ <= val );
        if( !insert ) ptr_->insert( m, val, IntervalAdditionNoValueCheck::getObject() ); //Into restpolynomial.
        return insert;
      }

  private:

    ipolynomial_type *ptr_;
    double            tol_;
  };

/*-------------------------------------------------------------------------------+



  Now we can define the TaylorModel_Impl member functions.



  ---------------------------------------------------------------------------------*/

/*

  Initialization of static members with default values.

*/

  unsigned int           TaylorModel_Impl::order_         = 2;
  double                 TaylorModel_Impl::sparsity_tol_  = 1e-25;
  const PolynomialRange *TaylorModel_Impl::poly_eval_     = &IntervalEvaluation::getObject();
  const DegreeBase      *TaylorModel_Impl::degree_        = &Total_Degree::getObject();
  unsigned int           TaylorModel_Impl::Mode           = 0;
  Variables             *TaylorModel_Impl::vtab_          = &Variables::getTable();

/*

  Record size.

*/

  const unsigned int _RECORD_SIZE_ = 10000;

/*

  This function sets the chosen mode.

*/

  void TaylorModel_Impl::set_mode    (unsigned int m)
  {
    Mode = m;

    switch( m )
    {
    case NORMAL_MODE:
    {
      break;
    }
    case REMAINDER_REC:
    {
      //
      // Reset the record counters in each function and operation.
      //
      static unsigned int *ptr_1 = init_record_cntr(0); //Gets pointer to first element.
      for(unsigned int i = 0; i < 20; i++) ptr_1[i] = 0;
      break;
    }
    case REMAINDER_PLAY:
    {
      //
      // Reset the play counters in each function and operation.
      //
      static unsigned int *ptr_2 = init_play_cntr(0); //Gets pointer to first element.
      for(unsigned int i = 0; i < 20; i++) ptr_2[i] = 0;
      break;
    }
    case NO_REMAINDER:
    {
      break;
    }
    default:
    {
      std::cerr << "Mode not available." << std::endl;
    }
    }
  }

/*

  This function reacts on messages concerning the modes.

*/

  void TaylorModel_Impl::mode_message(unsigned int m)
  {
    switch( m )
    {
    case SEQUENCE_BREAK:
    {
      //
      // Reset the play counters in each function and operation.
      //
      static unsigned int *ptr_ = init_play_cntr(0); //Gets pointer to first element.
      for(unsigned int i = 0; i < 20; i++) ptr_[i] = 0;
      break;
    }
    default:
    {
      std::cerr << "Can't understand your message." << std::endl;
    }
    }
  }

/*

  The following functions create arrays of integer numbers. They
  are used as counters in play modus and record modus.

*/

  unsigned int* TaylorModel_Impl::init_play_cntr(unsigned int i)
  {
    static unsigned int *play_cntr_ = new unsigned int[20];//Executed only at first call.

    return &(play_cntr_[i]); //Return adress of element.
  }

  unsigned int* TaylorModel_Impl::init_record_cntr(unsigned int i)
  {
    static unsigned int *record_cntr_ = new unsigned int[20];//Executed only at first call.

    return &(record_cntr_[i]); //Return adress of element.
  }

/*

  Definition of the static member function 'Const_TM'.
  The function creates a constant polynomial with zero
  remainder interval. Before construction the
  constant value will be checked from the
  'ValueCheckOnly' object, maybe it is smaller than the
  'sparsity_tol_'.

*/

  TaylorModel_Impl* TaylorModel_Impl::Const_TM(const double& d)
  {
    TaylorModel_Impl *result = new TaylorModel_Impl();

    ValueCheckOnly::getObject().setTol( sparsity_tol_ );

    (result->polynomial_).insert_if( Monom(), d, ValueCheckOnly::getObject() );

    return result;
  }

/*

  Definition of private member functions.

*/

//
//
//
  bool TaylorModel_Impl::check_var_data(additional_var_data& x, const additional_var_data& y)
  {
    /*
      Check if the development points and the domains
      of equal variables in 'x' and 'y' are identical and if so
      copy the others to the left hand side.
      The arguments 'x' and 'y' are arrays with pointers to
      data. So we must check if these pointers are equal.
    */

    unsigned int x_size = x.size(), y_size = y.size();

    /*
      First we discuss the case when one operand is a Taylor-Model with a constant
      polynomial like ZERO_TM or ONE_TM (because they contain an empty 'additional_var_data' object).
    */

    if( x_size == 0 && y_size == 0 ) return true;
    if(                 y_size == 0 ) return true;
    if( x_size == 0                 )
    {
      x = y; //Copy data.
      return true;
    }

    /*
      Now there is data available.
    */

    std::vector<int> ind(0); //For indirect indexing.

    if( x_size < y_size )
    {
      bool equal = true;
      for(unsigned int i = 0; i < x_size; i++)
      {
        if( equal && (x[i].operator ->()) )
        {
          ind.push_back(i); //Save index with x[i].operator ->() != 0.
          if( y[i].operator ->() )
            equal &= ( (x[i].operator ->()) == (y[i].operator ->()) ); //Look if pointers are equal.
        }
      }
      if( equal ) //All pointers are equal.
      {
        additional_var_data tmp(y);
        for(unsigned int i = 0; i < ind.size(); i++) tmp[ ind[i] ] = x[ ind[i] ];
        x.swap(tmp);//Swap the handles.
      }
      return equal;
    }
    else
    {
      bool equal = true;
      for(unsigned int i = 0; i < y_size; i++)
      {
        if( equal && (y[i].operator ->()) )
        {
          ind.push_back(i); //Save index with y[i].operator ->() != 0.
          if( x[i].operator ->() )
            equal &= ( (x[i].operator ->()) == (y[i].operator ->()) ); //Look if pointers are equal.
        }
      }
      if( equal ) for(unsigned int i = 0; i < ind.size(); i++) x[ ind[i] ] = y[ ind[i] ];
      return equal;
    }
  }

//
// This function searches the monom 'm' in the polynomial
// and returns the coefficient.
//
  double TaylorModel_Impl::coefficient_of(const Monom& m) const
  {
    const_polynomial_iterator pos = polynomial_.find(m);
    if( pos != polynomial_.end() ) return (*pos).value();
    else                           return 0.0;
  }

//
// The following function searches all monomials with order
// 'n' and returns a list.
//
  std::vector<Monom> TaylorModel_Impl::monom_of_order(unsigned int n) const
  {
    std::vector<Monom> result;

    const_polynomial_iterator
      curr = polynomial_.begin(),
      last = polynomial_.end();

    while( curr != last )
    {
      if( ((*curr).key()).expnt_sum() == n ) result.push_back( (*curr).key() );
      ++curr;
    }

    return result;
  }

//
// The next function calculates an interval enclosure of the
// Taylor-Model.
//
  Interval TaylorModel_Impl::eval() const
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(19);
    static unsigned int *record_cntr = init_record_cntr(19);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    /*
      Play mode. Use recorded information for calculating the remainder
      and leave this function.
    */

    Interval enclosure;

    if( Mode & REMAINDER_PLAY )
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in eval: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      enclosure = record[ (*play_cntr)++ ];

      return enclosure + rest_;
    }

    /*
      In the normal or record mode we calculate the enclosure.
    */

    enclosure = (*poly_eval_)( &polynomial_, &data_ );

    if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all calculated data.
    {
      record[ (*record_cntr)++ ] = enclosure;

      if( *record_cntr > _RECORD_SIZE_ )
      {
        std::cout << "ERROR in eval: *** Record overrun. *** " << std::endl;
        std::exit(1);
      }
    }

    return enclosure + rest_;
  }

//
// This function reads the coefficient of the given monom 'm',
// returns it and deletes the monom from the polynomial.
//
  double TaylorModel_Impl::read_coefficient_and_remove(const Monom& m)
  {
    const_polynomial_iterator pos = polynomial_.find(m);
    if( pos != polynomial_.end() )
    {
      double coeff = (*pos).value();
      polynomial_.erase(pos);
      return coeff;
    }
    else return 0.0;
  }

//
// Implementation of the Taylor-Model-Arithmetic based on Hashtables:
//
//    - Negation;
//
  TaylorModel_Impl* TaylorModel_Impl::operator - () const
  {
    /*
      TM = P + I   =>  -TM = (-P) + (-I)
    */

    TaylorModel_Impl *result = new TaylorModel_Impl( *this ); //Call to copy-ctor.

    /*
      Negation of polynomial.
    */

    polynomial_iterator
      curr = (result->polynomial_).begin(),
      last = (result->polynomial_).end();

    while( curr != last ) { (*curr).value() = -(*curr).value(); ++curr; }

    /*
      Negation of interval.
    */

    result->rest_ = -(result->rest_);

    return result;
  }

//
//    - Addition with seizing rounding errors;
//
  TaylorModel_Impl& TaylorModel_Impl::add_with_seizing_rounding_errors(const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM += TM --" << std::endl;
#endif

    /*
      Data needed for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(0);
    static unsigned int *record_cntr = init_record_cntr(0);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    /*
      First we have to check if the development points and the domains
      of the concerned variables are identical. This is also important
      for the record feature e.g. if the left hand side is a constant Taylor
      model that contains no data. With the check of the variable data we make
      sure that after this operation the data of the variables is contained in
      the result.
    */

    if( check_var_data( data_ , s.data_ ) )
    {
      /*
        Play mode. Use recorded information for calculating the remainder
        and leave this function.
      */

      if( Mode & REMAINDER_PLAY )
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in TM += TM: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

        rest_ += s.rest_;
        rest_ += record[ (*play_cntr)++ ];

        return *this;
      }

      /*
        In the normal or record mode we calculate the polynomial and the remainder.

        TM_1 = P_1 + I_1, TM_2 = P_2 + I_2    =>   TM_1 + TM_2 = (P_1 + P_2) + (I_1 + I_2)
      */

      //
      // Polynomial with interval coefficients.
      //
      ipolynomial_type
        rest_poly( polynomial_.tab_size() + s.polynomial_.tab_size(), MonomHasher2::getFunction() );

      /*
        Add the polynomials:
        Iterate trough the polynomial of 's' and add every monom to the polynomial
        of 'this'.
      */

      //
      // Resize hash table if necessary.
      //
      if( polynomial_.tab_size() < polynomial_.number_of_entries() + s.polynomial_.number_of_entries() )
        polynomial_.resize( polynomial_.tab_size() + s.polynomial_.tab_size() );

      //
      // What to do when inserting an already existing monom?
      // The answer gives the object 'IntervalAdditionValueCheckWithContainer'.
      // Here we set the needed information.
      //
      IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly    );//Set the container.
      IntervalAdditionValueCheckWithContainer::getObject().setTol( sparsity_tol_ );//Set the tolerance.

      const_polynomial_iterator
        curr = (s.polynomial_).begin(),
        last = (s.polynomial_).end();

      while( curr != last )
      { //                      Monom        Coefficient                 (see above)
        polynomial_.insert( (*curr).key(), (*curr).value(), IntervalAdditionValueCheckWithContainer::getObject() );
        ++curr;
      }

      /*
        Add the intervals.
      */

      rest_ += s.rest_;

      /*
        Bound the restpolynomial.
      */

      Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

      if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder in play mode.
      {
        record[ (*record_cntr)++ ] = B_REST;

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in TM += TM: *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }

      rest_ += B_REST;

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM += TM --" << std::endl;
#endif
    }
    else std::cerr << "Data of concerning variables does not correspond." << std::endl;

    return *this;
  }

//
//    - Addition without seizing rounding errors;
//
  TaylorModel_Impl& TaylorModel_Impl::add_without_seizing_rounding_errors(const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM += TM --" << std::endl;
#endif

    /*
      First we have to check if the development points and the domains
      of the concerned variables are identical.
    */

    if( check_var_data( data_ , s.data_ ) )
    {
      /*
        Only add the polynomials without considering
        the remainder term.

        TM_1 = P_1 + I_1, TM_2 = P_2 + I_2    =>   TM_1 + TM_2 = P_1 + P_2

        Iterate trough the polynomial of 's' and add every monom to the polynomial
        of 'this'.
      */

      //
      // Resize hash table if necessary.
      //
      if( polynomial_.tab_size() < polynomial_.number_of_entries() + s.polynomial_.number_of_entries() )
        polynomial_.resize( polynomial_.tab_size() + s.polynomial_.tab_size() );

      //
      // What to do when inserting an existing monom?
      // The answer gives the object 'RealAdditionValueCheck'.
      //
      RealAdditionValueCheck::getObject().setTol( sparsity_tol_ );

      const_polynomial_iterator
        curr = (s.polynomial_).begin(),
        last = (s.polynomial_).end();

      while( curr != last )
      { //                      Monom        Coefficient                 (see above)
        polynomial_.insert( (*curr).key(), (*curr).value(), RealAdditionValueCheck::getObject() );
        ++curr;
      }

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM += TM --" << std::endl;
#endif
    }
    else std::cerr << "Data of concerning variables does not correspond." << std::endl;

    return *this;
  }

//
//    - Subtraction with seizing rounding errors;
//
  TaylorModel_Impl& TaylorModel_Impl::subtract_with_seizing_rounding_errors(const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM -= TM --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(1);
    static unsigned int *record_cntr = init_record_cntr(1);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    /*
      First we have to check if the development points and the domains
      of the concerned variables are identical. This is also important
      for the record feature e.g. if the left hand side is a constant Taylor
      model that contains no data. With the check of the variable data we make
      sure that after this operation the data of the variables is contained in
      the result.
    */

    if( check_var_data( data_ , s.data_ ) )
    {
      /*
        Play mode. Use recorded information for calculating the remainder
        and leave this function.
      */

      if( Mode & REMAINDER_PLAY )
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in TM -= TM: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

        rest_ -= s.rest_;
        rest_ += record[ (*play_cntr)++ ];

        return *this;
      }

      /*
        In the normal or record mode we calculate the polynomial and the remainder.

        TM_1 = P_1 + I_1, TM_2 = P_2 + I_2    =>   TM_1 - TM_2 = (P_1 - P_2) + (I_1 - I_2)
      */

      //
      // Polynomial with interval coefficients.
      //
      ipolynomial_type rest_poly( polynomial_.tab_size() + s.polynomial_.tab_size(), MonomHasher2::getFunction() );

      /*
        Subtract the polynomials:
        Iterate trough the polynomial of 's' and subtract every monom from the
        polynomial of 'this'.
      */

      //
      // Resize hash table if necassary.
      //
      if( polynomial_.tab_size() < polynomial_.number_of_entries() + s.polynomial_.number_of_entries() )
        polynomial_.resize( polynomial_.tab_size() + s.polynomial_.tab_size() );

      //
      // What to do when inserting an existing monom?
      // The answer gives the object 'IntervalAdditionValueCheckWithContainer'.
      // Here we set the needed information.
      //
      IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly    );
      IntervalAdditionValueCheckWithContainer::getObject().setTol( sparsity_tol_ );

      const_polynomial_iterator
        curr = (s.polynomial_).begin(),
        last = (s.polynomial_).end();

      while( curr != last )//                +------- => Subtraction here -------+
      {                  //                |                                   |
        //                v                                   v
        polynomial_.insert( (*curr).key(), -(*curr).value(), IntervalAdditionValueCheckWithContainer::getObject() );
        ++curr;
      }

      /*
        Subtract the intervals.
      */

      rest_ -= s.rest_;

      /*
        Bound the restpolynomial.
      */

      Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

      if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = B_REST;

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in TM -= TM: *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }

      rest_ += B_REST;

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM -= TM --" << std::endl;
#endif
    }
    else std::cerr << "Data of concerning variables does not correspond." << std::endl;

    return *this;
  }

//
//    - Addition without seizing rounding errors;
//
  TaylorModel_Impl& TaylorModel_Impl::subtract_without_seizing_rounding_errors(const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM -= TM --" << std::endl;
#endif

    /*
      First we have to check if the development points and the domains
      of the concerned variables are identical.
    */

    if( check_var_data( data_ , s.data_ ) )
    {
      /*
        Only subtract the polynomials without considering
        the remainder term.

        TM_1 = P_1 + I_1, TM_2 = P_2 + I_2    =>    TM_1 - TM_2 = P_1 - P_2

        Iterate trough the polynomial of 's' and subtract every monom from the
        polynomial of 'this'.
      */

      //
      // Resize hash table if necassary.
      //
      if( polynomial_.tab_size() < polynomial_.number_of_entries() + s.polynomial_.number_of_entries() )
        polynomial_.resize( polynomial_.tab_size() + s.polynomial_.tab_size() );

      //
      // What to do when inserting an existing monom?
      // The answer gives the object 'RealAdditionValueCheck'.
      //
      RealAdditionValueCheck::getObject().setTol( sparsity_tol_ );

      const_polynomial_iterator
        curr = (s.polynomial_).begin(),
        last = (s.polynomial_).end();

      while( curr != last )//                +------- => Subtraction here -------+
      {                  //                |                                   |
        //                v                                   v
        polynomial_.insert( (*curr).key(), -(*curr).value(), RealAdditionValueCheck::getObject() );
        ++curr;
      }

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM -= TM --" << std::endl;
#endif
    }
    else std::cerr << "Data of concerning variables does not correspond." << std::endl;

    return *this;
  }

//
//    - Multiplikation with seizing rounding errors;
//
  TaylorModel_Impl& TaylorModel_Impl::multiply_with_seizing_rounding_errors(const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM *= TM --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(2);
    static unsigned int *record_cntr = init_record_cntr(2);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    /*
      First we have to check if the development points and the domains
      of the concerned variables are identical. This is also important
      for the record feature e.g. if the left hand side is a constant Taylor
      model that contains no data. With the check of the variable data we make
      sure that after this operation the data of the variables is contained in
      the result.
    */

    if( check_var_data( data_ , s.data_ ) )
    {
      /*
        Play mode. Use recorded information for calculating the remainder
        and leave this function.
      */

      if( Mode & REMAINDER_PLAY )
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in TM *= TM: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

        Interval
          BP = record[ (*play_cntr)++ ],
          BQ = record[ (*play_cntr)++ ],
          I  = rest_                   ,
          J  = s.rest_                 ,
          BR = record[ (*play_cntr)++ ];

        rest_  = (BP*J + I*(BQ+J)) & (BQ*I + J*(BP+I));
        rest_ += BR;

        return *this;
      }

      /*
        In the normal or record mode we calculate the polynomial and the remainder.

        TM_1 = P_1 + I_1, TM_2 = P_2 + I_2 => TM_1 * TM_2 = P_1*P_2 + (B(P_1)*I_2 + I_1*(B(P_2) + I_2))
      */

      /*
        Multiply the polynomials:
        Iterate trough the polynomial of 's' and multiply every monom with the
        polynomial of 'this'.
      */

      unsigned int  p_size =   polynomial_.tab_size();
      unsigned int sp_size = s.polynomial_.tab_size();

      //
      // Polynomial with interval coefficients.
      //
      ipolynomial_type iprod_poly( (p_size > sp_size)?p_size:sp_size, MonomHasher2::getFunction() );

      //
      // Polynomial with interval coefficients. Contains round off and cut off errors.
      //
      ipolynomial_type rest_poly( p_size * sp_size, MonomHasher2::getFunction() );

      const_polynomial_iterator
        rhs_curr,
        rhs_last = (s.polynomial_).end(),
        lhs_curr = polynomial_.begin(),
        lhs_last = polynomial_.end();

      /*
        x * Monom1 * y * Monom2 \in [x]*[y] * Monom3

        = mid([x]*[y]) * Monom3 + ([x]*[y] - mid([x]*[y])) * Monom3

        -----------+---------   ----------------+-----------------
        |                            |
        into Polynomial            into Restpolynomial
      */

      polynomial_type prod_poly( p_size * sp_size, MonomHasher2::getFunction() );

      while( lhs_curr != lhs_last )
      {
        Monom  lhs_m = (*lhs_curr).key();
        double lhs_v = (*lhs_curr).value();

        rhs_curr = (s.polynomial_).begin();

        while( rhs_curr != rhs_last )
        {
          //Product of monomials (=Monom3, see above).
          Monom rhs_m = (*rhs_curr).key();
          rhs_m *= lhs_m;

          //Product of coefficients (\in [x]*[y]).
          Interval icoeff = (*rhs_curr).value();
          icoeff *= lhs_v;

          //Check degree of monomial rhs_m (=Monom3).
          if( degree_->check( &rhs_m, order_ ) )
          {
            iprod_poly.insert( rhs_m, icoeff, IntervalAdditionNoValueCheck::getObject() );
          }
          else
          {
            rest_poly .insert( rhs_m, icoeff, IntervalAdditionNoValueCheck::getObject() );
          }

          ++rhs_curr;
        }

        ++lhs_curr;
      }

      const_ipolynomial_iterator
        curr = iprod_poly.begin(),
        last = iprod_poly.end();

      while( curr != last )
      {
        //Compute mid([x]*[y]).
        double mid_icoeff = (curr->value()).mid();

        if( mid_icoeff <= -sparsity_tol_ || sparsity_tol_ <= mid_icoeff )
        {
          prod_poly.insert( curr->key(), mid_icoeff                , DefaultAction<Monom,double> ::getObject() );
          rest_poly.insert( curr->key(), curr->value() - mid_icoeff, IntervalAdditionNoValueCheck::getObject() );
        }
        else
        {
          rest_poly.insert( curr->key(), curr->value()             , IntervalAdditionNoValueCheck::getObject() );
        }

        ++curr;
      }

      /*
        Compute the interval term.
      */

      //
      // Since TM_1 * TM_2 = TM_2 * TM_1 we should compute the intersection of
      // B(P_1)*I_2 + I_1*(B(P_2) + I_2) with B(P_2)*I_1 + I_2*(B(P_1) + I_1) for an improved enclosure.
      //
      Interval
        BP = (*poly_eval_)( &  polynomial_, &  data_ ),
        BQ = (*poly_eval_)( &s.polynomial_, &s.data_ ),
        I  =   rest_,
        J  = s.rest_;

      if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = BP;
        record[ (*record_cntr)++ ] = BQ;

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in TM *= TM: *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }

      rest_ = (BP*J + I*(BQ+J)) & (BQ*I + J*(BP+I));

      /*
        Bound the restpolynomial (which includes cut off and round off errors!)
        and add the result to I.
      */

      Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

      if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = B_REST;

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in TM *= TM: *** Remainder record overrun. *** " << std::endl;
          std::exit(1);
        }
      }

      rest_ += B_REST;

      polynomial_ = prod_poly;

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM *= TM --" << std::endl;
#endif
    }
    else std::cerr << "Data of concerning variables does not correspond." << std::endl;

    return *this;
  }

//
//    - Multiplication without seizing rounding errors;
//
  TaylorModel_Impl& TaylorModel_Impl::multiply_without_seizing_rounding_errors(const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM *= TM --" << std::endl;
#endif

    /*
      First we have to check if the development points and the domains
      of the concerned variables are identical.
    */

    if( check_var_data( data_ , s.data_ ) )
    {
      /*
        Only multiply the polynomials without considering
        the remainder term.

        TM_1 = P_1 + I_1, TM_2 = P_2 + I_2  => TM_1 * TM_2 = P_1*P_2

        Iterate trough the polynomial of 's' and multiply every monom with the
        polynomial of 'this'.
      */

      unsigned int  p_size = polynomial_.tab_size();
      unsigned int sp_size = s.polynomial_.tab_size();

      polynomial_type prod_poly( p_size * sp_size, MonomHasher2::getFunction() );

      //
      // When inserting already existing monoms we must check if the sum of the
      // coefficients has the correct size (see below).
      //
      RealAdditionValueCheck::getObject().setTol( sparsity_tol_ );

      const_polynomial_iterator
        rhs_curr,
        rhs_last = (s.polynomial_).end(),
        lhs_curr = polynomial_.begin(),
        lhs_last = polynomial_.end();

      while( lhs_curr != lhs_last )
      {
        Monom  lhs_m = (*lhs_curr).key();
        double lhs_v = (*lhs_curr).value();

        rhs_curr = (s.polynomial_).begin();

        while( rhs_curr != rhs_last )
        {
          //Product of monomials.
          Monom rhs_m = (*rhs_curr).key();
          rhs_m *= lhs_m;

          //Check degree of monomial 'rhs_m'.
          if( degree_->check( &rhs_m, order_ ) )
          {
            //Product of coefficients.
            double coeff = (*rhs_curr).value();
            coeff *= lhs_v;

            if( coeff <= -sparsity_tol_ || sparsity_tol_ <= coeff )
            {
              //
              // We call the 'insert_if' function, because it can be that 'rhs_m' doesn't
              // exist in 'prod_poly'. In this case we must check if the coefficient has the
              // correct size. The 'insert_if' function calls the
              // corresponding operator to check this.
              //
              prod_poly.insert_if( rhs_m, coeff, RealAdditionValueCheck::getObject() );
            }
          }

          ++rhs_curr;
        }

        ++lhs_curr;
      }

      polynomial_ = prod_poly;

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM *= TM --" << std::endl;
#endif
    }
    else std::cerr << "Data of concerning variables does not correspond." << std::endl;

    return *this;
  }

//
//    - Division with and without seizing rounding errors;
//
  TaylorModel_Impl& TaylorModel_Impl::operator /= (const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM /= TM --" << std::endl;
#endif

    TaylorModel_Impl *inverse_of_s = invert(s);

    (*this) *= (*inverse_of_s);

    delete inverse_of_s;

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM /= TM --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl* invert (const TaylorModel_Impl& s)
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(3);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(3);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    TaylorModel_Impl *result = 0;

    Interval B = (*TaylorModel_Impl::poly_eval_)( &s.polynomial_, &(s.data_) ) + s.rest_;

    if( ! B.contains(0.0) || (TaylorModel_Impl::Mode & REMAINDER_PLAY) ) //Zero in 's'?
    {
      /*
        Separate the polynomial of 's' into the constant and the non-constant part.
        At this point there exist a constant part, otherwise the Taylor Model 's'
        contains zero. But this was checked above.
      */

      double constant;
      TaylorModel_Impl *s_without_cnst = 0; //Taylor Model without constant part in polynomial.

      if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in invert(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        constant = record[ (*play_cntr)++ ].inf(); //The constant is saved as type interval. So we
        //have to use only one record array.

        result         = new TaylorModel_Impl( 1, Adapt::ZERO_INTERVAL(), s.data_);
        s_without_cnst = new TaylorModel_Impl( 1, s.rest_               , s.data_);
      }
      else //Other mode.
      {
        TaylorModel_Impl::const_polynomial_iterator
          pos = (s.polynomial_).find( Monom() ); //Find position of constant part ...
        constant = (*pos).value(); //... and read the constant.

        result         = new TaylorModel_Impl( (s.polynomial_).number_of_entries(), Adapt::ZERO_INTERVAL(), s.data_);
        s_without_cnst = new TaylorModel_Impl(s);       //Copy 's' ...
        (s_without_cnst->polynomial_).erase( Monom() ); //and remove constant part.

        if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
        {
          record[ (*record_cntr)++ ] = constant; //We save the constant as type interval.

          if( *record_cntr > _RECORD_SIZE_ )
          {
            std::cout << "ERROR in invert(TM): *** Record overrun. *** " << std::endl;
            std::exit(1);
          }
        }
      }

      /*
        Now compute the polynomial part of '1/s' with the horner scheme.
        The polynomial part is:

        1/constant * { 1 - s/constant + (s/constant)^2 - ... + (-1)^n * (s/constant)^n }

      */

      (*s_without_cnst) /= constant; // = s/constant.

      double coeff;
      if( (TaylorModel_Impl::order_ & 1) != 0 )
        coeff = -1.0; // = (-1)^n for odd n.
      else
        coeff =  1.0; // = (-1)^n for even n.

      (*result) += coeff;
      for(unsigned int i = TaylorModel_Impl::order_; i > 0; i-- )
      {
        (*result) *= (*s_without_cnst);
        coeff     *= -1.0;
        (*result) += coeff;
      }
      (*result) /= constant;

      /*
        Then we compute the remainder term. It is:

        1/constant * (-1)^(n+1) *  (s/constant)^(n+1) * 1/( 1 + theta * s/constant )^(n+2),  theta \in (0,1).

      */

      if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
      {
        Interval BS;
        if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
        {
          if( *record_cntr == 0 )
          {
            std::cout << "ERROR in invert(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
            std::exit(1);
          }

          BS = record[ (*play_cntr)++ ] + s_without_cnst->rest_;

          if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.
        }
        else //Other mode.
        {
          Interval
            bound_s_without_cnst = (*TaylorModel_Impl::poly_eval_)( &(s_without_cnst->polynomial_), &(s_without_cnst->data_) );

          BS = bound_s_without_cnst + s_without_cnst->rest_;

          if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
          {
            record[ (*record_cntr)++ ] = bound_s_without_cnst;

            if( *record_cntr > _RECORD_SIZE_ )
            {
              std::cout << "ERROR in invert(TM): *** Record overrun. *** " << std::endl;
              std::exit(1);
            }
          }
        }

        result->rest_ += power(-BS,TaylorModel_Impl::order_+1)
          * exp( -int(TaylorModel_Impl::order_+2) * ln(1.0 + Interval(0,1) * BS) )
          / constant;
      }

      //... that was it.

      delete s_without_cnst; //Delete copy of 's'.
    }
    else std::cerr << "invert(TM): Division by zero." << std::endl;

    return result;
  }

//
// Now there follow the arithmetic operations with type
// Interval and double.
//
  TaylorModel_Impl& TaylorModel_Impl::add_with_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM += IV --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(4);
    static unsigned int *record_cntr = init_record_cntr(4);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    if( Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in TM += IV: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      rest_ += (s - s.mid());
      rest_ += record[ (*play_cntr)++ ];

      return *this;
    }

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).
    */

    ipolynomial_type rest_poly(2, MonomHasher2::getFunction());

    /*
      Add the polynomials.
    */

    double mid_s = s.mid();

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'IntervalAdditionValueCheckWithContainer'.
    // Here we set the needed information.
    //
    IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly    );
    IntervalAdditionValueCheckWithContainer::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), mid_s, IntervalAdditionValueCheckWithContainer::getObject() );

    /*
      Add the intervals.
    */

    rest_ += (s - mid_s);

    /*
      Bound the restpolynomial (which includes round off and cut off errors!)
      and add the result to I.
    */

    Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

    if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
    {
      record[ (*record_cntr)++ ] = B_REST;

      if( *record_cntr > _RECORD_SIZE_ )
      {
        std::cout << "ERROR in TM += IV: *** Record overrun. *** " << std::endl;
        std::exit(1);
      }
    }

    rest_ += B_REST;

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM += IV --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::add_without_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM += IV --" << std::endl;
#endif

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).

      Only add the polynomials without considering
      the remainder term.
    */

    double mid_s = s.mid();

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'RealAdditionValueCheck'.
    //
    RealAdditionValueCheck::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), mid_s, RealAdditionValueCheck::getObject() );

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM += IV --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::subtract_with_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM -= IV --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(5);
    static unsigned int *record_cntr = init_record_cntr(5);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    if( Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in TM -= IV: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      rest_ -= (s - s.mid());
      rest_ += record[ (*play_cntr)++ ];

      return *this;
    }

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).
    */

    ipolynomial_type rest_poly(2, MonomHasher2::getFunction());

    /*
      Add the polynomials.
    */

    double mid_s = s.mid();

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'IntervalAdditionValueCheckWithContainer'.
    // Here we set the needed information.
    //
    IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly    );
    IntervalAdditionValueCheckWithContainer::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), -mid_s, IntervalAdditionValueCheckWithContainer::getObject() );

    /*
      Subtract the intervals.
    */

    rest_ -= (s - mid_s);

    /*
      Bound the restpolynomial (which includes round off and cut off errors!)
      and add the result to I.
    */

    Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

    if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
    {
      record[ (*record_cntr)++ ] = B_REST;

      if( *record_cntr > _RECORD_SIZE_ )
      {
        std::cout << "ERROR in TM -= IV: *** Record overrun. *** " << std::endl;
        std::exit(1);
      }
    }

    rest_ += B_REST;

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM -= IV --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::subtract_without_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM -= IV --" << std::endl;
#endif

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).
    */

    double mid_s = s.mid();

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'RealAdditionValueCheck'.
    // Here we set the needed information.
    //
    RealAdditionValueCheck::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), -mid_s, RealAdditionValueCheck::getObject() );

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM -= IV --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::multiply_with_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM *= IV --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(6);
    static unsigned int *record_cntr = init_record_cntr(6);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    if( Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in TM *= IV: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      Interval
        BP = record[ (*play_cntr)++ ],
        I  = rest_,
        sm = s-s.mid();

      rest_  = (BP*sm + I*s) & (s.mid()*I + sm*(BP+I));
      rest_ += record[ (*play_cntr)++ ];

      return *this;
    }

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).
      We use the algorithm of the TM *= TM operation with the polynomial
      multiplication of the TM *= double operation.
    */

    ipolynomial_type rest_poly( polynomial_.tab_size(), MonomHasher2::getFunction() );

    /*
      Multiply polynomial with mid(s).
    */

    polynomial_iterator
      prev = polynomial_.begin(),
      curr = prev,
      last = polynomial_.end();

    double s_mid = s.mid();
    bool previous_element_exist = false;
    while( curr != last )
    {
      Interval ival = (*curr).value();
      ival *= s_mid;

      double ival_mid = ival.mid();

      if( -sparsity_tol_ < ival_mid && ival_mid < sparsity_tol_ )
      {
        rest_poly.insert( (*curr).key(), ival, IntervalAdditionNoValueCheck::getObject() );
        polynomial_.erase( curr );

        if( previous_element_exist ) { curr = prev; ++curr; }
        else                           curr = polynomial_.begin();
      }
      else
      {
        (*curr).value() = ival_mid;
        rest_poly.insert( (*curr).key(), ival-ival_mid, IntervalAdditionNoValueCheck::getObject() );

        prev = curr;
        ++curr;
        previous_element_exist = true;
      }
    }

    /*
      Compute the interval term.
    */

    //
    // We use:
    // (P+I)*(mid(s)+(s-mid(s))) = P*mid(s) + B(P)*(s-mid(s)) + I*s = P*mid(s) + mid(s)*I + (s-mid(s))*(B(P)+I).
    //
    Interval
      BP = (*poly_eval_)( &polynomial_, &data_ ),
      I  = rest_,
      sm = s-s_mid;

    rest_ = (BP*sm + I*s) & (s_mid*I + sm*(BP+I));

    /*
      Bound the restpolynomial (which includes round off and cut off errors!)
      and add the result to I.
    */

    Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

    if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
    {
      record[ (*record_cntr)++ ] = BP;
      record[ (*record_cntr)++ ] = B_REST;

      if( *record_cntr > _RECORD_SIZE_ )
      {
        std::cout << "ERROR in TM *= IV: *** Record overrun. *** " << std::endl;
        std::exit(1);
      }
    }

    rest_ += B_REST;

    /*
      Adapt the size of the polynomial if necessary.
    */

    polynomial_.adapt_tab_size();

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM *= IV --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::multiply_without_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM *= IV --" << std::endl;
#endif

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).
      We use the algorithm of the TM *= TM operation with the polynomial
      multiplication of the TM *= double operation.
    */

    polynomial_iterator
      prev = polynomial_.begin(),
      curr = prev,
      last = polynomial_.end();

    double s_mid = s.mid();
    bool previous_element_exist = false;
    while( curr != last )
    {
      double val = (*curr).value();
      val *= s_mid;

      if( -sparsity_tol_ < val && val < sparsity_tol_ )
      {
        polynomial_.erase( curr );

        if( previous_element_exist ) { curr = prev; ++curr; }
        else                           curr = polynomial_.begin();
      }
      else
      {
        (*curr).value() = val;

        prev = curr;
        ++curr;
        previous_element_exist = true;
      }
    }

    /*
      Adapt the size of the polynomial if necessary.
    */

    polynomial_.adapt_tab_size();

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM *= IV --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::divide_with_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM /= IV --" << std::endl;
#endif

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).
    */

    if( !s.contains(0.0) ) //Zero in 's'?
    {
      TaylorModel_Impl *inverse_of_s = new TaylorModel_Impl( polynomial_.number_of_entries(), Adapt::ZERO_INTERVAL(), data_);

      /*
        Separate the polynomial of 's' into the constant and the non-constant part.
        At this point there exist a constant part, otherwise the Taylor Model 's'
        contains zero. But this was checked above.
      */

      double   constant      = s.mid(); //The constant part from 's'.
      Interval interval_part = s - constant; //"Taylor Model" without constant part in polynomial.

      /*
        Now compute the polynomial part of '1/s' with the horner scheme.
        The polynomial part is:

        1/constant * { 1 - s/constant + (s/constant)^2 - ... + (-1)^n * (s/constant)^n }

      */

      interval_part /= constant; // = s/constant.

      double coeff;
      if( (TaylorModel_Impl::order_ & 1) != 0 )
        coeff = -1.0; // = (-1)^n for odd n.
      else
        coeff =  1.0; // = (-1)^n for even n.

      inverse_of_s->rest_ += coeff;
      for(unsigned int i = TaylorModel_Impl::order_; i > 1; i-- )
      {
        inverse_of_s->rest_ *= interval_part;
        coeff     *= -1.0;
        inverse_of_s->rest_ += coeff;
      }
      inverse_of_s->rest_ /= constant;

      (*inverse_of_s) += 1.0; //Polynomial is 1/constant.
      (*inverse_of_s) /= constant;

      /*
        Then we compute the remainder term. It is:

        1/constant * (-1)^(n+1) *  (s/constant)^(n+1) * 1/( 1 + theta * s/constant )^(n+2),  theta \in (0,1).

      */

      Interval bound_s_without_cnst = interval_part; //This is a bound for 's/constant'.

      inverse_of_s->rest_ += power(-bound_s_without_cnst,TaylorModel_Impl::order_+1)
        / (constant * power(1.0 + Interval(0,1) * bound_s_without_cnst,TaylorModel_Impl::order_+2) );

      //... that was it.

      (*this) *= (*inverse_of_s); //'inverse_of_s' contains '1/s'.

      delete inverse_of_s;

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM /= IV --" << std::endl;
#endif

      return *this;
    }
    else std::cerr << "TM /= Interval: Division by zero." << std::endl;

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::divide_without_seizing_rounding_errors(const Interval& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM /= IV --" << std::endl;
#endif

    /*
      The interval 's' is equal to the Taylor Model mid(s) + (s-mid(s)).
    */

    if( !s.contains(0.0) ) //Zero in 's'?
    {
      TaylorModel_Impl *inverse_of_s = new TaylorModel_Impl( polynomial_.number_of_entries(), Adapt::ZERO_INTERVAL(), data_);

      /*
        Now compute the polynomial part of '1/s' with the horner scheme.
        The polynomial part is:

        1/constant * { 1 - s/constant + (s/constant)^2 - ... + (-1)^n * (s/constant)^n }

        But s is zero, so the polynomial part is 1/constant.
      */

      (*inverse_of_s) += 1.0;
      (*inverse_of_s) /= s.mid();;

      (*this) *= (*inverse_of_s); //'inverse_of_s' contains '1/s'.

      delete inverse_of_s;

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM /= IV --" << std::endl;
#endif

      return *this;
    }
    else std::cerr << "TM /= Interval: Division by zero." << std::endl;

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::add_with_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM += d --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(7);
    static unsigned int *record_cntr = init_record_cntr(7);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    if( Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in TM += double: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      rest_ += record[ (*play_cntr)++ ];

      return *this;
    }


    /*
      The value 'd' is equal to the Taylor Model d + [0,0].
    */

    ipolynomial_type rest_poly(2, MonomHasher2::getFunction());

    /*
      Add the polynomials.
    */

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'IntervalAdditionValueCheckWithContainer'.
    // Here we set the needed information.
    //
    IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly    );
    IntervalAdditionValueCheckWithContainer::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), d, IntervalAdditionValueCheckWithContainer::getObject() );

    /*
      Bound the restpolynomial (which includes round off and cut off errors!)
      and add the result to I.
    */

    Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

    if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
    {
      record[ (*record_cntr)++ ] = B_REST;

      if( *record_cntr > _RECORD_SIZE_ )
      {
        std::cout << "ERROR in TM += double: *** Record overrun. *** " << std::endl;
        std::exit(1);
      }
    }

    rest_ += B_REST;

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM += d --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::add_without_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM += d --" << std::endl;
#endif

    /*
      The value 'd' is equal to the Taylor Model d + [0,0].

      Only add the polynomials without considering
      the remainder term.
    */

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'RealAdditionValueCheck'.
    // Here we set the needed information.
    //
    RealAdditionValueCheck::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), d, RealAdditionValueCheck::getObject() );

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM += d --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::subtract_with_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM -= d --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(8);
    static unsigned int *record_cntr = init_record_cntr(8);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    if( Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in TM -= double: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      rest_ += record[ (*play_cntr)++ ];

      return *this;
    }

    /*
      The value 'd' is equal to the Taylor Model d + [0,0].
    */

    ipolynomial_type rest_poly(2, MonomHasher2::getFunction());

    /*
      Subtract the polynomials.
    */

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'IntervalAdditionValueCheckWithContainer'.
    // Here we set the needed information.
    //
    IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly    );
    IntervalAdditionValueCheckWithContainer::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), -d, IntervalAdditionValueCheckWithContainer::getObject() );

    /*
      Bound the restpolynomial (which includes round off and cut off errors!)
      and add the result to I.
    */

    Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

    if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
    {
      record[ (*record_cntr)++ ] = B_REST;

      if( *record_cntr > _RECORD_SIZE_ )
      {
        std::cout << "ERROR in TM -= double: *** Record overrun. *** " << std::endl;
        std::exit(1);
      }
    }

    rest_ += B_REST;

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM -= d --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::subtract_without_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM -= d --" << std::endl;
#endif

    /*
      The value 'd' is equal to the Taylor Model d + [0,0].

      Subtract the polynomials without seizing the rounding errors.
      The remainder term is not considered.
    */

    //
    // What to do when inserting an existing monom?
    // The answer gives the object 'RealAdditionValueCheck'.
    // Here we set the needed information.
    //
    RealAdditionValueCheck::getObject().setTol( sparsity_tol_ );

    polynomial_.insert_if( Monom(), -d, RealAdditionValueCheck::getObject() );

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM -= d --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::multiply_with_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM *= d --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(9);
    static unsigned int *record_cntr = init_record_cntr(9);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    if( Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in TM *= double: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      rest_ *= d;
      rest_ += record[ (*play_cntr)++ ];

      return *this;
    }

    /*
      The value 'd' is equal to the Taylor Model d + [0,0].
    */

    ipolynomial_type rest_poly( polynomial_.tab_size(), MonomHasher2::getFunction() );

    /*
      Multiply polynomial with 'd'.
    */

    polynomial_iterator
      prev = polynomial_.begin(),
      curr = prev,
      last = polynomial_.end();

    bool previous_element_exist = false;
    while( curr != last )
    {
      Interval ival = (*curr).value();
      ival *= d;

      double ival_mid = ival.mid();

      if( -sparsity_tol_ < ival_mid && ival_mid < sparsity_tol_ )
      {
        rest_poly.insert( (*curr).key(), ival, IntervalAdditionNoValueCheck::getObject() );
        polynomial_.erase( curr );

        if( previous_element_exist ) { curr = prev; ++curr; }
        else                           curr = polynomial_.begin();
      }
      else
      {
        (*curr).value() = ival_mid;
        rest_poly.insert( (*curr).key(), ival-ival_mid, IntervalAdditionNoValueCheck::getObject() );

        prev = curr;
        ++curr;
        previous_element_exist = true;
      }
    }

    /*
      Multiply interval term.
    */

    rest_ *= d;

    /*
      Bound the restpolynomial (which includes round off and cut off errors!)
      and add the result to I.
    */

    Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

    if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
    {
      record[ (*record_cntr)++ ] = B_REST;

      if( *record_cntr > _RECORD_SIZE_ )
      {
        std::cout << "ERROR in TM *= double: *** Record overrun. *** " << std::endl;
        std::exit(1);
      }
    }

    rest_ += B_REST;

    /*
      Adapt the size of the polynomial if necessary.
    */

    polynomial_.adapt_tab_size();

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM *= d --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::multiply_without_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM *= d --" << std::endl;
#endif

    /*
      The value 'd' is equal to the Taylor Model d + [0,0].

      Multiply polynomial with 'd' without considering
      the remainder term.
    */

    polynomial_iterator
      prev = polynomial_.begin(),
      curr = prev,
      last = polynomial_.end();

    bool previous_element_exist = false;
    while( curr != last )
    {
      double val = (*curr).value();
      val *= d;

      if( -sparsity_tol_ < val && val < sparsity_tol_ )
      {
        polynomial_.erase( curr );
        if( previous_element_exist ) { curr = prev; ++curr; }
        else                           curr = polynomial_.begin();
      }
      else
      {
        (*curr).value() = val;
        prev = curr;
        ++curr;
        previous_element_exist = true;
      }
    }

    polynomial_.adapt_tab_size();

#ifdef HASHINFO
    polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- TM *= d --" << std::endl;
#endif

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::divide_with_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM /= d --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(10);
    static unsigned int *record_cntr = init_record_cntr(10);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    if( Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in TM /= double: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      rest_ /= d;
      rest_ += record[ (*play_cntr)++ ];

      return *this;
    }

    /*
      The value 'd' is equal to the Taylor Model d + [0,0].
    */

    if( d != 0.0 )
    {
      ipolynomial_type rest_poly( polynomial_.tab_size(), MonomHasher2::getFunction() );

      /*
        Divide polynomial with 'd'.
      */

      polynomial_iterator
        prev = polynomial_.begin(),
        curr = prev,
        last = polynomial_.end();

      bool previous_element_exist = false;
      while( curr != last )
      {
        Interval ival = (*curr).value();
        ival /= d;

        double ival_mid = ival.mid();

        if( -sparsity_tol_ < ival_mid && ival_mid < sparsity_tol_ )
        {
          rest_poly.insert( (*curr).key(), ival, IntervalAdditionNoValueCheck::getObject() );
          polynomial_.erase( curr );

          if( previous_element_exist ) { curr = prev; ++curr; }
          else                           curr = polynomial_.begin();
        }
        else
        {
          (*curr).value() = ival_mid;
          rest_poly.insert( (*curr).key(), ival-ival_mid, IntervalAdditionNoValueCheck::getObject() );

          prev = curr;
          ++curr;
          previous_element_exist = true;
        }
      }

      /*
        Divide interval term.
      */

      rest_ /= d;

      /*
        Bound the restpolynomial (which includes round off and cut off errors!)
        and add the result to I.
      */

      Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

      if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = B_REST;

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in TM /= double: *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }

      rest_ += B_REST;

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM /= d --" << std::endl;
#endif
    }
    else std::cerr << "TM /= double: Division by zero." << std::endl;

    return *this;
  }

  TaylorModel_Impl& TaylorModel_Impl::divide_without_seizing_rounding_errors(const double& d)
  {
#ifdef HASHINFO
    std::cout << "Start -- TM /= d --" << std::endl;
#endif

    /*
      The value 'd' is equal to the Taylor Model d + [0,0].
    */

    if( d != 0.0 )
    {
      /*
        Divide polynomial with 'd'.
      */

      polynomial_iterator
        prev = polynomial_.begin(),
        curr = prev,
        last = polynomial_.end();

      bool previous_element_exist = false;
      while( curr != last )
      {
        double val  = (*curr).value();
        val /= d;

        if( -sparsity_tol_ < val && val < sparsity_tol_ )
        {
          polynomial_.erase( curr );

          if( previous_element_exist ) { curr = prev; ++curr; }
          else                           curr = polynomial_.begin();
        }
        else
        {
          (*curr).value() = val;

          prev = curr;
          ++curr;
          previous_element_exist = true;
        }
      }

      /*
        Adapt the size of the polynomial if necessary.
      */

      polynomial_.adapt_tab_size();

#ifdef HASHINFO
      polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- TM /= d --" << std::endl;
#endif
    }
    else std::cerr << "TM /= double: Division by zero." << std::endl;

    return *this;
  }

//
// The following function only subtracts the polynomial part, so the
// remainder term keeps unchanged.
//
  TaylorModel_Impl& TaylorModel_Impl::subtract_polynomial_part_of(const TaylorModel_Impl& s)
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = init_play_cntr(11);
    static unsigned int *record_cntr = init_record_cntr(11);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    /*
      First we have to check if the development points and the domains
      of the concerned variables are identical. This is also important
      for the record feature e.g. if the left hand side is a constant Taylor
      model that contains no data. With the check of the variable data we make
      sure that after this operation the data of the variables is contained in
      the result.
    */

    if( check_var_data( data_ , s.data_ ) )
    {
      /*
        Play mode. Use recorded information for calculating the remainder
        and leave this function.
      */

      if( Mode & REMAINDER_PLAY )
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in subtract_polynomial_part_of: *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

        rest_ += record[ (*play_cntr)++ ];

        return *this;
      }

      /*
        In the normal or record mode we calculate the polynomial and the remainder.
      */

      //
      // Polynomial with interval coefficients. Contains round off errors.
      //
      ipolynomial_type rest_poly( polynomial_.tab_size(), MonomHasher2::getFunction() );

      /*
        Subtract the polynomials:
        Iterate trough the polynomial of 's' and subtract every monom from the
        polynomial of 'this'.
      */

      //
      // What to do when inserting an existing monom?
      // The answer gives the object 'IntervalAdditionValueCheckWithContainer'.
      // Here we set the needed information.
      //
      IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly    );
      IntervalAdditionValueCheckWithContainer::getObject().setTol( sparsity_tol_ );

      const_polynomial_iterator
        curr = (s.polynomial_).begin(),
        last = (s.polynomial_).end();

      while( curr != last )
      {
        polynomial_.insert( (*curr).key(), -(*curr).value(), IntervalAdditionValueCheckWithContainer::getObject() );
        ++curr;
      }

      /*
        Bound the restpolynomial.
      */

      Interval B_REST = (*poly_eval_)( &rest_poly, &data_ );

      if( Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = B_REST;

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in subtract_polynomial_part_of: *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }

      rest_ += B_REST;
    }

    /*
      Adapt the size of the polynomial if necessary.
    */

    polynomial_.adapt_tab_size();

    return *this;
  }

/*

  Definition of friend functions

*/

  TaylorModel_Impl* sqr (const TaylorModel_Impl& s)
  {
#ifdef HASHINFO
    std::cout << "Start -- sqr(TM) --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(12);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(12);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    TaylorModel_Impl *result;

    if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in sqr(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

      result = new TaylorModel_Impl( 1, s.rest_, s.data_);

      result->rest_ += 2*record[ (*play_cntr)++ ]; // 2*B(Q)+I
      result->rest_ *= s.rest_;                    //(2*B(Q)+I)*I
      result->rest_ += record[ (*play_cntr)++ ];

      return result; //Leave function.
    }

    unsigned int sp_size = (s.polynomial_).tab_size();
    result = new TaylorModel_Impl( sp_size * sp_size, s.rest_, s.data_);

    /*
      TM_1 = P + I  => TM_1^2 = P*P + (2*B(P) + I)*I

      Iterate trough the polynomial of 's' and multiply every monom with the
      polynomial of 'this'.
    */

    //
    // Polynomial with interval coefficients.
    //
    TaylorModel_Impl::ipolynomial_type iprod_poly( sp_size, MonomHasher2::getFunction() );

    //
    // Polynomial with interval coefficients. Contains round off errors and cut off errors.
    //
    TaylorModel_Impl::ipolynomial_type rest_poly( sp_size * sp_size, MonomHasher2::getFunction() );

    TaylorModel_Impl::const_polynomial_iterator
      rhs_curr,
      rhs_last = (s.polynomial_).end(),
      lhs_curr = (s.polynomial_).begin(),
      lhs_last = rhs_last;

    /*
      x * Monom1 * y * Monom2 \in [x]*[y] * Monom3

      = mid([x]*[y]) * Monom3 + ([x]*[y] - mid([x]*[y])) * Monom3

      -----------+---------   ----------------+-----------------
      |                            |
      into Polynomial            into Restpolynomial
    */

    while( lhs_curr != lhs_last )
    {
      Monom  lhs_m = (*lhs_curr).key();
      double lhs_v = (*lhs_curr).value();

      rhs_curr = (s.polynomial_).begin();

      while( rhs_curr != rhs_last )
      {
        //Product of monomials (=Monom3, see above).
        Monom rhs_m = (*rhs_curr).key();
        rhs_m *= lhs_m;

        //Product of coefficients (\in [x]*[y]).
        Interval icoeff = (*rhs_curr).value();
        icoeff *= lhs_v;

        //Check degree of monomial rhs_m (=Monom3).
        if( TaylorModel_Impl::degree_->check( &rhs_m, TaylorModel_Impl::order_ ) )
        {
          iprod_poly.insert( rhs_m, icoeff, IntervalAdditionNoValueCheck::getObject() );
        }
        else
        {
          rest_poly .insert( rhs_m, icoeff, IntervalAdditionNoValueCheck::getObject() );
        }

        ++rhs_curr;
      }

      ++lhs_curr;
    }

    TaylorModel_Impl::const_ipolynomial_iterator
      curr = iprod_poly.begin(),
      last = iprod_poly.end();

    while( curr != last )
    {
      //Compute mid([x]*[y]).
      double mid_icoeff = (curr->value()).mid();

      if( mid_icoeff <= -TaylorModel_Impl::sparsity_tol_ || TaylorModel_Impl::sparsity_tol_ <= mid_icoeff )
      {
        result->polynomial_.insert( curr->key(), mid_icoeff                , DefaultAction<Monom,double>::getObject() );
        rest_poly.          insert( curr->key(), curr->value() - mid_icoeff, IntervalAdditionNoValueCheck::getObject() );
      }
      else
      {
        rest_poly.insert( curr->key(), curr->value(), IntervalAdditionNoValueCheck::getObject() );
      }

      ++curr;
    }

    /*
      Compute the interval term.
    */

    if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
    {
      //
      // Compute 2*B(Q) and add it to I.
      //
      Interval
        BSP = (*TaylorModel_Impl::poly_eval_)( &(s.polynomial_), &(s.data_) );

      result->rest_ += 2*BSP;

      //
      // Now multiply 2*B(Q)+I with I.
      //
      result->rest_ *= s.rest_;

      /*
        Bound the restpolynomial (which includes cut off and round off errors!)
        and add the result to I.
      */

      Interval B_REST = (*TaylorModel_Impl::poly_eval_)( &rest_poly, &s.data_ );

      if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = BSP;
        record[ (*record_cntr)++ ] = B_REST;

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in sqr(TM): *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }

      result->rest_ += B_REST;
    }

    /*
      Adapt the size of the polynomial if necessary.
    */

    result->polynomial_.adapt_tab_size();

#ifdef HASHINFO
    result->polynomial_.analyze_table_occupancy();
    std::cout << "Ende -- sqr(TM) --" << std::endl;
#endif

    return result;
  }

  TaylorModel_Impl* sqrt(const TaylorModel_Impl& s)
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(13);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(13);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    TaylorModel_Impl *result = 0;

    Interval B = (*TaylorModel_Impl::poly_eval_)( &(s.polynomial_), &(s.data_) ) + s.rest_;

    if( 0.0 < B.inf() || (TaylorModel_Impl::Mode & REMAINDER_PLAY) ) //Positive 's'?
    {
      /*
        Separate the polynomial of 's' into the constant and the non-constant part.
        At this point there exist a constant part, otherwise the Taylor Model 's'
        contains zero. But this was checked above.
      */

      double constant;
      TaylorModel_Impl *s_without_cnst = 0; //Taylor Model without constant part in polynomial.

      if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in sqrt(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        constant = record[ (*play_cntr)++ ].inf(); //The constant is saved as type interval. So we
        //have to use only one record array.

        result         = new TaylorModel_Impl( 1, Adapt::ZERO_INTERVAL(), s.data_);
        s_without_cnst = new TaylorModel_Impl( 1, s.rest_               , s.data_);
      }
      else //Other mode.
      {
        TaylorModel_Impl::const_polynomial_iterator
          pos = (s.polynomial_).find( Monom() ); //Find the position of the constant part ...
        constant = (*pos).value();                  //... and read the constant part from 's'.

        result         = new TaylorModel_Impl( (s.polynomial_).number_of_entries(), Adapt::ZERO_INTERVAL(), s.data_);
        s_without_cnst = new TaylorModel_Impl(s);       //Copy 's' ...
        (s_without_cnst->polynomial_).erase( Monom() ); //... and remove constant part.

        if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
        {
          record[ (*record_cntr)++ ] = constant; //We save the constant as type interval.

          if( *record_cntr > _RECORD_SIZE_ )
          {
            std::cout << "ERROR in sqrt(TM): *** Record overrun. *** " << std::endl;
            std::exit(1);
          }
        }
      }

      /*
        Now compute the polynomial part of 'sqrt(s)' with the horner scheme.
        The polynomial part is:

        sqrt(constant) * { 1 + s/(2*constant) - 1/2! (s/(2*constant))^2 + ... + (-1)^(n-1) (2n-3)!!/n! (s/(2*constant))^n }
      */

      (*s_without_cnst) /= constant; // = s/constant.
      (*s_without_cnst) /= 2.0;      // = s/(2*constant).

      if( TaylorModel_Impl::order_ > 1 )
      {
        double coeff;
        if( (TaylorModel_Impl::order_ & 1) != 0 )
          coeff =  1.0; // = (-1)^(n-1) for odd n.
        else
          coeff = -1.0; // = (-1)^(n-1) for even n.

        (*result) += coeff
          * Factorials::getTable().odd_factorial(2*TaylorModel_Impl::order_ - 3)
          / Factorials::getTable().factorial(TaylorModel_Impl::order_);
        for(unsigned int i = TaylorModel_Impl::order_-1; i > 1; i-- )
        {
          (*result) *= (*s_without_cnst);
          coeff     *= -1.0;
          (*result) += coeff
            * Factorials::getTable().odd_factorial(2*i - 3)
            / Factorials::getTable().factorial(i);
        }
      }
      else (*result) += 1.0;

      (*result) *= (*s_without_cnst);
      (*result) += 1.0;
      (*result) *= (*s_without_cnst);
      (*result) += 1.0;
      (*result) *= sqrt(Interval(constant));

      /*
        Then we compute the remainder term. It is:

        sqrt(constant) * (-1)^n * (2n-1)!!/(n+1)! * (s/(2*constant))^(n+1) * 1/( 1 + theta * s/constant )^(n + 1/2),  theta \in (0,1).
      */

      if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
      {
        Interval BS;
        if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
        {
          if( *record_cntr == 0 )
          {
            std::cout << "ERROR in sqrt(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
            std::exit(1);
          }

          BS = record[ (*play_cntr)++ ] + s_without_cnst->rest_;

          if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.
        }
        else //Other mode.
        {
          Interval
            bound_s_without_cnst = (*TaylorModel_Impl::poly_eval_)( &(s_without_cnst->polynomial_), &(s_without_cnst->data_) );

          BS = bound_s_without_cnst + s_without_cnst->rest_;

          if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
          {
            record[ (*record_cntr)++ ] = bound_s_without_cnst;

            if( *record_cntr > _RECORD_SIZE_ )
            {
              std::cout << "ERROR in sqrt(TM): *** Record overrun. *** " << std::endl;
              std::exit(1);
            }
          }
        }

        result->rest_ -= sqrt(Interval(constant))
          * Factorials::getTable().odd_factorial(2*TaylorModel_Impl::order_ - 1)
          * power(-BS,TaylorModel_Impl::order_+1)
          * exp( -(TaylorModel_Impl::order_ + 0.5) * ln( 1.0 + Interval(0,2) * BS ) )
          / Factorials::getTable().factorial(TaylorModel_Impl::order_ + 1);
      }

      //... that was it.

      delete s_without_cnst; //Delete copy of 's'.
    }
    else std::cerr << "sqrt(TM): Function not defined." << std::endl;

    return result;
  }

  TaylorModel_Impl* invsqrt(const TaylorModel_Impl& s)
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(14);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(14);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    TaylorModel_Impl *result = 0;

    Interval B = (*TaylorModel_Impl::poly_eval_)( &(s.polynomial_), &(s.data_) ) + s.rest_;

    if( 0.0 < B.inf() || (TaylorModel_Impl::Mode & REMAINDER_PLAY) ) //Positive 's'?
    {
      /*
        Separate the polynomial of 's' into the constant and the non-constant part.
        At this point there exist a constant part, otherwise the Taylor Model 's'
        contains zero. But this was checked above.
      */

      double constant;
      TaylorModel_Impl *s_without_cnst = 0; //Taylor Model without constant part in polynomial.

      if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in invsqrt(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        constant = record[ (*play_cntr)++ ].inf(); //The constant is saved as type interval. So we
        //have to use only one record array.

        result         = new TaylorModel_Impl( 1, Adapt::ZERO_INTERVAL(), s.data_);
        s_without_cnst = new TaylorModel_Impl( 1, s.rest_               , s.data_);
      }
      else //Other mode.
      {
        TaylorModel_Impl::const_polynomial_iterator
          pos = (s.polynomial_).find( Monom() ); //Find position of constant part ...
        constant = (*pos).value(); //... and read the constant.

        result         = new TaylorModel_Impl( (s.polynomial_).number_of_entries(), Adapt::ZERO_INTERVAL(), s.data_);
        s_without_cnst = new TaylorModel_Impl(s);       //Copy 's' ...
        (s_without_cnst->polynomial_).erase( Monom() ); //and remove constant part.

        if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
        {
          record[ (*record_cntr)++ ] = constant; //We save the constant as type interval.

          if( *record_cntr > _RECORD_SIZE_ )
          {
            std::cout << "ERROR in invsqrt(TM): *** Record overrun. *** " << std::endl;
            std::exit(1);
          }
        }
      }

      /*
        Now compute the polynomial part of '1/sqrt(s)' with the horner scheme.
        The polynomial part is:

        1/sqrt(constant) * { 1 - s/(2*constant) + 3!!/2! (s/(2*constant))^2 - ... + (-1)^n (2n-1)!!/n! (s/(2*constant))^n }
      */

      (*s_without_cnst) /= constant; // = s/constant.
      (*s_without_cnst) /= 2.0;      // = s/(2*constant).

      if( TaylorModel_Impl::order_ > 0 )
      {
        double coeff;
        if( (TaylorModel_Impl::order_ & 1) != 0 )
          coeff = -1.0; // = (-1)^n for odd n.
        else
          coeff =  1.0; // = (-1)^n for even n.

        (*result) += coeff
          * Factorials::getTable().odd_factorial(2*TaylorModel_Impl::order_ - 1)
          / Factorials::getTable().factorial(TaylorModel_Impl::order_);
        for(unsigned int i = TaylorModel_Impl::order_-1; i > 0; i-- )
        {
          (*result) *= (*s_without_cnst);
          coeff     *= -1.0;
          (*result) += coeff
            * Factorials::getTable().odd_factorial(2*i - 1)
            / Factorials::getTable().factorial(i);
        }
        (*result) *= (*s_without_cnst);
      }

      (*result) += 1.0;
      Interval onedivsqrtc = 1/sqrt(Interval(constant));
      (*result) *= onedivsqrtc;

      /*
        Then we compute the remainder term. It is:

        1/sqrt(constant) * (-1)^(n+1) * (2n+1)!!/(n+1)! * (s/(2*constant))^(n+1) * 1/( 1 + theta * s/constant )^(n + 3/2),
        theta \in (0,1).
      */

      if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
      {
        Interval BS, RP;
        if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
        {
          BS = record[ (*play_cntr)++ ] + s_without_cnst->rest_; //vorher wurde "s.rest_" addiert.;
          RP = record[ (*play_cntr)++ ];

          if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.
        }
        else //Other mode.
        {
          Interval
            bound_s_without_cnst = (*TaylorModel_Impl::poly_eval_)( &(s_without_cnst->polynomial_), &(s_without_cnst->data_) );

          BS = bound_s_without_cnst + s_without_cnst->rest_;
          RP = onedivsqrtc * Factorials::getTable().odd_factorial(2*TaylorModel_Impl::order_ + 1)
            / Factorials::getTable().    factorial(  TaylorModel_Impl::order_ + 1);

          if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
          {
            record[ (*record_cntr)++ ] = bound_s_without_cnst;
            record[ (*record_cntr)++ ] = RP;

            if( *record_cntr > _RECORD_SIZE_ )
            {
              std::cout << "ERROR in invsqrt(TM): *** Record overrun. *** " << std::endl;
              std::exit(1);
            }
          }
        }

        result->rest_ += RP * power(-BS,TaylorModel_Impl::order_+1)
          * exp( -(TaylorModel_Impl::order_ + Interval(3)/2.0) * ln( 1.0 + Interval(0,2) * BS ) );
      }

      //... that was it.

      delete s_without_cnst; //Delete copy of 's'.
    }
    else std::cerr << "invsqrt(TM): Function not defined." << std::endl;

    return result;
  }

  TaylorModel_Impl* exp (const TaylorModel_Impl& s)
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(15);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(15);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    TaylorModel_Impl *result = 0;

    /*
      Separate the polynomial of 's' into the constant and the non-constant part.
    */

    Interval constant;
    bool   s_copy = false;
    const TaylorModel_Impl *s_without_cnst = &s; //Taylor Model without constant part in polynomial.

    if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in exp(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      constant = record[ (*play_cntr)++ ]; //The constant is saved as type interval.

      result         = new TaylorModel_Impl( 1, Adapt::ZERO_INTERVAL(), s.data_);
      s_without_cnst = new TaylorModel_Impl( 1, s.rest_               , s.data_);
    }
    else //Other mode.
    {
      TaylorModel_Impl::const_polynomial_iterator
        pos = (s.polynomial_).find( Monom() ); //Find position of constant part ...

      constant = (*pos).value(); //... an read it. It may be zero!

      if( pos != (s.polynomial_).end() ) //There is a constant part in polynomial.
      {
        //
        // We must use a temporary pointer because 's_without_cnst' is const.
        //
        TaylorModel_Impl *tmp = new TaylorModel_Impl(s); //Copy 's'.
        (tmp->polynomial_).erase( Monom() ); //Remove constant part from the copy of 's'.

        s_copy         = true;
        s_without_cnst = tmp;
      }

      result   = new TaylorModel_Impl( (s.polynomial_).number_of_entries(), Adapt::ZERO_INTERVAL(), s.data_);

      if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = constant; //We save the constant as type interval.

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in exp(TM): *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }
    }

    /*
      Now compute the polynomial part of 'exp(s)' with the horner scheme.
      The polynomial part is:

      exp(constant) * { 1 + s + 1/2! * s^2 + ... + 1/n! * s^n }
    */

    (*result) += 1.0/Factorials::getTable().factorial(TaylorModel_Impl::order_);
    for(unsigned int i = TaylorModel_Impl::order_; i > 0; i-- )
    {
      (*result) *= (*s_without_cnst);
      (*result) += 1.0/Factorials::getTable().factorial(i-1);
    }
    (*result) *= exp(constant);

    /*
      Then we compute the remainder term. It is:

      exp(constant) * 1/(n+1)! * s^(n+1) * exp( theta * s ),  theta \in (0,1).
    */

    if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
    {
      Interval BS;
      if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in exp(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        BS = record[ (*play_cntr)++ ] + s_without_cnst->rest_;

        if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.
      }
      else //Other mode.
      {
        Interval
          bound_s_without_cnst = (*TaylorModel_Impl::poly_eval_)( &(s_without_cnst->polynomial_), &(s_without_cnst->data_) );

        BS = bound_s_without_cnst + s_without_cnst->rest_;

        if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
        {
          record[ (*record_cntr)++ ] = bound_s_without_cnst;

          if( *record_cntr > _RECORD_SIZE_ )
          {
            std::cout << "ERROR in exp(TM): *** Record overrun. *** " << std::endl;
            std::exit(1);
          }
        }
      }

      result->rest_ += exp(constant + Interval(0,1) * BS)
        / Factorials::getTable().factorial(TaylorModel_Impl::order_+1)
        * power(BS,TaylorModel_Impl::order_+1);
    }

    //... that was it.

    if( s_copy ) delete s_without_cnst; //Delete copy of 's'.

    return result;
  }

  TaylorModel_Impl* log (const TaylorModel_Impl& s)
  {
    std::cerr<<"log(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* sin (const TaylorModel_Impl& s)
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(16);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(16);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    TaylorModel_Impl *result = 0;

    /*
      Separate the polynomial of 's' into the constant and the non-constant part.
    */

    double constant;
    bool s_copy = false;
    const TaylorModel_Impl *s_without_cnst = &s; //Taylor Model without constant part in polynomial.

    if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in sin(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      constant = record[ (*play_cntr)++ ].inf(); //The constant is saved as type interval.

      result         = new TaylorModel_Impl( 1, Adapt::ZERO_INTERVAL(), s.data_);
      s_without_cnst = new TaylorModel_Impl( 1, s.rest_               , s.data_);
    }
    else //Other mode.
    {
      TaylorModel_Impl::const_polynomial_iterator
        pos = (s.polynomial_).find( Monom() ); //Find position of constant part ...

      constant = (*pos).value(); //... an read it. It may be zero!

      if( pos != (s.polynomial_).end() ) //There is a constant part in polynomial.
      {
        //
        // We must use a temporary pointer because 's_without_cnst' is const.
        //
        TaylorModel_Impl *tmp = new TaylorModel_Impl(s); //Copy 's'.
        (tmp->polynomial_).erase( Monom() ); //Remove constant part from the copy of 's'.

        s_copy         = true;
        s_without_cnst = tmp;
      }

      result = new TaylorModel_Impl( (s.polynomial_).number_of_entries(), Adapt::ZERO_INTERVAL(), s.data_);

      if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = constant; //We save the constant as type interval.

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in sin(TM): *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }
    }

    /*
      Now compute the polynomial part of 'sin(s)' with the horner scheme.
      The polynomial part is:

      sin(constant) + cos(constant) * s - 1/2! * sin(constant) * s^2 - 1/3! * cos(constant) * s^3 + ...

    */

    // TODO: korrekt? (Sk)
    Interval sin_c = ::sin(constant), cos_c = ::cos(constant);

    int k_mod_4 = (TaylorModel_Impl::order_ % 4);
    if     ( k_mod_4 == 0 ) (*result) += sin_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);
    else if( k_mod_4 == 1 ) (*result) += cos_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);
    else if( k_mod_4 == 2 ) (*result) -= sin_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);
    else if( k_mod_4 == 3 ) (*result) -= cos_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);

    for(unsigned int i = TaylorModel_Impl::order_; i > 0; i-- )
    {
      (*result) *= (*s_without_cnst);

      k_mod_4 = (i-1) % 4;
      if     ( k_mod_4 == 0 ) (*result) += sin_c/Factorials::getTable().factorial(i-1);
      else if( k_mod_4 == 1 ) (*result) += cos_c/Factorials::getTable().factorial(i-1);
      else if( k_mod_4 == 2 ) (*result) -= sin_c/Factorials::getTable().factorial(i-1);
      else if( k_mod_4 == 3 ) (*result) -= cos_c/Factorials::getTable().factorial(i-1);
    }

    /*
      Then we compute the remainder term. It is:

      1/(n+1)! * s^(n+1) * J,

      /                                      /
      |  -J_0, k mod 4 = 1,2                 |  cos(constant + theta * s), k is even
      with J = <                          and   J_0 = <                                            theta \in (0,1).
      |   J_0, else                          |  sin(constant + theta * s), else
      \                                      \
    */

    if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
    {
      Interval BS;
      if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in sin(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        BS = record[ (*play_cntr)++ ] + s_without_cnst->rest_;

        if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.
      }
      else //Other mode.
      {
        Interval
          bound_s_without_cnst = (*TaylorModel_Impl::poly_eval_)( &(s_without_cnst->polynomial_), &(s_without_cnst->data_) );

        BS = bound_s_without_cnst + s_without_cnst->rest_;

        if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
        {
          record[ (*record_cntr)++ ] = bound_s_without_cnst;

          if( *record_cntr > _RECORD_SIZE_ )
          {
            std::cout << "ERROR in sin(TM): *** Record overrun. *** " << std::endl;
            std::exit(1);
          }
        }
      }

      Interval remainder = power(BS,TaylorModel_Impl::order_+1)
        / Factorials::getTable().factorial(TaylorModel_Impl::order_+1);

      k_mod_4 = (TaylorModel_Impl::order_ % 4);
      if     ( k_mod_4 == 0 ) remainder *=  cos( constant + Interval(0,1) * BS );
      else if( k_mod_4 == 1 ) remainder *= -sin( constant + Interval(0,1) * BS );
      else if( k_mod_4 == 2 ) remainder *= -cos( constant + Interval(0,1) * BS );
      else if( k_mod_4 == 3 ) remainder *=  sin( constant + Interval(0,1) * BS );

      result->rest_ += remainder;
    }

    //... that was it.

    if( s_copy ) delete s_without_cnst; //Delete copy of 's'.

    return result;
  }

  TaylorModel_Impl* cos (const TaylorModel_Impl& s)
  {
    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(17);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(17);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    TaylorModel_Impl *result = 0;

    /*
      Separate the polynomial of 's' into the constant and the non-constant part.
    */

    double constant;
    bool s_copy = false;
    const TaylorModel_Impl *s_without_cnst = &s; //Taylor Model without constant part in polynomial.

    if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
    {
      if( *record_cntr == 0 )
      {
        std::cout << "ERROR in cos(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
        std::exit(1);
      }

      constant = record[ (*play_cntr)++ ].inf(); //The constant is saved as type interval.

      result         = new TaylorModel_Impl( 1, Adapt::ZERO_INTERVAL(), s.data_);
      s_without_cnst = new TaylorModel_Impl( 1, s.rest_               , s.data_);
    }
    else //Other mode.
    {
      TaylorModel_Impl::const_polynomial_iterator
        pos = (s.polynomial_).find( Monom() ); //Find position of constant part ...

      constant = (*pos).value(); //... an read it. It may be zero!

      if( pos != (s.polynomial_).end() ) //There is a constant part in polynomial.
      {
        //
        // We must use a temporary pointer because 's_without_cnst' is const.
        //
        TaylorModel_Impl *tmp = new TaylorModel_Impl(s); //Copy 's'.
        (tmp->polynomial_).erase( Monom() ); //Remove constant part from the copy of 's'.

        s_copy         = true;
        s_without_cnst = tmp;
      }

      result = new TaylorModel_Impl( (s.polynomial_).number_of_entries(), Adapt::ZERO_INTERVAL(), s.data_);

      if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
      {
        record[ (*record_cntr)++ ] = constant; //We save the constant as type interval.

        if( *record_cntr > _RECORD_SIZE_ )
        {
          std::cout << "ERROR in cos(TM): *** Record overrun. *** " << std::endl;
          std::exit(1);
        }
      }
    }

    /*
      Now compute the polynomial part of 'cos(s)' with the horner scheme.
      The polynomial part is:

      cos(constant) - sin(constant) * s - 1/2! * cos(constant) * s^2 + 1/3! * sin(constant) * s^3 + ...

    */

    // TODO: Korrekt? (SK)
    Interval sin_c = ::sin(constant), cos_c = ::cos(constant);

    int k_mod_4 = (TaylorModel_Impl::order_ % 4);
    if     ( k_mod_4 == 0 ) (*result) += cos_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);
    else if( k_mod_4 == 1 ) (*result) -= sin_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);
    else if( k_mod_4 == 2 ) (*result) -= cos_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);
    else if( k_mod_4 == 3 ) (*result) += sin_c/Factorials::getTable().factorial(TaylorModel_Impl::order_);

    for(unsigned int i = TaylorModel_Impl::order_; i > 0; i-- )
    {
      (*result) *= (*s_without_cnst);

      k_mod_4 = (i-1) % 4;
      if     ( k_mod_4 == 0 ) (*result) += cos_c/Factorials::getTable().factorial(i-1);
      else if( k_mod_4 == 1 ) (*result) -= sin_c/Factorials::getTable().factorial(i-1);
      else if( k_mod_4 == 2 ) (*result) -= cos_c/Factorials::getTable().factorial(i-1);
      else if( k_mod_4 == 3 ) (*result) += sin_c/Factorials::getTable().factorial(i-1);
    }

    /*
      Then we compute the remainder term. It is:

      1/(n+1)! * s^(n+1) * J,

      /                                      /
      |  -J_0, n mod 4 = 0,1                 |  sin(constant + theta * s), n is even
      with J = <                          and   J_0 = <                                            theta \in (0,1).
      |   J_0, else                          |  cos(constant + theta * s), else
      \                                      \
    */

    if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
    {
      Interval BS;
      if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
      {
        if( *record_cntr == 0 )
        {
          std::cout << "ERROR in cos(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
          std::exit(1);
        }

        BS = record[ (*play_cntr)++ ] + s_without_cnst->rest_;

        if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.
      }
      else //Other mode.
      {
        Interval
          bound_s_without_cnst = (*TaylorModel_Impl::poly_eval_)( &(s_without_cnst->polynomial_), &(s_without_cnst->data_) );

        BS = bound_s_without_cnst + s_without_cnst->rest_;

        if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
        {
          record[ (*record_cntr)++ ] = bound_s_without_cnst;

          if( *record_cntr > _RECORD_SIZE_ )
          {
            std::cout << "ERROR in cos(TM): *** Record overrun. *** " << std::endl;
            std::exit(1);
          }
        }
      }

      Interval remainder = power(BS,TaylorModel_Impl::order_+1)
        / Factorials::getTable().factorial(TaylorModel_Impl::order_+1);

      k_mod_4 = (TaylorModel_Impl::order_ % 4);
      if     ( k_mod_4 == 0 ) remainder *= -sin( constant + Interval(0,1) * BS );
      else if( k_mod_4 == 1 ) remainder *= -cos( constant + Interval(0,1) * BS );
      else if( k_mod_4 == 2 ) remainder *=  sin( constant + Interval(0,1) * BS );
      else if( k_mod_4 == 3 ) remainder *=  cos( constant + Interval(0,1) * BS );

      result->rest_ += remainder;
    }

    //... that was it.

    if( s_copy ) delete s_without_cnst; //Delete copy of 's'.

    return result;
  }

  TaylorModel_Impl* tan  (const TaylorModel_Impl& s)
  {
    std::cerr<<"tan(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* cot  (const TaylorModel_Impl& s)
  {
    std::cerr<<"cot(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* asin (const TaylorModel_Impl& s)
  {
    std::cerr<<"asin(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* acos (const TaylorModel_Impl& s)
  {
    std::cerr<<"acos(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* atan (const TaylorModel_Impl& s)
  {
    std::cerr<<"atan(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* acot (const TaylorModel_Impl& s)
  {
    std::cerr<<"acot(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* sinh (const TaylorModel_Impl& s)
  {
    std::cerr<<"sinh(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* cosh (const TaylorModel_Impl& s)
  {
    std::cerr<<"cosh(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* tanh (const TaylorModel_Impl& s)
  {
    std::cerr<<"tanh(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* coth (const TaylorModel_Impl& s)
  {
    std::cerr<<"coth(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* asinh(const TaylorModel_Impl& s)
  {
    std::cerr<<"asinh(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* acosh(const TaylorModel_Impl& s)
  {
    std::cerr<<"acosh(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* atanh(const TaylorModel_Impl& s)
  {
    std::cerr<<"atanh(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* acoth(const TaylorModel_Impl& s)
  {
    std::cerr<<"acoth(TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* pow  (const TaylorModel_Impl& s,const TaylorModel_Impl& t)
  {
    std::cerr<<"pow(TM,TM): Noch nicht implementiert"<<std::endl;
    return (TaylorModel_Impl*)(0);
  }

  TaylorModel_Impl* integrate (const TaylorModel_Impl& s,unsigned int var_code)
  {
#ifdef HASHINFO
    std::cout << "Start -- integrate(TM) --" << std::endl;
#endif

    /*
      Data for record feature.
    */

    static unsigned int *play_cntr   = TaylorModel_Impl::init_play_cntr(18);
    static unsigned int *record_cntr = TaylorModel_Impl::init_record_cntr(18);
    static Interval     *record      = new Interval[_RECORD_SIZE_];

    unsigned int p_size = 1;
    if( ! (TaylorModel_Impl::Mode & REMAINDER_PLAY) ) p_size = (s.polynomial_).tab_size();

    TaylorModel_Impl *result = new TaylorModel_Impl( p_size, Adapt::ZERO_INTERVAL(), s.data_);

    /*
      First we have to check if the Taylor Model contains information about the
      variable in which we should integrate.
    */

    if( 0 < var_code && var_code <= s.data_.size() && (s.data_[var_code-1]).operator ->() != 0 )
    {
      //
      // Contains round off and cut off errors.
      //
      TaylorModel_Impl::ipolynomial_type rest_poly( (s.polynomial_).tab_size(), MonomHasher2::getFunction() );

      /*
        Integrate the polynomial.
      */

      ValueCheckOnlyWithContainer::getObject().setPtr( &rest_poly                      );
      ValueCheckOnlyWithContainer::getObject().setTol( TaylorModel_Impl::sparsity_tol_ );

      TaylorModel_Impl::const_polynomial_iterator
        curr = (s.polynomial_).begin(),
        last = (s.polynomial_).end();

      while( curr != last )
      {
        //Integrate the monomial.
        Monom m = (*curr).key();
        unsigned int expnt = m.integrate_regarding_to(var_code);

        //Compute the coefficient.
        Interval icoeff = (*curr).value();
        icoeff /= expnt;

        //Check degree of integrated monomial.
        if( TaylorModel_Impl::degree_->check( &m, TaylorModel_Impl::order_ ) )
        {
          double mid_icoeff = icoeff.mid();

          (result->polynomial_).insert_if( m, mid_icoeff, ValueCheckOnlyWithContainer::getObject() );
          rest_poly.insert( m, icoeff - mid_icoeff, IntervalAdditionNoValueCheck::getObject() );
        }
        else
        {
          rest_poly.insert( m, icoeff, IntervalAdditionNoValueCheck::getObject() );//Cut off error.
        }

        ++curr;
      }

      if( ! (TaylorModel_Impl::Mode & NO_REMAINDER) )
      {
        /*
          Integrate the interval term.
        */

        result->rest_ += s.rest_ * ( (*(result->data_[var_code-1])).domain_ - (*(result->data_[var_code-1])).devel_point_ );

        /*
          Bound the restpolynomial (which includes round off and cut off errors!)
          and add the result to I.
        */

        Interval BR;
        if( TaylorModel_Impl::Mode & REMAINDER_PLAY ) //Play mode. Use recorded information for calculating the remainder.
        {
          if( *record_cntr == 0 )
          {
            std::cout << "ERROR in integrate(TM): *** Use PLAY Mode after RECORD Mode. ***" << std::endl;
            std::exit(1);
          }

          if( *play_cntr == *record_cntr ) *play_cntr = 0;//End of recorded data reached. Start from beginning.

          BR = record[ (*play_cntr)++ ];
        }
        else //Other mode.
        {
          BR = (*TaylorModel_Impl::poly_eval_)( &rest_poly, &(result->data_) );

          if( TaylorModel_Impl::Mode & REMAINDER_REC ) //REMAINDER_REC Mode. Save all data for calculating the remainder.
          {
            record[ (*record_cntr)++ ] = BR;

            if( *record_cntr > _RECORD_SIZE_ )
            {
              std::cout << "ERROR in integrate(TM): *** Record overrun. *** " << std::endl;
              std::exit(1);
            }
          }
        }

        result->rest_ += BR;
      }

      result->polynomial_.adapt_tab_size();

#ifdef HASHINFO
      result->polynomial_.analyze_table_occupancy();
      std::cout << "Ende -- integrate(TM) --" << std::endl;
#endif
    }
    else std::cerr << "INTEGRATE: The given TM contains no information about the respective variable" << std::endl;

    return result;
  }

  TaylorModel_Impl* derivate  (const TaylorModel_Impl& s,unsigned int var_code)
  {
    TaylorModel_Impl *result = new TaylorModel_Impl( (s.polynomial_).tab_size(),Adapt::ZERO_INTERVAL(),s.data_);

    /*
      First we have to check if the Taylor Model contains information about the
      variable in which we should derivate.
    */

    if( 0 < var_code && var_code <= s.data_.size() && (s.data_[var_code-1]).operator ->() != 0 )
    {
      //
      // Contains round off and cut off errors.
      //
      TaylorModel_Impl::ipolynomial_type rest_poly( (s.polynomial_).tab_size(), MonomHasher2::getFunction() );

      /*
        Derivate the polynomial.
      */

      ValueCheckOnlyWithContainer::getObject().setPtr( &rest_poly                      );
      ValueCheckOnlyWithContainer::getObject().setTol( TaylorModel_Impl::sparsity_tol_ );

      TaylorModel_Impl::const_polynomial_iterator
        curr = (s.polynomial_).begin(),
        last = (s.polynomial_).end();

      while( curr != last )
      {
        //Derivate the monomial.
        Monom m = (*curr).key();
        unsigned int expnt = m.derivate_regarding_to(var_code);

        if( expnt > 0 )
        {
          //Compute the coefficient.
          Interval icoeff = (*curr).value();
          icoeff *= expnt;

          //Check degree of derivated monomial.
          if( TaylorModel_Impl::degree_->check( &m, TaylorModel_Impl::order_ ) )
          {
            double mid_icoeff = icoeff.mid();

            (result->polynomial_).insert_if( m, mid_icoeff, ValueCheckOnlyWithContainer::getObject() );
            rest_poly.insert( m, icoeff - mid_icoeff, IntervalAdditionNoValueCheck::getObject() );//Round off error.
          }
          else
          {
            //Order of m is higer than 'order_'.
            rest_poly.insert( m, icoeff, IntervalAdditionNoValueCheck::getObject() );//Cut off error.
          }
        }

        ++curr;
      }

      /*
        Bound the restpolynomial (which includes round off and cut off errors!)
        and add the result to I.
      */

      result->rest_ += (*TaylorModel_Impl::poly_eval_)( &rest_poly, &(result->data_) );
    }

    return result;
  }

  TaylorModel_Impl* substitute(const TaylorModel_Impl& s,unsigned int var_code,const Interval& iv)
  {
#ifdef HASHINFO
    std::cout << "Start -- substitute(TM) --" << std::endl;
#endif

    /*
      To substitute a variable also means to delete it from the Taylor Model. So the
      information about the development point and the domain will be deleted. From now on
      this Taylor Model can operate with other Taylor Models that contain the subtituted
      variable with other information of them.
    */

    if( 0 < var_code && var_code <= s.data_.size() && (s.data_[var_code-1]).operator ->() != 0 )
    {
      if( iv.subset( (*(s.data_[var_code-1])).domain_ ) ) //Interval 'iv' is a subset of '...domain_'.
      {
        TaylorModel_Impl *result = new TaylorModel_Impl( (s.polynomial_).tab_size(),s.rest_,s.data_);

        Interval sub = iv - (*(s.data_[var_code-1])).devel_point_;

        //
        // Contains round off errors. No cut off errors occur during the substitution,
        // because the polynomial order does not grow.
        //
        TaylorModel_Impl::ipolynomial_type rest_poly( (s.polynomial_).tab_size(), MonomHasher2::getFunction() );

        /*
          Substitute variable with code 'var_code' in every monom.
        */

        IntervalAdditionValueCheckWithContainer::getObject().setPtr( &rest_poly                      );
        IntervalAdditionValueCheckWithContainer::getObject().setTol( TaylorModel_Impl::sparsity_tol_ );

        TaylorModel_Impl::const_polynomial_iterator
          curr = (s.polynomial_).begin(),
          last = (s.polynomial_).end();

        while( curr != last )
        {
          //The monomial.
          Monom m = (*curr).key();
          unsigned int expnt = m.substitute(var_code);

          //Compute the coefficient.
          Interval icoeff = (*curr).value();
          icoeff *= power( sub, expnt );

          double mid_icoeff = icoeff.mid();

          (result->polynomial_).insert_if( m, mid_icoeff, IntervalAdditionValueCheckWithContainer::getObject() );
          rest_poly.insert( m, icoeff - mid_icoeff, IntervalAdditionNoValueCheck::getObject() );//Round off error.

          ++curr;
        }

        /*
          Bound the restpolynomial (which includes round off and cut off errors!)
          and add the result to I.
        */

        result->rest_ += (*TaylorModel_Impl::poly_eval_)( &rest_poly, &(result->data_) );

        /*
          Delete the information of the substituded variable.
        */

        result->data_[var_code-1] = ReferenceCounter<VarDataListNode>();

        /*
          Adapt the size of the polynomial if necessary.
        */

        result->polynomial_.adapt_tab_size();

#ifdef HASHINFO
        result->polynomial_.analyze_table_occupancy();
        std::cout << "Ende -- substitute(TM) --" << std::endl;
#endif

        return result;
      }
      else std::cerr << "SUBSTITUTE: The given interval should be a subset of the respective domain" << std::endl;
    }

    return new TaylorModel_Impl(s);
  }

  TaylorModel_Impl* power(const TaylorModel_Impl& s,int n)
  {
    /*
      The cases n = 0,1,2 are discussed in 'power(TaylorModel,int)'
      (see file 'taylormodel.cpp').
    */

    TaylorModel_Impl *w, *z, *u = TaylorModel_Impl::Const_TM(1.0);

    if( TaylorModel_Impl::Mode & REMAINDER_PLAY )
    {
      z = new TaylorModel_Impl(1,s.rest_,s.data_);
    }
    else //Other modes.
    {
      z = new TaylorModel_Impl(s);
    }

    while( n > 0 )
    {
      if( (n & 1) == 0 )//even exponent
      {
        w = z;        //Point to the content of 'z'.
        z = sqr( *z );//Now 'z' points to a new content.
        delete w;     //Delete the old content of 'z'.
        n >>= 1;//Division by 2, realized with a right shift
      }
      else
      {
        (*u) *= (*z);
        n -= 1;
      }
    }
    delete z;

    return u;
  }

//
// And last but not least the output operators.
//
//
  std::ostream& operator <<  (std::ostream& os, const TaylorModel_Impl& s)
  {
    os << "\n" << "Polynomial = ";

    TaylorModel_Impl::const_polynomial_iterator
      curr = (s.polynomial_).begin(),
      last = (s.polynomial_).end();

    if( curr == last ) os << "0"; //Polynomial is zero.
    else
    {
      bool sign = false;
      while( curr != last )
      {
        double coeff = (*curr).value();
        if( sign )
        {
          if( coeff > 0 ) os << " + ";
          else
          {
            os << " - ";
            coeff = -coeff;
          }
        }
        else sign = true;

        //if( coeff != 1 || (*curr).key().is_const() ) os << coeff;

        os << coeff << " " << (*curr).key();

        ++curr;
      }
    }

    os << std::endl;

    os << "Interval = " << s.rest_ << std::endl;

    return os;
  }

  std::string&  operator <<  (std::string& st, const TaylorModel_Impl& s)
  {
    std::ostringstream out;
    out << s;
    return (st = out.str());
  }

/*

  End of File: taylormodel_impl.cpp

*/
}

