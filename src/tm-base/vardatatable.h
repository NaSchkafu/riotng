/*

 File: vardatatable.h, 2004/08/18

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

#ifndef VARDATATABLE_H_INCLUDED
#define VARDATATABLE_H_INCLUDED

#include <vector>
#include "adaptintval.h"
#include "vartable.h"

using Adapt::Interval;

/*
  This Tabular was developed for use with reference counters. That means, when a 
  List Node will be deletet by a reference counter object it will chain
  the previous and next ListNode. Therefor a Base_ListNode is needed.
*/

class VarDataListNode
{
  friend class VarData;

  public:

  double   devel_point_;
  Interval domain_;

  unsigned int *cnt; //For reference counting. Will be destroyed by the ReferenceCounter object.

  VarDataListNode(double dp, Interval dm, VarDataListNode *p, VarDataListNode *n) 
    : devel_point_(dp), domain_(dm), cnt(new unsigned int(0)), prev_(p), next_(n)
    {}
  ~VarDataListNode() //Chains the previous and the next VarDataListNode before deleting itself.
    {
      if( prev_ ) prev_->next_ = next_;
      if( next_ ) next_->prev_ = prev_;
    }

 private:

  /*
    Data.
  */

  VarDataListNode *prev_, *next_;
};

/*
  The VarData is implementet as Singleton.
*/

class VarData
{
 private:

  /*
    To make it a Singleton.
  */

  VarData() : Tab_(100,(VarDataListNode*)(0)) {}
  VarData(const VarData&);
  ~VarData() //Deletes the base VarDataListNodes. The others will be deletet
    {        //by the reference counter objects.
      for(unsigned int i = 0; i < Tab_.size(); i++) 
	{
	  if( Tab_[i] ) 
	    {
	      delete Tab_[i]->cnt; //Delete unused counter of base node.
	      delete Tab_[i];      //Delete base node.
	    }
	}
    }
  VarData& operator = (const VarData&);

 public:

  static VarData& getTable()
    {
      static VarData singleton;
      return singleton;
    }

  VarDataListNode* insert(unsigned int var_code, const double& dp, const Interval& dm)
    {
      if( var_code > Tab_.size() ) Tab_.resize( 2*Tab_.size(), (VarDataListNode*)(0) );

      if( Tab_[ var_code - 1] == 0 ) //Insert a base node.
	Tab_[ var_code - 1 ] = new VarDataListNode(0.0,Interval(0.0),(VarDataListNode*)(0),(VarDataListNode*)(0));

      VarDataListNode *cur = Tab_[ var_code - 1 ]->next_;
      while( cur )
	{
	  if( cur->devel_point_ == dp && cur->domain_ == dm ) return cur;
	  cur = cur->next_;
	}
      
      //Insert a new VarDataListNode at the front behind the base VarDataListNode.
      cur = new VarDataListNode(dp, dm, Tab_[ var_code - 1 ], Tab_[ var_code - 1 ]->next_ );
      Tab_[ var_code - 1 ]->next_ = cur;
      if( cur->next_ ) (cur->next_)->prev_ = cur;
      return cur;
    }

 private:

  /*
    Data.
  */

  std::vector<VarDataListNode*> Tab_;
};

#endif

/*

  End of File: vardatatable.h

*/
