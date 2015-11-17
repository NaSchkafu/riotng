/*

  File: hashtable.h, 2003/06/24

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

#ifndef HASHTABLE_H_INCLUDED
#define HASHTABLE_H_INCLUDED

#include <vector>
#include <functional>
#include <utility>

#include "primetable.h"

/*-----------------------------------------------------------------------------+
  |                  Hash Table and Hash Table Iterator                         |
  +-----------------------------------------------------------------------------*/


/*
  Data to store in a hash table. Consists of a key and value.
*/
namespace riot
{

  template<typename K, typename T>
    class HashTableData
  {
    K key_;
    T value_;

  public:

    // Constructors

  HashTableData()                       : value_(0) {}
  HashTableData(const K& k, const T& v) : key_(k), value_(v) {}
  HashTableData(const HashTableData& x) : key_(x.key_), value_(x.value_) {}

    ~HashTableData() {}

    HashTableData& operator = (const HashTableData& x)
      {
        key_   = x.key_;
        value_ = x.value_;

        return *this;
      }

    // Access

    const K&   key() const { return   key_; }
    const T& value() const { return value_; }
    T&       value()       { return value_; }
  };

/*
  An Element of the hash table consisting of a pointer to the following
  element and an data-object.
*/

  template<typename K, typename T>
    struct HashTableEntry
    {
      HashTableData<K,T>  data;
      HashTableEntry     *next;

      // Constructors

    HashTableEntry() : next((HashTableEntry*)(0)) {}
    HashTableEntry(const K& key, const T& value, HashTableEntry* n)
    : data(key, value), next(n) {}
    HashTableEntry(const HashTableEntry& x)
    : data(x.data), next((HashTableEntry*)(0)) {}

      ~HashTableEntry() {}

      HashTableEntry& operator = (const HashTableEntry& x)
        {
          data = x.data;
          next = (HashTableEntry*)(0);

          return *this;
        }
    };

/*
  Next we define iterator objects to walk trough a
  hash table.
*/

  template<typename K, typename T>
    class HashTableConstIterator;

  template<typename K, typename T>
    class HashTableIterator
  {
    friend class HashTableConstIterator<K,T>;

  public:

    // Introduce better readable names

    typedef HashTableIterator<K,T>             iterator;

    typedef unsigned int                       size_type;
    typedef HashTableData<K,T>                 value_type;
    typedef HashTableData<K,T>&                reference;
    typedef HashTableData<K,T>*                pointer;
    typedef HashTableEntry<K,T>                Node;

    typedef std::vector<HashTableEntry<K,T>* > Tabular;

    // Constructors

  HashTableIterator()
    : _M_node_((Node*)(0)), _M_htab_((Tabular*)(0)), tab_position_(0), tab_size_(0) {}
  HashTableIterator(Node *p, Tabular *t, size_type tp)
    : _M_node_(p), _M_htab_(t), tab_position_(tp)
    {
      tab_size_ = (*_M_htab_).size();
    }
  HashTableIterator(const iterator& x)
    : _M_node_(x._M_node_), _M_htab_(x._M_htab_), tab_position_(x.tab_position_), tab_size_(x.tab_size_) {}

    // Operators for comparison

    bool operator == (const iterator& x) const { return _M_node_ == x._M_node_; }
    bool operator != (const iterator& x) const { return _M_node_ != x._M_node_; }

    bool is_equal_to(Node* x) const { return _M_node_ == x; }

    // Access

    reference operator * () const { return   (*_M_node_).data ; }
    pointer   operator ->() const { return &((*_M_node_).data); }

    // Increment operators

    iterator& operator ++ ();
    iterator  operator ++ (int)
    {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    // Access size of table

    size_type ht_pos() const { return tab_position_; }

  private:

    Node            *_M_node_;
    Tabular         *_M_htab_;
    size_type        tab_position_;
    size_type        tab_size_;
  };

  template<typename K, typename T>
    typename HashTableIterator<K,T>::iterator& HashTableIterator<K,T>::operator ++ ()
  {
    if( (_M_node_ = _M_node_->next) ) return *this; //_M_node_->next is not NULL

    while( ! (_M_node_ = (*_M_htab_)[ (++tab_position_) %= tab_size_ ]) )
    {}

    return *this;
  }

  template<typename K, typename T>
    class HashTableConstIterator
  {
  public:

    // Introduce better readable names

    typedef HashTableIterator<K,T>             iterator;
    typedef HashTableConstIterator<K,T>        const_iterator;

    typedef unsigned int                       size_type;
    typedef HashTableData<K,T>                 value_type;
    typedef const HashTableData<K,T>&          reference;
    typedef const HashTableData<K,T>*          pointer;
    typedef HashTableEntry<K,T>                Node;

    typedef std::vector<HashTableEntry<K,T>* > Tabular;

    // Constructors

  HashTableConstIterator()
    : _M_node_((Node*)(0)), _M_htab_((Tabular*)(0)), tab_position_(0), tab_size_(0) {}
  HashTableConstIterator(const Node *p, const Tabular *t, size_type tp)
    : _M_node_(p), _M_htab_(t), tab_position_(tp)
    {
      tab_size_ = (*_M_htab_).size();
    }
  HashTableConstIterator(const iterator& x)
    : _M_node_(x._M_node_), _M_htab_(x._M_htab_), tab_position_(x.tab_position_), tab_size_(x.tab_size_) {}

    // Operators for comparison

    bool operator == (const const_iterator& x) const { return _M_node_ == x._M_node_; }
    bool operator != (const const_iterator& x) const { return _M_node_ != x._M_node_; }

    bool is_equal_to(Node* x) const { return _M_node_ == x; }

    // Access

    reference operator * () const { return   (*_M_node_).data ; }
    pointer   operator ->() const { return &((*_M_node_).data); }

    // Increment operators

    const_iterator& operator ++ ();
    const_iterator  operator ++ (int)
    {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    // Access size of table

    size_type ht_pos() const { return tab_position_; }

  private:

    const Node      *_M_node_;
    const Tabular   *_M_htab_;
    size_type        tab_position_;
    size_type        tab_size_;
  };

  template<typename K, typename T>
    typename HashTableConstIterator<K,T>::const_iterator& HashTableConstIterator<K,T>::operator ++ ()
  {
    if( (_M_node_ = _M_node_->next) ) return *this; //_M_node_->next is not NULL

    while( ! (_M_node_ = (*_M_htab_)[ (++tab_position_) %= tab_size_ ]) )
    {}

    return *this;
  }

/*-----------------------------------------------------------------------+
  | Action when inserting elements with same key.                         |
  | Implemented with Barton and Nackman Trick.                            |
  +-----------------------------------------------------------------------*/

  template<typename K, typename T, typename Leaf_Type>
    class BasicAction
  {
  public:

    //Delegate to leaf.
    bool operator () (HashTableData<K,T>& x, const HashTableData<K,T>& y) { return (static_cast<Leaf_Type&>(*this))(x,y); }
    bool operator () (const K& k, const T& v)                             { return (static_cast<Leaf_Type&>(*this))(k,v); }
  };

  template<typename K, typename T>
    class DefaultAction : public BasicAction<K,T,DefaultAction<K,T> >
  {
  private:

    /*
      To make it a Singleton.
    */

    DefaultAction() {}
    ~DefaultAction() {}
    DefaultAction(const DefaultAction&);
    DefaultAction& operator = (const DefaultAction&);

  public:

    static DefaultAction& getObject()
    {
      static DefaultAction singleton;
      return singleton;
    }

    bool operator () (HashTableData<K,T>& x, const HashTableData<K,T>& y) { return true; }
    bool operator () (const K& k, const T& v)                             { return true; }
  };


/*--------------------------------------------+
  |                                            |
  |          Hashtable-Implementation          |
  |                                            |
  +--------------------------------------------*/


  template<typename K, typename T, typename H, typename EQ=std::equal_to<K> >
    class HashTable
    {
    public:

    // Introduce better readable names

    typedef HashTableIterator<K,T>        iterator;
    typedef HashTableConstIterator<K,T>   const_iterator;

    typedef unsigned int                  size_type;
    typedef K                             key_type;
    typedef T                             value_type;
    typedef H                             hasher;
    typedef EQ                            key_equal;
    typedef HashTableEntry<K,T>           Node;

    /*
      Constructors
    */

    HashTable(size_type s, H& h, const EQ& e=EQ())
    : entries_(0), tab_size_(next_prime(s)), hash_(&h), eq_(e), ht_(tab_size_,(Node*)(0))
    {
      set_load();
      ht_[ tab_size_-1 ] = (_M_node_ = new Node);
    }
    HashTable(const HashTable& x)
    : max_load_(x.max_load_), grow_(x.grow_), entries_(x.entries_), tab_size_(x.tab_size_), hash_(x.hash_), ht_(tab_size_,(Node*)(0))
    {
      copy_table(x.ht_);
    }

    ~HashTable() { clear(); }

    HashTable& operator = (const HashTable& x);

    /*
      Iteration through the table
    */

    iterator       begin();
    const_iterator begin() const;
    iterator         end()       { return       iterator(_M_node_,&ht_,tab_size_-1); }
    const_iterator   end() const { return const_iterator(_M_node_,&ht_,tab_size_-1); }

    /*
      Elementary funtions
    */

    iterator       find(const key_type&);
    const_iterator find(const key_type&) const;

    template<typename OP> std::pair<iterator,bool>
    insert(const key_type& key, const value_type& value,
           BasicAction<key_type,value_type,OP>& op)
    {
      resize(entries_ + 1);
      return insert_noresize(key,value,op);
    }
    template<typename OP> std::pair<iterator,bool>
    insert_if(const key_type& key, const value_type& value,
              BasicAction<key_type,value_type,OP>& op)
    {
      resize(entries_ + 1);
      return insert_if_noresize(key,value,op);
    }

    size_type erase(const key_type&);
    void erase(const iterator&);
    void erase(const const_iterator&);

    /*
      State
    */

    size_type number_of_entries() const    { return entries_;          }
    size_type tab_size()          const    { return tab_size_;         }
    void set_load(float m=0.7,float g=1.6) { max_load_ = m; grow_ = g; }
    void analyze_table_occupancy() const;
    void adapt_tab_size();
    void resize(size_type);

    private:

    template<typename OP> std::pair<iterator,bool>
    insert_noresize   (const key_type&, const value_type&, BasicAction<key_type,value_type,OP>& op);
    template<typename OP> std::pair<iterator,bool>
    insert_if_noresize(const key_type&, const value_type&, BasicAction<key_type,value_type,OP>& op);

    void clear();
    void copy_table(const std::vector<Node*>&);
    size_type next_prime(const double& ul)
    {
      return primtab_->find_SmallestPrimeGreaterThan(ul);
    }

    /*
      Data.
    */

    float max_load_; //condition: entries <= max_load_ * tab_size_
    float grow_;     //for resizing

    size_type entries_;
    size_type tab_size_;

    hasher    *hash_; //Pointer to hashing object.
    key_equal  eq_;

    Node *_M_node_;

    std::vector<Node*> ht_; //Hash Table.

    static PrimeNumbers *primtab_;
    };

  template<typename K, typename T, typename H, typename EQ>
    PrimeNumbers *HashTable<K,T,H,EQ>::primtab_ = &PrimeNumbers::getTable();

  template<typename K, typename T, typename H, typename EQ>
    HashTable<K,T,H,EQ>& HashTable<K,T,H,EQ>::operator = (const HashTable<K,T,H,EQ>& x)
    {
      if( &x == this ) return *this;

      hash_ = x.hash_;

      clear();
      ht_.resize(x.tab_size_,(Node*)(0));
      copy_table(x.ht_);

      max_load_ = x.max_load_;
      grow_     = x.grow_;
      entries_  = x.entries_;
      tab_size_ = x.tab_size_;

      return *this;
    }

  template<typename K, typename T, typename H, typename EQ>
    typename HashTable<K,T,H,EQ>::iterator HashTable<K,T,H,EQ>::begin()
  {
    if( entries_ )
    {
      Node *tmp;
      int tab_position = -1;

      while( ! (tmp = ht_[++tab_position]) ) {} //Search for an entry.
      return iterator(tmp,&ht_,tab_position);
    }
    else return end();
  }

  template<typename K, typename T, typename H, typename EQ>
    typename HashTable<K,T,H,EQ>::const_iterator HashTable<K,T,H,EQ>::begin() const
  {
    if( entries_ )
    {
      Node *tmp;
      int tab_position = -1;

      while( ! (tmp = ht_[++tab_position]) ) {} //Search for an entry.
      return const_iterator(tmp,&ht_,tab_position);
    }
    else return end();
  }

/*
  Die 'find' Funktion sucht einen Eintrag mit dem Schluessel 'key' und
  liefert einen Iterator darauf zurück. Ist der gesuchte Eintrag nicht vorhanden
  zeigt der Iterator hinter das letzte Element der Hash-Tabelle.
*/

  template<typename K, typename T, typename H, typename EQ>
    typename HashTable<K,T,H,EQ>::const_iterator HashTable<K,T,H,EQ>::find(const key_type& key) const
  {
    const size_type tab_position = (*hash_)(key,tab_size_);
    Node *tmp = ht_[ tab_position ];
    while( tmp )
    {
      if( eq_(key,tmp->data.key()) ) return const_iterator(tmp,&ht_,tab_position);
      tmp = tmp->next;
    }

    return const_iterator(_M_node_,&ht_,tab_size_-1); //Nothing found.
  }

  template<typename K, typename T, typename H, typename EQ>
    typename HashTable<K,T,H,EQ>::iterator HashTable<K,T,H,EQ>::find(const key_type& key)
  {
    const size_type tab_position = (*hash_)(key,tab_size_);
    Node *tmp = ht_[ tab_position ];
    while( tmp )
    {
      if( eq_(key,tmp->data.key()) ) return iterator(tmp,&ht_,tab_position);
      tmp = tmp->next;
    }

    return iterator(_M_node_,&ht_,tab_size_-1); //nothing found
  }

/*
  Die 'erase(key)' Funktion loescht alle Eintraege mit dem Schluessel 'key',
  und liefert die Anzahl der geloeschten Eintraege zurück.
*/

  template<typename K, typename T, typename H, typename EQ>
    typename HashTable<K,T,H,EQ>::size_type HashTable<K,T,H,EQ>::erase(const key_type& key)
  {
    const size_type n = (*hash_)(key,tab_size_);
    Node *first = ht_[n];
    size_type erased = 0;

    if( first )
    {
      Node *current = first;
      Node *next    = current->next;
      while( next )
      {
        if( eq_( next->data.key(), key ) )
        {
          current->next = next->next;
          delete next;
          next = current->next;
          ++erased;
          --entries_;
        }
        else
        {
          current = next;
          next    = current->next;
        }
      }
      if( eq_( first->data.key(), key ) )
      {
        ht_[n] = first->next;
        delete first;
        ++erased;
        --entries_;
      }
    }
    return erased;
  }

/*
  Die 'erase(iterator)' Funktion loescht den Eintrag, welcher vom Iterator
  referenziert wird, aus der Hash Tabelle. Der Iterator wird dadurch ungültig!
*/

  template<typename K, typename T, typename H, typename EQ>
    void HashTable<K,T,H,EQ>::erase(const typename HashTable<K,T,H,EQ>::iterator& it)
  {
    if( it != end() )
    {
      const size_type tab_position = it.ht_pos();
      Node *current = ht_[ tab_position ];

      if( it.is_equal_to(current) )
      {
        ht_[ tab_position ] = current->next;
        delete current;
        --entries_;
      }
      else
      {
        Node *next = current->next;
        while( next )
        {
          if( it.is_equal_to(next) )
          {
            current->next = next->next;
            delete next;
            --entries_;
            break;
          }
          else
          {
            current = next;
            next    = current->next;
          }
        }
      }
    }
  }

  template<typename K, typename T, typename H, typename EQ>
    void HashTable<K,T,H,EQ>::erase(const typename HashTable<K,T,H,EQ>::const_iterator& it)
  {
    if( it != end() )
    {
      const size_type tab_position = it.ht_pos();
      Node *current = ht_[ tab_position ];

      if( it.is_equal_to(current) )
      {
        ht_[ tab_position ] = current->next;
        delete current;
        --entries_;
      }
      else
      {
        Node *next = current->next;
        while( next )
        {
          if( it.is_equal_to(next) )
          {
            current->next = next->next;
            delete next;
            --entries_;
            break;
          }
          else
          {
            current = next;
            next    = current->next;
          }
        }
      }
    }
  }

/*
  Die Funktion 'insert_noresize(key,value,op)' sucht nach einem Eintrag mit dem
  Schluessel 'key' in der Hash Tabelle. Ist 'key' nicht vorhanden
  wird das Element eingefuegt. Ist 'key' vorhanden wird das Funktionsobjekt 'op'
  auf den Knoten angewendet (Aufruf von 'op(node->data,data)') und liefert 'true'
  wenn alles in Ordnung ist und 'false' wenn der Knoten mit
  Schluessel 'key' geloescht werden soll.
*/

  template<typename K, typename T, typename H, typename EQ>
    template<typename OP>
    std::pair<typename HashTable<K,T,H,EQ>::iterator,bool>
    HashTable<K,T,H,EQ>::insert_noresize(const typename HashTable<K,T,H,EQ>::key_type& key,
                                         const typename HashTable<K,T,H,EQ>::value_type& value,
                                         BasicAction<K,T,OP>& op                                    )
  {
    const size_type tab_position = (*hash_)(key,tab_size_);

    Node *first = ht_[ tab_position ];

    if( first ) //The pointer 'first' is not null.
    {
      if( eq_( first->data.key(), key ) )
      {
        if( op( first->data, HashTableData<K,T>(key,value) ) )
          return std::pair<iterator, bool>(iterator(first,&ht_,tab_position), false);
        else
        {
          //remove
          ht_[ tab_position ] = first->next;
          delete first;
          --entries_;
          return std::pair<iterator, bool>(end(), false);
        }
      }

      for(Node *prev = first, *current = first->next; current; prev = current, current = current->next )
      {
        if( eq_( current->data.key(), key ) )
        {
          if( op( current->data, HashTableData<K,T>(key,value) ) )
            return std::pair<iterator, bool>(iterator(current,&ht_,tab_position), false);
          else
          {
            //remove
            prev->next = current->next;
            delete current;
            --entries_;
            return std::pair<iterator, bool>(end(), false);
          }
        }
      }
      //not found
    }

    Node *tmp = new Node(key,value,first);
    ht_[ tab_position ] = tmp;
    ++entries_;
    return std::pair<iterator, bool>(iterator(tmp,&ht_,tab_position), true);
  }

/*
  Die Funktion 'insert_if_noresize(key,value,op)' sucht nach einem Eintrag mit dem
  Schluessel 'key' in der Hash Tabelle. Ist 'key' nicht vorhanden
  wird das Funktionsobjekt 'op' befragt, ob das Element eingefuegt werden soll oder nicht.
  Ist 'key' vorhanden wird das Funktionsobjekt 'op'
  auf den Knoten angewendet (Aufruf von 'op(node->data,data)') und liefert 'true'
  wenn alles in Ordnung ist und 'false' wenn der Knoten mit
  Schluessel 'key' geloescht werden soll.
*/

  template<typename K, typename T, typename H, typename EQ>
    template<typename OP>
    std::pair<typename HashTable<K,T,H,EQ>::iterator,bool>
    HashTable<K,T,H,EQ>::insert_if_noresize(const typename HashTable<K,T,H,EQ>::key_type& key,
                                            const typename HashTable<K,T,H,EQ>::value_type& value,
                                            BasicAction<K,T,OP>& op                                     )
  {
    const size_type tab_position = (*hash_)(key,tab_size_);

    Node *first = ht_[ tab_position ];

    if( first ) //The pointer 'first' is not null.
    {
      if( eq_( first->data.key(), key ) )
      {
        if( op( first->data, HashTableData<K,T>(key,value) ) )
          return std::pair<iterator, bool>(iterator(first,&ht_,tab_position), false);
        else
        {
          //remove
          ht_[ tab_position ] = first->next;
          delete first;
          --entries_;
          return std::pair<iterator, bool>(end(), false);
        }
      }

      for(Node *prev = first, *current = first->next; current; prev = current, current = current->next )
      {
        if( eq_( current->data.key(), key ) )
        {
          if( op( current->data, HashTableData<K,T>(key,value) ) )
            return std::pair<iterator, bool>(iterator(current,&ht_,tab_position), false);
          else
          {
            //remove
            prev->next = current->next;
            delete current;
            --entries_;
            return std::pair<iterator, bool>(end(), false);
          }
        }
      }
      //not found
    }

    if( op( key, value ) )
    {
      Node *tmp = new Node(key,value,first);
      ht_[ tab_position ] = tmp;
      ++entries_;
      return std::pair<iterator, bool>(iterator(tmp,&ht_,tab_position), true);
    }
    else return std::pair<iterator, bool>(end(), false);
  }

/*
  Die 'resize' Funktion vergroessert bei Bedarf die Hash Tabelle.
*/

  template<typename K, typename T, typename H, typename EQ>
    void HashTable<K,T,H,EQ>::resize(size_type num_entries_hint)
  {
    if( num_entries_hint > tab_size_ * max_load_ )
    {
      const size_type n = next_prime( tab_size_ * grow_ );
      if( n > tab_size_ )
      {
#ifdef HASHINFO
        std::cout << "Hashtable resized" << std::endl;
#endif

        std::vector<Node*> tmp(n,(Node*)(0));

        try
        {
          tmp[ n - 1 ] = _M_node_; //Copy the _M_node_ separatly.

          for(size_type i = 0; i < tab_size_; i++)
          {
            Node *curr = ht_[i];
            while( curr && curr != _M_node_ )
            {
              size_type new_tab_position = (*hash_)( curr->data.key(), n );
              ht_[i] = curr->next;
              curr->next = tmp[ new_tab_position ];
              tmp[ new_tab_position ] = curr;
              curr = ht_[i];
            }
          }

          tmp.swap(ht_);//Swap the handles.
          tab_size_ = n;
        }
        catch(...)
        {
          for(size_type i = 0; i < tmp.size(); i++)
          {
            while( tmp[i] )
            {
              Node *next = tmp[i]->next;
              delete tmp[i];
              tmp[i] = next;
            }
          }
          throw;
        }
      }
    }
  }

/*
  Die 'clear' Funktion löscht alle Eintraege der Tabelle und bringt
  diese auf einen definierten Anfangszustand.
*/

  template<typename K, typename T, typename H, typename EQ>
    void HashTable<K,T,H,EQ>::clear()
  {
    for(unsigned int i=0; i < tab_size_; i++)
    {
      Node *current = ht_[i];
      while( current )
      {
        Node *next = current->next;
        delete current;
        current = next;
      }
      ht_[i] = (Node*)(0);
    }
  }

/*
  Die 'copy_table' Funktion kopiert alle Listen der uebergebenen Tabelle.
  Dies setzt voraus, dass die Zieltabelle diesselbe Größe besitzt.
*/

  template<typename K, typename T, typename H, typename EQ>
    void HashTable<K,T,H,EQ>::copy_table(const std::vector<Node*>& xht)
  {
    //When we copy the Hash Table we have to get the adress of the '_M_node_' in the new Hash Table.
    //So we copy the last list separatly (this saves 'x.tab_size_-1' if-conditions).

    Node *x_current, *current;
    unsigned int size = xht.size();

    for(unsigned int i = 0; i < size-1; i++)
    {
      x_current = xht[i];

      if( x_current )
      {
        Node *current = (ht_[i] = new Node(*x_current));
        x_current = x_current->next;

        while( x_current )
        {
          current->next = new Node(*x_current);
          current   =   current->next;
          x_current = x_current->next;
        }
      }
    }

    //now copy the last list separatly

    x_current =   xht[ size-1 ];
    current   = ( ht_[ size-1 ] = new Node(*x_current) );
    x_current = x_current->next;

    while( x_current )
    {
      current->next = new Node(*x_current);
      current   =   current->next;
      x_current = x_current->next;
    }
    _M_node_ = current;
  }

/*
  Die Funktion 'adapt_tab_size' passt die Größe der Hashtabelle auf die
  Anzahl der tatsächlichen Einträge an.
*/

  template<typename K, typename T, typename H, typename EQ>
    void HashTable<K,T,H,EQ>::adapt_tab_size()
  {
    const size_type n = next_prime( entries_ * grow_ );

    if( n < tab_size_ )
    {
#ifdef HASHINFO
      std::cout << "Hashtable adapted" << std::endl;
#endif

      std::vector<Node*> tmp(n,(Node*)(0));

      tmp[ n - 1 ] = _M_node_; //Copy the _M_node_ separatly.

      for(size_type i = 0; i < tab_size_; i++)
      {
        Node *curr = ht_[i];
        while( curr && curr != _M_node_ )
        {
          size_type new_tab_position = (*hash_)( curr->data.key(), n );
          ht_[i] = curr->next;
          curr->next = tmp[ new_tab_position ];
          tmp[ new_tab_position ] = curr;
          curr = ht_[i];
        }
      }

      tmp.swap(ht_);//Swap the handles.
      tab_size_ = n;
    }
  }

/*
  Die folgende Funktion stellt Daten zur Hash Tabelle zusammen
  und gibt diese in uebersichtlicher Form aus.
*/

  template<typename K, typename T, typename H, typename EQ>
    void HashTable<K,T,H,EQ>::analyze_table_occupancy() const
  {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Hashtable:\n";
    std::cout << "size                 : " << tab_size_ << "\n";
    std::cout << "# terms              : " << entries_ << "\n";
    std::cout << "load                 : " << (double(entries_)/double(tab_size_)*100) << " percent\n";

    int number_of_lists = 0;
    int max_number_of_list_entries = 0;
    for(long unsigned int i = 0; i < tab_size_; i++)
    {
      if( ht_[i] )
      {
        ++number_of_lists;
        int list_entries = 0;
        Node *current = ht_[i];
        while( current && ++list_entries ) current = current->next;
        if( max_number_of_list_entries < list_entries ) max_number_of_list_entries = list_entries;
      }
    }

    std::cout << "# used table-entries : " << number_of_lists << "\n";
    std::cout << "average list-length  : " << (double(entries_)/double(number_of_lists)) << "\n";
    std::cout << "maximum list-length  : " << max_number_of_list_entries << "\n";
    std::cout << "---------------------------------------------" << std::endl;
  }
}

#endif //Header

/*

  End of File: hashtable.h

*/
