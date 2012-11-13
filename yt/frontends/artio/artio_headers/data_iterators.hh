#ifndef __DATA_ITERATORS_HH
#define __DATA_ITERATORS_HH

#include <iterator>
/*template<typename _Container,typename _InputIterator>
class conditional_iterator
: public _Container::iterator
{
	
	
};*/

//! TBD
template<typename _Container,typename _InputIterator>
class conditional_back_insert_iterator
: public std::iterator<std::output_iterator_tag, void, void, void, void>
{
	protected:
		_Container* container;			//!< TBD
		_InputIterator cond_iterator;	//!< TBD
	
	public:
		//! TBD
		typedef _Container		container_type;
		
		//! TBD
		conditional_back_insert_iterator(_Container& __x,const _InputIterator& __i, int __prealloc=0) 
		:	container(&__x), cond_iterator(__i) 
		{ 
			if( __prealloc > 0 )
				__x.reserve( __x.size()+__prealloc);
		}
		
		/*conditional_back_insert_iterator( const conditional_back_insert_iterator& __it )
		: container( __it.container ), cond_iterator( __it.cond_iterator )
		{ }*/
		
		//! TBD
		conditional_back_insert_iterator&
		operator=( typename _Container::const_reference __value )
		{
			if( *cond_iterator )
				container->push_back(__value);
			++cond_iterator;
			return *this;
		}
		
		//! TBD
		conditional_back_insert_iterator&
		operator*()
		{ return *this; }
		
		//! TBD
		conditional_back_insert_iterator&
		operator++()
		{ return *this; }
		
		//! TBD
		conditional_back_insert_iterator&
		operator++(int)
		{ return *this; }
};


//! TBD
template<typename _Container, typename _InputIterator>
inline conditional_back_insert_iterator<_Container, _InputIterator>
conditional_back_inserter(_Container& __x, const _InputIterator& __i, int __prealloc=0)
{ return conditional_back_insert_iterator<_Container,_InputIterator>(__x,__i,__prealloc); }




#endif
