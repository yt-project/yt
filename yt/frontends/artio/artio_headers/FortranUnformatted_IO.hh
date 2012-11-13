/*
	FortranUnformatted_IO.hh
	This file contains a C++ class for access to FORTRAN unformatted files
			
	Copyright (C) 2008  Oliver Hahn, ojha@gmx.de

    It is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __FORTRAN_UNFORMATTED_HH
#define __FORTRAN_UNFORMATTED_HH

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <iterator>

#define DEFAULT_ADDLEN 4

//! A class to perform IO on FORTRAN unformatted files
/*! FortranUnformatted provides sequential read access to FORTRAN
    unformatted files.
 */
class FortranUnformatted{
	
protected:
	std::string   m_filename;		//!< the file name 
	std::fstream  m_ifs;			//!< STL ifstream object
	int           m_addlen;			//!< number of bytes in pre-/suffix data of FORTRAN unformatted data
		
public:
	
	//! constructor for FortranUnformatted
	/*! simple constructor for FortranUnformatted
	 * @param filename the name  of the FORTRAN unformatted file to be opened for IO
	 * @param mode a combination of std::ios_base::openmode determining the mode in which the files are openend
	 * @param addlen the number of bytes which are pre- & postpended to unformatted arrays giving their size (default=4)
	 */
	explicit FortranUnformatted( std::string filename, std::ios_base::openmode mode = std::ios_base::in, int addlen=DEFAULT_ADDLEN )
		: m_filename( filename ), 
			m_ifs( m_filename.c_str(), mode|std::ios::binary ), 
			m_addlen( addlen )
	{ 
		if(!m_ifs.good() || !m_ifs.is_open() )
			throw std::runtime_error("FortranUnformatted : unable to open file \'"
						+m_filename+"\'for read access");
		m_ifs.exceptions ( std::fstream::eofbit | std::fstream::failbit 
						| std::fstream::badbit );
	}
	
	//! read data from FORTRAN unformatted file
	/*! read a scalar from FORTRAN unformatted file. If the data is not scalar 
		(i.e. an array) or the function is called with a template parameter of 
		a different size than the stored data, a std::runtime_error is thrown.
	 \param r reference to the return value
	 */
	template< typename T > void read( T& r )
	{
		unsigned n1,n2;
		
		try{
			m_ifs.read( (char*)&n1, m_addlen );
			m_ifs.read( (char*)&r, sizeof(T) );
			m_ifs.read( (char*)&n2, m_addlen );
		}catch(...){
			throw;
		}

		if( n1 != sizeof(T) || n1 != n2 ){
			throw std::runtime_error("FortranUnformatted::read : invalid type"\
						" conversion when reading FORTRAN unformatted.");
		}
	}
	
	//! write a single variable, arbitrary type to the Fortran unformatted file
	/*!
	 * @param w the variable, arbitrary type
	 */
	template< typename T > void write( T& w )
	{
		unsigned n1 = sizeof(T);
		
		try{
			m_ifs.write( (char*)&n1, m_addlen );
			m_ifs.write( (char*)&w, n1 );
			m_ifs.write( (char*)&n1, m_addlen );
		}catch(...){
			throw;
		}
	}
	
	
	//! write a range of a container data object given by iterators
	/*! 
	 * @param begin iterator pointing to the beginning of the range
	 * @param end iterator pointing to the end of the range
	 */
	template< typename _InputIterator >
	void write( _InputIterator begin, _InputIterator end )
	{
		_InputIterator it(begin);
		unsigned nelem = std::distance(begin,end);
		unsigned sz = sizeof(*begin);
		unsigned totsz = sz*nelem;
		
		try{
			m_ifs.write( (char*)&totsz, m_addlen );
			while( it!=end ){
				m_ifs.write( (char*)&(*it), sz );
				++it;
			}
			m_ifs.write( (char*)&totsz, m_addlen );
		}catch(...){
			throw;
		}
	}
	
	//! read masked data
	/*! this function reads data but discards all those for which the iterator
	 *  mask is not set to a 'true' bit. It is necessary that the mask iterator
	 *  can be increased the same number of times as there is data to be read.
	 * @param mask an iterator which is increased for each element read from the file and determines rejection or not
	 * @param data an output iterator to which the data is written
	 */
	template< typename basetype, typename _InputIterator, typename _OutputIterator >
	_InputIterator read( _InputIterator mask, _OutputIterator data )
	{
		_InputIterator oldmask = mask;
		std::vector<basetype> temp;
		typename std::vector<basetype>::iterator temp_it;
		
		unsigned n1,n2;
		try{
			m_ifs.read( (char*)&n1, m_addlen );
			temp.resize(n1/sizeof(basetype));
			m_ifs.read( (char*)&temp[0], n1 );
			
			for( temp_it = temp.begin(); temp_it!=temp.end(); ++temp_it ){
				//... copy data if masked - this also performs a type conversion if needed ...//
				if( *mask )
					*data = *temp_it;
				oldmask = mask++;
				if( mask == oldmask ) break;
			}
			//std::copy(temp.begin(),temp.end(), data);
			
			m_ifs.read( (char*)&n2, m_addlen );
		
		}catch(...){
			throw std::runtime_error("FortranUnformatted::read : error reading FORTRAN unformatted.");
		}
		if( n1!=n2 )
			throw std::runtime_error("FortranUnformatted::read : invalid type conversion when"\
						" reading FORTRAN unformatted.");
		
		return mask;
	}
	
	
	//! read data
	/*! this function reads data from a Fortran unformatted file and writes it to
	 *  an output iterator.
	 * @param data an output iterator to which the data is written
	 */
	template< typename basetype, typename _OutputIterator >
	_OutputIterator read( _OutputIterator data )
	{
		std::vector<basetype> temp;
		typename std::vector<basetype> temp_it;
		
		unsigned n1,n2;
		try{
			m_ifs.read( (char*)&n1, m_addlen );
			temp.resize(n1/sizeof(basetype));
			m_ifs.read( (char*)&temp[0], n1 );
			//... copy data - this also performs a type conversion if needed ...//
			std::copy(temp.begin(),temp.end(), data);
			m_ifs.read( (char*)&n2, m_addlen );
		
		}catch(...){
			throw std::runtime_error("FortranUnformatted::read : error reading FORTRAN unformatted.");
		}
		if( n1!=n2 )
			throw std::runtime_error("FortranUnformatted::read : invalid type conversion when"\
						" reading FORTRAN unformatted.");
		
		
		return data;
	}
	
	//! check if beyond end-of-file
	bool eof( void )
	{
		bool bcheck;
		try{
			bcheck = (m_ifs.peek()==EOF);
		}catch(...){
			m_ifs.clear();
			return true;
		}
		return bcheck;
	}
		
	//! skip ahead in FORTRAN unformatted file
	/*! skips n datasets ahead in FORTRAN unformatted file. Equivalent to 
	 *	n READ(X) without argument in FORTRAN
	 *  @param n number of entries to skip ahead (scalar or arrays)
	 */
	void skip_n( unsigned n )
	{
		unsigned ndone = 0;
		while( ndone < n ){
			int n1, n2;
			try{
				m_ifs.read( (char*)&n1, m_addlen );
				m_ifs.seekg( n1, std::ios::cur );
				m_ifs.read( (char*)&n2, m_addlen );
			}catch(...){
				throw std::runtime_error("FortranUnformatted::skip_n : error seeking in FORTRAN unformatted file.");
			}
			++ndone;
		}
	}
	
	//! just a std::streampos
	typedef std::streampos streampos;
	
	//! wrapper for the std::ifstream::tellg() function
	streampos tellg( void )
	{ return m_ifs.tellg(); }
	
	//! wrapper for the std::ifstream::seekg() function
	void seekg( streampos pos )
	{ m_ifs.seekg( pos ); /*, std::ios::beg );*/ }
	
	
	//! jump to dataset in FORTRAN unformatted files
	/*! moves past the nth WRITE(X) operation from the beginning in a FORTRAN 
	 *	unformatted file 
	 * @param n number of entries to skip (scalar or arrays)
	 * @sa skip_n()
	 */
	void skip_n_from_start( unsigned n )
	{
		m_ifs.seekg(0,std::ios::beg);
		skip_n( n );
	}
	
	//! skip ahead in FORTRAN unformatted file
	/*! skips to next dataset in FORTRAN unformatted file. Equivalent to a 
	 * READ(X) without argument in FORTRAN
	 * @sa skip_n
	 */
	void skip( void )
	{
		skip_n(1);
	}

	//! rewind file
	/*! rewinds the file to the beginning 
	 */
	void rewind( void )
	{
		m_ifs.seekg(0,std::ios::beg);
	}

	
	
	/*template< typename T >
	friend FortranUnformatted& operator>> (FortranUnformatted &is, T& data);

	template< typename T >
	friend FortranUnformatted& operator>> (FortranUnformatted &is, std::vector<T>& data);*/
};



#endif //__FORTRAN_UNFORMATTED_HH
