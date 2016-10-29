/******************************************************************/
/* igluRange.h                                                    */
/* -----------------------                                        */
/*                                                                */
/* Defines a simple 1D templated range class.  This is similar to */
/*    IGLUInterval, except it is templated, and doesn't come with */
/*    as many associated methods.  Mainly it's an encapsulation   */
/*    for passing ranges around.                                  */
/*                                                                */
/* Chris Wyman (12/06/2009)                                       */
/******************************************************************/


#ifndef IGLU_RANGE
#define IGLU_RANGE

#include "igluMath.h"

namespace iglu {

template<class T> 
class IGLURange {
public:
	// Constructor and destructor.
	IGLURange( T iMin, T iMax );
	IGLURange( const IGLURange<T> &iVal );
	~IGLURange( void ) { }

	// Assignment constructors
	inline IGLURange<T>& operator=(const IGLURange<T>& iVal);

	// Accessor methods
	T GetMin( void ) const	{ return m_min; }
	T GetMax( void ) const	{ return m_max; }

	// A pointer to a IGLURange<T> could have type IGLURange<T>::Ptr
	typedef IGLURange<T> *Ptr;

private:
	// The min/max along interval
	T m_min, m_max;
};



//////////////////////////////////////////////////////
//  Define the inline functions here...   
//////////////////////////////////////////////////////

template<class T> 
inline IGLURange<T>::IGLURange( T iMin, T iMax ) : 
	m_min(iMin), m_max(iMax) 
{ 
}

template<class T> 
inline IGLURange<T>::IGLURange( const IGLURange<T> &iVal ) : 
	m_min(iVal.GetMin()), m_max(iVal.GetMax()) 
{ 
}

template<class T> 
inline IGLURange<T>& IGLURange<T>::operator=(const IGLURange<T>& iVal)
{
	m_min = iVal.m_min; 
	m_max = iVal.m_max; 
	return *this; 
}
	

}

#endif