/******************************************************************/
/* igluiInterval.h                                                */
/* -----------------------                                        */
/*                                                                */
/* Defines a 1D interval an assorted operations.                  */
/*                                                                */
/* Chris Wyman (12/06/2009)                                       */
/******************************************************************/


#ifndef IGLUI_INTERVAL
#define IGLUI_INTERVAL

#include "igluMath.h"

namespace iglu {

class IGLUInterval {
public:
	// Constructor and destructor.
	inline IGLUInterval( float iMin, float iMax );
	inline IGLUInterval( const IGLUInterval &iVal );
	inline ~IGLUInterval( void ) { }

	// Assignment constructors
	inline IGLUInterval& operator=(const IGLUInterval& iVal);

	// Check if this is really a valid interval
	inline bool IsValid( void ) const;

	// Check if the interval is degenerate (i.e., a point) 
	inline bool IsDegenerate( void ) const;

	// Is a value inside or outside the interval?
	inline bool IsInside ( float value ) const;
	inline bool IsOutside( float value ) const;

	// Do two intervals overlap?
	inline bool HasOverlap( const IGLUInterval &iVal ) const;

	// Get an overlapping interval
	inline IGLUInterval GetOverlap( const IGLUInterval &iVal );

	// Accessor methods
	inline float GetMin( void )	const;
	inline float GetMax( void )	const;

	// Modifiers
	inline void SetInterval( float nMin, float nMax );

	// A pointer to a IGLUInterval could have type IGLUInterval::Ptr
	typedef IGLUInterval *Ptr;

private:
	// The min/max along interval
	float iMin, iMax;
};



//////////////////////////////////////////////////////
//  Define the inline functions here...   
//////////////////////////////////////////////////////


inline IGLUInterval::IGLUInterval( float iMin, float iMax ) : 
	iMin(iMin), iMax(iMax) 
{ 
}

inline IGLUInterval::IGLUInterval( const IGLUInterval &iVal ) : 
	iMin(iVal.GetMin()), iMax(iVal.GetMax()) 
{ 
}

inline IGLUInterval& IGLUInterval::operator=(const IGLUInterval& iVal)
{
	iMin = iVal.iMin; 
	iMax = iVal.iMax; 
	return *this; 
}
	
inline bool IGLUInterval::IsValid( void ) const					
{ 
	return ( iMin <= iMax ); 
}

inline bool IGLUInterval::IsDegenerate( void ) const				
{ 
	return ( iMin == iMax ); 
}

inline bool IGLUInterval::IsInside ( float value ) const
{ 
	return ( value >= iMin ) && ( value <= iMax ); 
}

inline bool IGLUInterval::IsOutside( float value ) const			
{ 
	return ( value < iMin ) || ( value > iMax ); 
}

inline bool IGLUInterval::HasOverlap( const IGLUInterval &iVal ) const 
{ 
	return ( iVal.iMin <= iMax ) && ( iVal.iMax >= iMin ); 
}

inline IGLUInterval IGLUInterval::GetOverlap( const IGLUInterval &iVal ) 
{ 
	return IGLUInterval( max( iMin, iVal.iMin ), min( iMax, iVal.iMax ) ); 
}

inline float IGLUInterval::GetMin( void )	const					
{ 
	return iMin; 
}

inline float IGLUInterval::GetMax( void )	const					
{ 
	return iMax; 
}

inline void IGLUInterval::SetInterval( float nMin, float nMax )	
{ 
	iMin = nMin; 
	iMax = nMax; 
}


}

#endif