/******************************************************************/
/* igluInt.h                                                      */
/* -----------------------                                        */
/*                                                                */
/* There's nearly always values in your program that need to      */
/*     change according to some user input.  Often these values   */
/*     require all sorts of logic to update, set, get, or         */
/*     change in response to user interface.                      */
/* This class, and others like it, inherit from UIVariable and    */
/*     try to avoid spreading UI interface logic out over the     */
/*     entire code, and instead encapsulate it in these classes.  */
/*     The additional advantage here, is that derived classes     */
/*     from UIVariable can all be treated identically, avoiding   */
/*     special-casing depending on the variable type.             */
/*                                                                */
/* Realize that this class adds overhead, and shouldn't be used   */
/*     for all variables.  Just those that might vary with input  */
/*     or, perhaps, those defined in an input or script file.     */
/*                                                                */
/* Chris Wyman (02/09/2012)                                       */
/******************************************************************/


#ifndef IGLUINT_H
#define IGLUINT_H

#include <cstdio>
#include <cstdlib>
#include "igluVariable.h"
#include "iglu/igluRange.h"

class Fl_Widget;

namespace iglu {

class IGLUInt : public IGLUVariable {
public:
	explicit IGLUInt();
	explicit IGLUInt( int value, int _min=INT_MIN, int _max=INT_MAX, char *name=0 );
	explicit IGLUInt( int value, const IGLURange<int> &range, int delta=1, char *name=0 );
	virtual ~IGLUInt() {}

	// If you used a default constructor, you can initialize using this method
	void Initialize( int value, const IGLURange<int> &range, int delta=1, char *name=0 );

	// Get and set values
	//     Note:  We can get/set values with typical class accessor methods
	void SetValue( int newValue, bool force = false, bool executeCallback=true );
	inline int GetValue( void ) const								{ return m_value; }
	inline int *GetValuePtr( void )									{ return &m_value; }

	// Get the range for the UI variable
	inline int GetMaxValue( void ) const                            { return m_max; }
	inline int GetMinValue( void ) const                            { return m_min; }

	// Get delta value.  This is NOT used in the traditional keyboard handler, but
	//   only in conjunction with FLTK widgets
	inline int GetDelta( void ) const                               { return m_delta; }

	// Set the max and min values
	inline void SetValueRange( int newMin, int newMax )				{ m_max=newMax; m_min=newMin; m_value=iglu::clamp(m_value, m_min, m_max);}

	//     Note:  We can also use UIInt just like a standard bool, thanks
	//            the the magic of overloaded assignment and casting operators
	inline IGLUInt& operator=( int newValue )						{ SetValue( newValue ); return *this; }
	inline operator const int()										{ return m_value; }

	// This method updates the current variable, if appropriate
	virtual bool UpdateFromInput( unsigned int inputKey );

	// This method adds a key for the variable to respond to
	void AddKeyResponse( unsigned int key, int response );

	// This method is a callback that allows FLTK widget to modify IGLUInt values
	static void FLTKCallback( Fl_Widget *, void * );

	// A pointer to a IGLUInt could have type IGLUInt::Ptr
	typedef IGLUInt *Ptr;
protected:
	int m_value, m_min, m_max, m_delta;
	IGLUArray1D<int> m_responses;

private:
	// Explicitly outlaw copy-constructors or UIVariable assignment
	//     by declaring undefined methods.
	IGLUInt( const IGLUInt& );
	IGLUInt& operator=( const IGLUInt& );
};

}// end namespace

#endif
