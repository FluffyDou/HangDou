/******************************************************************/
/* igluFloat.h	                                                  */
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


#ifndef IGLUFLOAT_H
#define IGLUFLOAT_H

#include <cstdio>
#include <cstdlib>
#include "igluVariable.h"
#include "iglu/igluRange.h"

class Fl_Widget;

namespace iglu {

class IGLUFloat : public IGLUVariable {
public:
	explicit IGLUFloat();
	explicit IGLUFloat( float value, float _min=FLT_MIN, float _max=FLT_MAX, char *name=0 );
	explicit IGLUFloat( float value, const IGLURange<float> &range, float delta=0.1, char *name=0 );
	virtual ~IGLUFloat() {}

	// If you used a default constructor, you can initialize using this method
	void Initialize( float value, const IGLURange<float> &range, float delta=0.1, char *name=0 );

	// Get and set values
	//     Note:  We can get/set values with typical class accessor methods
	void SetValue( float newValue, bool force = false, bool executeCallback=true );							
	inline float GetValue( void ) const								{ return m_value; }
	inline float *GetValuePtr( void )								{ return &m_value; }

	// Get the range for the UI variable
	inline float GetMaxValue( void ) const                          { return m_max; }
	inline float GetMinValue( void ) const                          { return m_min; }

	// Get delta value.  This is NOT used in the traditional keyboard handler, but
	//   only in conjunction with FLTK widgets
	inline float GetDelta( void ) const                             { return m_delta; }

	// Set the max and min values
	inline void SetValueRange( float newMin, float newMax )			{ m_max=newMax; m_min=newMin; m_value=iglu::clamp(m_value, m_min, m_max); }

	//     Note:  We can also use UIInt just like a standard bool, thanks
	//            the the magic of overloaded assignment and casting operators
	inline IGLUFloat& operator=( float newValue )					{ SetValue( newValue ); return *this; }
	inline operator const float()									{ return m_value; }

	// This method updates the current variable, if appropriate
	virtual bool UpdateFromInput( unsigned int inputKey );

	// This method adds a key for the variable to respond to
	inline void AddKeyResponse( unsigned int key, float response )	{ m_keys.Add( key ); m_responses.Add( response ); }

	// This method is a callback that allows FLTK widget to modify IGLUFloat values
	static void FLTKCallback( Fl_Widget *, void * );

	// A pointer to a IGLUFloat could have type IGLUFloat::Ptr
	typedef IGLUFloat *Ptr;
protected:
	float m_value, m_max, m_min, m_delta;
	IGLUArray1D<float> m_responses;

private:
	// Explicitly outlaw copy-constructors or UIVariable assignment
	//     by declaring undefined methods.
	IGLUFloat( const IGLUFloat& );
	IGLUFloat& operator=( const IGLUFloat& );
};


} // end iglu namespace

#endif
