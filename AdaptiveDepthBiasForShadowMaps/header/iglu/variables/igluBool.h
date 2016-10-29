/******************************************************************/
/* igluBool.h                                                     */
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


#ifndef IGLUBOOL_H
#define IGLUBOOL_H

#include <cstdio>
#include <cstdlib>
#include "igluVariable.h"

namespace iglu {

class IGLUBool : public IGLUVariable {
public:
	explicit IGLUBool()                           : IGLUVariable(0), m_value(false) {}
	explicit IGLUBool( bool value, char *name=0 ) : IGLUVariable(name), m_value(value) {}
	virtual ~IGLUBool() {}

	// If you used a default constructor, you can initialize using this method
	void Initialize( bool value, char *name=0 );

	// Get and set values
	//     Note:  We can get/set values with typical class accessor methods
	void SetValue( bool newValue, bool executeCallback=true );
	inline bool GetValue( void ) const					{ return m_value; }
	inline bool *GetValuePtr( void )					{ return &m_value; }


	//     Note:  We can also use UIBool just like a standard bool, thanks
	//            the the magic of overloaded assignment and casting operators
	inline IGLUBool& operator=( bool newValue )			{ SetValue(newValue); return *this; }
	inline operator const bool()						{ return m_value; }

	// This method updates the current variable, if appropriate
	virtual bool UpdateFromInput( unsigned int inputKey );

	// This method adds a key for the variable to respond to
	inline void AddKeyResponse( unsigned int key )		{ m_keys.Add( key ); }

	// This method is a callback that allows FLTK widget to modify IGLUBool values
	static void FLTKCallback( Fl_Widget *, void * );

	// This method is a callback for a button bound to this IGLUBool.  There is no
	//     equivalent in other IGLUVariable types.  A button sets the IGLUBool to
	//     true, with the expectation that it will be reset to false within the 
	//     program (not by the user).
	static void FLTKButtonCallback( Fl_Widget *, void * );

	// A pointer to a IGLUBool could have type IGLUBool::Ptr
	typedef IGLUBool *Ptr;
protected:
	bool m_value;

private:
	// Explicitly outlaw copy-constructors or UIVariable assignment
	//     by declaring undefined methods.
	IGLUBool( const IGLUBool& );
	IGLUBool& operator=( const IGLUBool& );
};

} // end iglu namespace

#endif
