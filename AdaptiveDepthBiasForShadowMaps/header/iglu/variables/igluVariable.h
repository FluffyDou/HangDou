/******************************************************************/
/* igluVariable.h                                                 */
/* -----------------------                                        */
/*                                                                */
/* There's nearly always values in your program that need to      */
/*     change according to some user input.  Often these values   */
/*     require all sorts of logic to update, set, get, or         */
/*     change in response to user interface.                      */
/* This class tries to avoid spreading this logic out over the    */
/*     entire code, and instead encapsulate it all in a class.    */
/*     The additional advantage here, is that derived classes     */
/*     from IGLUVariable can all be treated identically, avoiding */
/*     special-casing depending on the variable type.             */
/*                                                                */
/* Realize that this class adds overhead, and shouldn't be used   */
/*     for all variables.  Just those that might vary with input  */
/*     or, perhaps, those defined in an input or script file.     */
/*                                                                */
/* Chris Wyman (02/09/2012)                                       */
/******************************************************************/


#ifndef IGLUVARIABLE_H
#define IGLUVARIABLE_H

#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "../igluArray1D.h"

#pragma warning( disable : 4996 )

namespace iglu {

//class IGLUCallback<T>;

typedef void (*CallbackPtr_VoidP)(void *);
typedef void (*CallbackPtr_VoidP_VoidP)(void *, void*);

class IGLUVariable {
public:
	IGLUVariable( char *name=0 ) : m_name(0), m_valueCallback(0), m_valueCallback2(0), m_callback2Data(0), m_isActive(true)
	                             { if (name) m_name = strdup(name); }
	virtual ~IGLUVariable()      
	{ 
		if(m_name) free(m_name); 
		for (uint i=0;i<m_callbacks.Size();i++) delete m_callbacks[i]; 
	}

	// This returns a true/false depending on if this variable 
	//     is setup to respond to the specified keystroke. 
	// PLEASE NOTE:  This is not necessarily a cheap call and
	//     likely duplicates work in UpdateFromInput()
	bool RespondsToInput( unsigned int inputKey );

	// This method updates the current variable, if appropriate
	virtual bool UpdateFromInput( unsigned int inputKey ) = 0;

	// Get/Set the variable name.  Beware, the returned name may be NULL.
	inline void SetName( const char *newName )	 { if(m_name) free(m_name); m_name = strdup(newName); }
	inline char *GetName( void )				 { return m_name; }


	// Add a callback to the list of callbacks that are executed when this
	//   variable changes.
	void AddCallback( IGLUCallback *callback )   { m_callbacks.Add( callback ); }

	// Setup callback(s) that are invoked whenever the variable value changes.
	//   -> Exactly one of each type may be setup.
	//   -> Both callbacks get a pointer to the IGLUVariable that changed.
	//   -> The second callback's 1st parameter is set to be 'firstInput'
	void SetCallback( CallbackPtr_VoidP callback );
	void SetCallback( CallbackPtr_VoidP_VoidP callback, void *firstInput );

	// Activate or deactivate all user-interface widgets associated with this variable
	void Activate( void );
	void Deactivate( void );
	inline bool IsActive( void ) const           { return m_isActive; }

	// This method sets the FLTK valuator (if one exists) that is tied to
	//    this UI varaible.  This allows settor method to update not only
	//    the variable's value, but also the displayed value in the FLTK widget.
	inline void AddFLTKValuator( void *fltkVal ) { m_fltkValuators.Add( fltkVal ); }

	// A pointer to a IGLUVariable could have type IGLUVariable::Ptr
	typedef IGLUVariable *Ptr;
protected:
	// Name of the variable, possibly NULL.  No effect on functioning of class
	char *m_name;	

	// Determines if widgets affiliated with this variable are active.  
	//   NOTE:  Inactive variables can still be changed and used inside C code,
	//          just not via the user interface widgets!
	bool m_isActive;

	// A list of keys this variable responds to.
	IGLUArray1D<unsigned int> m_keys;

	// This method is not good to use, as it doesn't really clean up.
	//     AddKeyResponse() is defined individually by each derived class,
	//     since the number of inputs may vary by variable type.
	// Using AddKeyResponse() & RemoveKeyResponse() repeatedly is
	//     asking for a huge memory leak.  (i.e., AddKeyReponse 
	//     increases the dynamic array size(s), but RemoveKey does not 
	//     reduce it.)  That's why this function is hidden.
	void RemoveKeyResponse( unsigned int key );

	// The FLTK widget that controls this variable (if it exists)
	IGLUArray1D<void *> m_fltkValuators;

	// A method that should be called when the value changes
	CallbackPtr_VoidP           m_valueCallback;
	CallbackPtr_VoidP_VoidP     m_valueCallback2;
	void                       *m_callback2Data;
	IGLUArray1D<IGLUCallback *> m_callbacks;

	// These methods actually call the callbacks appropriatly.  (XXX) Ideally, these
	//    protected methods should not be needed, but we're calling things from
	//    static members.  Evidently accessing member function pointers from static
	//    members is allowed, but sometimes accessing other member variables from
	//    the same static methods is *not* allowed.  This gets around this.
	void CallFirstCallback( void )    { if (m_valueCallback)  (*m_valueCallback)((void *)this ); }
	void CallSecondCallback( void )   { if (m_valueCallback2) (*m_valueCallback2)(m_callback2Data, (void *)this ); }

	void ExecuteCallbacks( void );
private:
	// Explicitly outlaw copy-constructors or UIVariable assignment
	//     by declaring undefined methods.
	IGLUVariable( const IGLUVariable& );
	IGLUVariable& operator=( const IGLUVariable& );
};


} // End namespace

// Define values used for input keyas.
#ifndef KEY_UNKNOWN 

// Key for all unknown keys
#define KEY_UNKNOWN         0

// Control characters (these are ASCII)
#define KEY_ESCAPE			27
#define KEY_TAB				9
#define KEY_RETURN			13
#define KEY_BACKSPACE       8
#define KEY_DELETE			127

// Special (usually non-ASCII keys)
#define KEY_UP_ARROW		1001
#define KEY_DOWN_ARROW		1002
#define KEY_LEFT_ARROW		1003
#define KEY_RIGHT_ARROW		1004
#define KEY_F1				1005
#define KEY_F2				1006
#define KEY_F3				1007
#define KEY_F4				1008
#define KEY_F5				1009
#define KEY_F6				1010
#define KEY_F7				1011
#define KEY_F8				1012
#define KEY_F9				1013
#define KEY_F10				1014
#define KEY_F11				1015
#define KEY_F12				1016
#define KEY_INSERT			1017
#define KEY_HOME			1019
#define KEY_END				1020
#define KEY_PAGE_UP			1021
#define KEY_PAGE_DOWN		1022

#define KEY_NO_MODIFIER			0x00000000
#define KEY_UNMODIFIED_MASK		0x0FFFFFFF
#define KEY_MODIFIER_MASK       0x70000000
#define KEY_MODIFIER_SHIFT		0x10000000
#define KEY_MODIFIER_CONTROL	0x20000000
#define KEY_MODIFIER_ALT		0x40000000

#define MOUSE_BASE               5000
#define MOUSE_NOBUTTON_POS_X	 5000
#define MOUSE_NOBUTTON_NEG_X	 5001
#define MOUSE_NOBUTTON_POS_Y	 5002
#define MOUSE_NOBUTTON_NEG_Y	 5003
#define MOUSE_LBUTTON_POS_X      5004
#define MOUSE_LBUTTON_NEG_X      5005
#define MOUSE_LBUTTON_POS_Y      5006
#define MOUSE_LBUTTON_NEG_Y      5007
#define MOUSE_MBUTTON_POS_X      5008
#define MOUSE_MBUTTON_NEG_X      5009
#define MOUSE_MBUTTON_POS_Y      5010
#define MOUSE_MBUTTON_NEG_Y      5011
#define MOUSE_RBUTTON_POS_X      5012
#define MOUSE_RBUTTON_NEG_X      5013
#define MOUSE_RBUTTON_POS_Y      5014
#define MOUSE_RBUTTON_NEG_Y      5015

#define MOUSE_LBUTTON            5025
#define MOUSE_MBUTTON            5026
#define MOUSE_RBUTTON            5027

#endif


#endif