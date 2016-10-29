/*****************************************************************************
** igluCallback.h                                                           **
** -------------                                                            **
**                                                                          **
** A 'callback' class to store state needed to update IGLU objects based on **
**    updates to IGLU variables by the user-interface (e.g., a FLTK slider) **
**                                                                          **
** The callback class is templated.  The template variable 'T' should be    **
**    one a class type derived from IGLUVariable (e.g., IGLUInt).  The idea **
**    is that when this particular variable changes values, some function   **
**    or class method will be called to handle the results of the update.   **
**    In order to handle the update, it needs to know about the variable.   **
**    This means, whatever callback function is called, it's last parameter **
**    will be a pointer to the modified IGLUVariable.                       **
**                                                                          **
** Unfortunately, 
**    
**                                                                          **
** Chris Wyman (4/10/2012)                                                  **
*****************************************************************************/

#ifndef __IGLU_CALLBACK_H
#define __IGLU_CALLBACK_H


namespace iglu {

// We'll define this later
class IGLUVariable;


//////////////////////////////////////////////////////////////////////////////
// All callbacks derive from this abstract class, and must provide an 
//    Execute() method.  The execute method will be called with the 
//    IGLUVariable that has been modified.

class IGLUCallback
{
public:
	IGLUCallback()  {}
	virtual ~IGLUCallback() {}

	// Call Execute() to run the selected callback.
	virtual void Execute( IGLUVariable *data ) const = 0;
};



//////////////////////////////////////////////////////////////////////////////
// Simple callback, in case you have a need to respond to a variable change
//     but don't need its value (e.g., essentially respond by 'you changed
//     the value' without needing the value)

class IGLUSimpleCallback : public IGLUCallback
{
public:
	IGLUSimpleCallback( void (*fnPtr)(void) ) : IGLUCallback(), m_callbackFn(fnPtr) {}
	virtual ~IGLUSimpleCallback() {}

	// Call our callback
	virtual void Execute( IGLUVariable * ) const { (*m_callbackFn)(); }

private:
	void (*m_callbackFn)(void);
};


//////////////////////////////////////////////////////////////////////////////
// Simple callback that passes a specified bit of fixed data to the callback.
//    

class IGLUStaticDataCallback : public IGLUCallback
{
public:
	IGLUStaticDataCallback( void (*fnPtr)(int), int val ) : IGLUCallback(), m_callbackFn(fnPtr), m_val(val) {}
	virtual ~IGLUStaticDataCallback() {}

	// Call our callback
	virtual void Execute( IGLUVariable * ) const { (*m_callbackFn)( m_val ); }

private:
	int m_val;
	void (*m_callbackFn)(int);
};

//////////////////////////////////////////////////////////////////////////////
// Callback that requires only the new value of the IGLUVariable to handle.
//    

class IGLUVariableCallback : public IGLUCallback
{
public:
	IGLUVariableCallback( void (*fnPtr)(IGLUVariable *) ) : IGLUCallback(), m_callbackFn(fnPtr) {}
	virtual ~IGLUVariableCallback() {}

	// Call our callback
	virtual void Execute( IGLUVariable *data ) const { (*m_callbackFn)( data ); }

private:
	void (*m_callbackFn)(IGLUVariable *);
};

//////////////////////////////////////////////////////////////////////////////
// A callback that executes a member function from a UpdateClass object.
//     This particular member function should take a single input, a pointer
//     to the IGLUVariable that changed.

template <class UpdateClass>
class IGLUMemberCallback : public IGLUCallback 
{
public:
	IGLUMemberCallback( UpdateClass &updateObj, void (UpdateClass::*fnPtr)(IGLUVariable *) ) 
		: IGLUCallback(), m_callbackFn(fnPtr), m_objToUpdate(updateObj) {}
	virtual ~IGLUMemberCallback() {}

	// Call our callback
	virtual void Execute( IGLUVariable *data ) const { (m_objToUpdate.*m_callbackFn)( data ); }

private:
	UpdateClass &m_objToUpdate;
	void (UpdateClass::*m_callbackFn)(IGLUVariable *);
};






} // end namespace


#endif
