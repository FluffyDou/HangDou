/****************************************************************************/
/* igluCountedPtr.h                                                         */
/* -----------------------                                                  */
/*                                                                          */
/* This is an IGLU helper file that defines a number of smart-pointer types */
/*    used throughout IGLU.                                                 */
/*                                                                          */
/* PLEASE do not change code inside this file unless you have a good grasp  */
/*    of how smart pointers work and are implemented.  There are a fair     */
/*    number of online resources discussing smart pointers, and the classes */
/*    IGLUCountedPtr<> and IGLUWrapperCountedPtr<> are (respectively) the   */
/*    developer-visible smart pointer type and an internal wrapper that     */
/*    should not be directly used in code.  While it is always hard to say  */
/*    code is 100% bug free, these two classes are straightforward smart    */
/*    pointer classes and really should have few (if any) bugs and should   */
/*    not need any new functionality.  If you need to change those two      */
/*    classes, think again (and again) before changing anything!  And talk  */
/*    it over with someone before modify a line in this file!               */
/*                                                                          */
/* The third class is a specialized smart pointer that allows a couple      */
/*    array operators [] to be applied directly to the pointer type rather  */
/*    than first dereferencing it.  This code is more likely to have one or */
/*    two obscure logic bugs, but if you think that is the case, please     */
/*    talk it over first.  (Also the third class could feasibly be moved    */
/*    out of this file eventually, as it is somewhat specialized and other, */
/*    different variations may be desirable.)                               */
/*                                                                          */
/* If you don't believe me about not changing the code on a whim, this file */
/*    has taken a good 5-10 days of my life to get right and debug.  Many   */
/*    minor changes, if you don't understand them fully, will introduce     */
/*    immediate compilation errors and/or segfaults throughout all of IGLU  */
/*    or will cause mysterious segfaults only in strange usage cases.       */
/*                                                                          */
/* Chris Wyman (01/12/2011)                                                 */
/****************************************************************************/


#ifndef IGLU_COUNTED_PTR_H
#define IGLU_COUNTED_PTR_H

#include <assert.h>

namespace iglu {


// This is the wrapper for a traditional pointer that allows us to create a smart/counter pointer
//     and do garbage collection once all reference to it have gone away.  You should never use
//     on directly in most developer code.
// A IGLUWrapperCountedPtr<T>* is the pointer that smart pointers actually pass around.  That way
//     when the smart pointer in question is deleted, we don't actually delete the memory, we just 
//     delete one reference to this class, which stores the real pointer.
template <class T>
class IGLUWrapperCountedPtr {
public:
	// Copy constructor from an original object.  Assuming non-NULL input pointer, we have a reference of one
	IGLUWrapperCountedPtr( T *origPtr ) : m_ptr(origPtr)   { m_count = (origPtr ? 1 : 0); }

	// When this is deleted, by design, m_count and m_ptr should both be 0.  If not, there's a bug in this class.
	~IGLUWrapperCountedPtr()                               { assert(m_count==0 && m_ptr==0); }

	// Somebody made a new copy to this object.  Remember that someone else knows about us.
	inline unsigned int AddReference( void )  { return ++m_count; }

	// Somebody who pointed to us is going away.  Decrease the number of people who know about us.
	//    If that was the last reference to this object, we should delete our internal memory
	inline unsigned int FreeReference( void ) { assert(m_count > 0); 
	                                            if (m_count == 1 && m_ptr) { delete m_ptr; m_ptr = 0; }
	                                            return --m_count; }

	// Is the pointer we are wrapped around NULL?  That would be silly.  But sillier things have happened.
	inline bool IsNull( void ) const          { return (m_ptr==0 || m_count==0); }

	// OK, they want to actually use this pointer to DO something.  Since they don't know about the real
	//    pointer, we need to access it for the -> or * operators for them.
	inline T*   operator->()             { assert(m_ptr); return m_ptr; }
	inline T&   operator*(void)          { assert(m_ptr); return *m_ptr; }

	// Compare two wrapped pointers to see if they point to the same thing
	bool operator==(const IGLUWrapperCountedPtr &cmp) const  { return (m_ptr == cmp.m_ptr); }
private:
	T             *m_ptr;
	unsigned  int  m_count;
};



// This is the smart pointer class.  
//    You can use "IGLUCountedPtr<T> myPtr;" exactly as you could "T *myPtr;"
//    However, this smart pointer will garbage collect after all references are gone.
template <class T>
class IGLUCountedPtr {
public:
	// Default constructor
	IGLUCountedPtr() : m_wrapPtr(0) {}

	// Constructor from an original object.  Unless input was NULL, we now have a reference of one in our counted object
	IGLUCountedPtr( T *origPtr ) : m_wrapPtr(0)       { if (origPtr) m_wrapPtr = new IGLUWrapperCountedPtr<T>( origPtr ); }

	// Copy constructor from an existing wrapped pointer.  If the copied pointer is non-null, add
	//    one to its count, then make a copy of its pointer.
	IGLUCountedPtr( const IGLUCountedPtr &copy )      { if (copy.m_wrapPtr) copy.m_wrapPtr->AddReference(); 
	                                                     m_wrapPtr = copy.m_wrapPtr; }

	// Copy operator using assignment.  Check if the two sides of the assignment are already the same.
	//    If not, then if we *had* a pointer, we lose one reference to that so decrement its reference.
	//    Then the pointer on the right hand side of the assignment gets an additional reference.
	void operator= ( const IGLUCountedPtr &copy )	  { if (m_wrapPtr && copy.m_wrapPtr && ( (*m_wrapPtr) == *(copy.m_wrapPtr) ) ) return;
		                                                RemoveReference(); 
	                                                    if (copy.m_wrapPtr) copy.m_wrapPtr->AddReference();
														m_wrapPtr = copy.m_wrapPtr; }

	// Copy operator using assignment.  Check if the two sides of the assignment are already the same.
	//    If not, then if we *had* a pointer, we lose one reference to that so decrement its reference.
	//    Then the pointer on the right hand side of the assignment gets an additional reference.
	void operator= ( T *origPtr )	                  { RemoveReference(); 
	                                                    if (origPtr)   m_wrapPtr = new IGLUWrapperCountedPtr<T>( origPtr ); }

	// OK.  We're deleting our pointer, so we have one less reference to what it points at.
	~IGLUCountedPtr()                                 { RemoveReference(); } 

	// Perhaps we want to check if the pointer is null
	inline bool IsNull( void ) const                  { return (m_wrapPtr==0 || m_wrapPtr->IsNull()); }

	// Some operator overloading for comparisons
	inline bool operator==( const IGLUCountedPtr<T> &cmp ) const { return (m_wrapPtr == cmp.m_wrapPtr); }
	inline bool operator==( const T *cmp ) const                 { return (m_wrapPtr == cmp); }
	inline bool operator!=( const IGLUCountedPtr<T> &cmp ) const { return (m_wrapPtr != cmp.m_wrapPtr); }
	inline bool operator!=( const T *cmp ) const                 { return (m_wrapPtr != cmp); }

	// Ok.  Now if we actually want to use this, we need dereference operators:
	inline T*   operator->()                    { assert(m_wrapPtr); return m_wrapPtr->operator->(); }
	inline T&   operator*(void)                 { assert(m_wrapPtr); return m_wrapPtr->operator*(); }
private:
	// The wrapper around the real pointer that counts the number of active references.
	IGLUWrapperCountedPtr<T> *m_wrapPtr;

	// An internal utility that removes our current reference to the wrapped pointer.  We use this a 
	//    couple times, and it's a bit tricky.  If our pointer is NULL, do nothing.  If we have
	//    an object referenced, reduce the references by 1.  If we were the last person looking at it,
	//    we need to actually go ahead and delete the thing.
	inline void RemoveReference( void )         { if (m_wrapPtr) { unsigned int count = m_wrapPtr->FreeReference();
	                                                               if (count==0) delete m_wrapPtr; } }

};




// This is a fancy smart pointer that allows our odd syntactic sugar that applies a [] operator to a pointer.
// How it works:  
//     Design a class T that can have two operators that return types 'type1' and 'type2'
//     Declare "IGLUPtr<T,type1,type2> myPtr = 0;"             (The '= 0' is optional)
//     Assign values "myPtr = new T();"                        (Assigning directly from base type T)
//     Assign values "IGLUPtr<T,type1,type2> myPtr2 = myPtr;"  (Assign from another pointer)
//     Use as real pointers (e.g., myPtr->Func() accesses the Func() method in T)
//     Use an array operator (e.g., myPtr[const char *] accesses an overloaded [const char*] operator in T returning type1)
//     Use an array operator (e.g., myPtr[int] accesses an overloaded [int] operator in T returning type2)
//     If you don't have corresponding [] operators, that template type does not matter (as long as the
//         user doesn't try to *call* it, in which case he'll get a boatload of templating errors.
//     Ideally this wrapper class would be specialized to each IGLU class with appropriate functionality,
//         However, this is currently not high on my priority list since this works pretty well.
//
//     The template pointer mess is mostly hidden from the user by IGLU's per-class IGLUClass::Ptr typedef,
//         so the developer actually doesn't even need to know what type of pointer is being used, a "real"
//         one or some fancy smart one. 
template <class T, class ArrRet1, class ArrRet2>
class IGLUPtr {
public:
	// Default constructor
	IGLUPtr() : m_ptr(), m_cheaterPtr(0)                         {}

	// Constructor from an existing pointer
	IGLUPtr( T *origPtr ) : m_ptr(), m_cheaterPtr(0)             { m_ptr = origPtr; 
	                                                               m_cheaterPtr = IsNull() ? 0 : m_ptr.operator->(); }

	// Copy constructor from an existing pointer. 
	IGLUPtr( const IGLUPtr<T,ArrRet1,ArrRet2> &copy ) : m_ptr(copy.m_ptr) { m_cheaterPtr = IsNull() ? 0 : m_ptr.operator->(); } 

	// Copy operator using assignment.  
	void operator= ( const IGLUPtr<T,ArrRet1,ArrRet2> &copy )    { m_ptr = copy.m_ptr; 
	                                                               m_cheaterPtr = IsNull() ? 0 : m_ptr.operator->(); }	  

	// Copy operator using assignment. 
	void operator= ( T *copyPtr )	                             { m_ptr = copyPtr; 
	                                                               m_cheaterPtr = IsNull() ? 0 : m_ptr.operator->(); }

	// OK.  We're deleting our IGLUPtr (probably because it left scope).  This destructor
	//     automatically calls m_ptr.~IGLUCounterPtr(), so we don't need to delete m_ptr.
	//     And m_cheaterPtr is *only* a reference to the internals of m_ptr, so leave it alone!
	~IGLUPtr()                                                   {} 

	// Perhaps we want to check if the internal pointer is null
	inline bool IsNull( void ) const                             { return m_ptr.IsNull(); }

	// Some operator overloading for comparisons
	inline bool operator==( const IGLUPtr<T,ArrRet1,ArrRet2> &cmp ) const { return (m_ptr == cmp.m_ptr); }
	inline bool operator==( const T *cmp ) const                          { return (m_ptr == cmp); }
	inline bool operator!=( const IGLUPtr<T,ArrRet1,ArrRet2> &cmp ) const { return (m_ptr != cmp.m_ptr); }
	inline bool operator!=( const T *cmp ) const                          { return (m_ptr != cmp); }

	// Ok.  Now if we actually want to use this fancy pointer, we'll need dereference operators:
	//    In all cases, if the internal pointer is NULL, we'll fail an assertion.  If you're reading
	//    this to debug an error, that means your pointer is NULL and you're trying to dereference it.
	//    (Generally, this is not an error in IGLU!  Fix your pointer.)
	inline T*        operator->()                                { assert(m_cheaterPtr); return m_cheaterPtr; }
	inline T&        operator*(void)                             { assert(m_cheaterPtr); return *m_cheaterPtr; }
	inline ArrRet2   operator[](int var)                         { assert(m_cheaterPtr); return m_cheaterPtr->operator[](var); }
	inline ArrRet1   operator[](const char *var)                 { assert(m_cheaterPtr); return m_cheaterPtr->operator[](var); }
	
private:
	// This it the glorified pointer to the user's data structure.  This is a counted smart pointer,
	//       allowing IGLUPtr<> to be used like an actual, real pointer...  except when you go applying
	//       array [] operators to it!  It's counted, because if IGLUPtr<> is a local variable, it gets
	//       deleted when it leaves scope, which would delete m_ptr, which would delete T*.  And that's
	//       just not what we want.  We want the pointer to be persistant until all references are gone,
	//       at which point m_ptr will garbage collect itself and delete the internal T* object.
	IGLUCountedPtr<T> m_ptr;

	// NOTE: This cheater is to reduce the template/virtual method overhead when the developer
	//       accesses the wrapped, counted pointer with ->, *, or [] operators.  DO NOT 
	//       use this member for anything else except when handling developer queries to
	//       the underlying actual templated data structure T*.   For EVERYTHING else,
	//       including "simple" tasks like checking if NULL, query the smart pointer: m_ptr.  
	//       If you ignore this warning without a GREAT deal of care, this whole templated mess
	//       could fall down on your face.  Since IGLU depends on IGLUPtr<> for a couple important
	//       classes, you will get a ton of compilation errors or random segfaults!
	T* m_cheaterPtr;

};


// End namespace iglu
}

#endif
