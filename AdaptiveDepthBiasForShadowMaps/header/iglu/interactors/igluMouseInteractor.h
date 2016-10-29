/************************************************************************/
/* igluMouseInteractor.h                                                */
/* ---------------------                                                */
/*                                                                      */
/* This is a base class for all mouse-interactors (i.e., things that    */
/*     update variables based on mouse clicks, movement, and releases)  */
/*                                                                      */
/* The key methods are:                                                 */
/*     SetOnClick( int x, int y )     // When user clicks mouse         */
/*     UpdateOnMotion( int x, int y ) // Updates during mouse motion    */
/*     Release( void )                // When user releases mouse       */
/*                                                                      */
/* Some interactors may need to keep track of the window size for       */
/* currect interaction, which is handled by:                            */
/*     ResizeInteractorWindow( int width, int height )                  */
/* This method has a default to sent internal m_width & m_height values */
/* but may be redefined if more complex tracking is required.           */
/*                                                                      */
/* Only advanced interactors may need to change methods tracking if     */
/* changes have occured (HasChanged()), executing on re-initialization  */
/* of the interactor (Reset()), and settors that directly up the        */
/* internal matrix (SetMatrix()).                                       */
/*                                                                      */
/*                                                                      */
/* Chris Wyman (1/24/2012)                                              */
/************************************************************************/


#ifndef IGLU_MOUSE_INTERACTOR_H
#define IGLU_MOUSE_INTERACTOR_H

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>

#include "iglu/vectors/vec3.h"

namespace iglu {

class vec4;
class IGLUMatrix4x4;

class IGLUMouseInteractor
{
public:

	// sets up an interactor with window size of width x height 
	IGLUMouseInteractor( int width, int height );
	virtual ~IGLUMouseInteractor()                                         {}

	/************************************************************************
	**  These three methods must be defined for every interactor!          **
	************************************************************************/

	// functions for updating the interactor due to mouse input.           
	//    x & y are the window coordinates of the current mouse position.
	//    The return value of UpdateOnMotion() is usually ignored, but should
	//    return something related how far we've moved (could be in pixels or
	//    degrees, for instance).
	virtual void SetOnClick( int x, int y )                              = 0;
	virtual float UpdateOnMotion( int x, int y )                         = 0;
	virtual void Release( void )                                         = 0;

	/************************************************************************
	**  These X methods MAY be redefined for every interactor!          **
	**     Be careful, however, as the base class methods perform useful   **
	**     operations that may need to be replicated in the derived class. **
	************************************************************************/

	// Computations may rely on knowing correct window width & height      
	virtual void ResizeInteractorWindow( int width, int height );         

	// checks if the interactor has been updated since the most recent of: 
	//   a) interactor initialization or b) the last call to HasChanged()  
	//   Please note: HasChanged() destructively checks the internal value 
	//   (i.e., it sets the internal hasChanged variable to false!)        
	virtual bool HasChanged( void );

	// Reset the interactor to its initial state... 
	virtual void Reset( void );


	/************************************************************************
	**  The following methods should only need redefinition, if a derived  **
	**     class needs to know about all changes to the internal matrices. **
	************************************************************************/
	virtual void SetMatrix( float *newMat );
	virtual void SetMatrix( const IGLUMatrix4x4 &newMat );

	/************************************************************************
	**  The following methods should not need redefinition, as they only   **
	**     get matrix values stored inside this base class.                **
	************************************************************************/
	
	// Get the interactor matrix 
	void GetMatrix( float *matrix );
	inline const IGLUMatrix4x4 &GetMatrix( void ) const          { return m_matrix; }

	// Get the interactor's inverse 
	void GetInverseMatrix( GLfloat *matrix );
	inline const IGLUMatrix4x4 &GetInverseMatrix( void ) const   { return m_inverseMatrix; }

	// Apply the current interactor matrix to a vec4 
	void ApplyMatrix( float inVec[4], float result[4] );
	vec4 ApplyMatrix( const vec4 &vec );


	// A pointer to an interactor should be IGLUMouseInteractor::Ptr
	typedef IGLUMouseInteractor *Ptr;

protected:
	int m_hasChanged;
	int m_currentlyTracking;
	int m_width, m_height;
	IGLUMatrix4x4 m_matrix;
	IGLUMatrix4x4 m_inverseMatrix;
};


// End namespace iglu
}

#endif
