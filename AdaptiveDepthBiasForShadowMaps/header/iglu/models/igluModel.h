/******************************************************************/
/* igluModel.h                                                    */
/* -----------------------                                        */
/*                                                                */
/* Defines a base class that all renderable geometry inherits.    */
/*     This provides a base functionality developers can expect   */
/*     from draw-able geometry.                                   */
/*                                                                */
/* Chris Wyman (10/24/2011)                                       */
/******************************************************************/

#ifndef IGLU_MODEL_H
#define IGLU_MODEL_H

#include "iglu/igluShaderProgram.h"
#include "iglu/igluMatrix4x4.h"

namespace iglu {

class IGLUModel 
{
public:
	// Trivial constructor/destructor
	IGLUModel() : m_transform( IGLUMatrix4x4::Identity() ), m_userControlID(0) {}
	virtual ~IGLUModel()                                   {}

	// Draw routines returning either IGLU_NO_ERROR or an error code
	//    -> One enables a specific shader
	//    -> One assumes the developer enabled an appropriate shader
	virtual int Draw( IGLUShaderProgram::Ptr &shader ) = 0;
	virtual int Draw( void ) = 0;

	// Specifies what data is available in the renderable primitive
    virtual bool HasVertices ( void ) const         { return true; }
	virtual bool HasTexCoords( void ) const         { return false; }
	virtual bool HasNormals  ( void ) const         { return false; }
	virtual bool HasMatlID   ( void ) const         { return false; }
    virtual bool HasObjectID ( void ) const		    { return false; }

	// Each model might have an independent model matrix internally
	virtual const IGLUMatrix4x4 &GetModelMatrix( void ) const                { return m_transform; }
	virtual void                 SetModelMatrix( const IGLUMatrix4x4 &mat )  { m_transform = mat; }

	// If there are user controls over this model/geometry, we might need to check that
	//    during rendering and apply appropriate matrices.  That's left up to the developer.
	//    This just allows storing an ID to remember if this geometry is user controllable
	virtual void  SetUserControlID( int ctrlID )    { m_userControlID = ctrlID; }
	virtual int   IsUserControlled( void ) const    { return m_userControlID; }

	// A pointer to a IGLUModel could have type IGLUModel::Ptr
	typedef IGLUModel *Ptr;

protected:
	IGLUMatrix4x4 m_transform;
	int m_userControlID;
};




// End namespace iglu
}



#endif
