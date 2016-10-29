/******************************************************************/
/* igluTransformFeedback.h                                        */
/* -----------------------                                        */
/*                                                                */
/* A class for creating and using OpenGL's transform feedback     */
/*     mechanism.  To output values to a GPU memory buffer after  */
/*     the vertex/geometry GLSL stages (and prior to rasterizing) */
/*                                                                */
/* Chris Wyman (01/16/2012)                                       */
/******************************************************************/

#ifndef __IGLU__TRANSFORM_FEEDBACK_H_
#define __IGLU__TRANSFORM_FEEDBACK_H_

#include <GL/glew.h>
#include "iglu/igluBuffer.h"

namespace iglu { 

class IGLUTransformFeedback {
public:
	IGLUTransformFeedback( bool discardFragments=false );
	~IGLUTransformFeedback();

	// Routines to bind / unbind the feedback object
	void Bind( void );
	void Unbind( void );

	// Start using the transform feedback
	int Begin( GLenum primitiveType );
	int End( void );

	// While using transform feedback, if one wants to stop/pause without
	//    restarting from the beginning of the buffer, use Pause()/Restart().
	//    Note there are some limitations to most changes to transform feedback
	//    state when paused.
	int Pause( void );
	int Resume( void );

	// Bind buffers to the transform feedback object.  If outputting using GL_INTERLEAVED_ATTRIBS,
	//    it is likely only one buffer (on attachIdx=0) is needed.  If separating attributes into
	//    different buffers, one buffer must be attached for each attribute to be output.  If
	int AttachBuffer( const IGLUBuffer::Ptr &buf, GLuint attachIdx=0 );
	int AttachSubBuffer( const IGLUBuffer::Ptr &buf, GLintptr startOffsetBytes, GLsizeiptr sizeBytes, GLuint attachIdx=0 );

	// Accessors to see if the buffer is mapped to a host pointer
	inline bool IsBound( void ) const					{ return m_isBound; }
	inline bool IsActive( void ) const					{ return m_isActive; }
	inline bool IsPaused( void ) const					{ return m_isPaused; }

	// Accessors to get the OpenGL state handle
	inline GLuint GetTransformID( void ) const          { return m_transformID; }
	
	// A pointer to a IGLUVertexArray should have type IGLUVertexArray::Ptr
	typedef IGLUTransformFeedback *Ptr;

private:

	// The OpenGL ID of the transform feedback object
	GLuint m_transformID;

	// Check if this transform feedback object is currently bound
	bool m_isBound;

	// Check if this transform feedback is currently active/paused
	bool m_isActive, m_isPaused;

	// This class allows two ways to enable/bind the object.  It can be
	//     done explicitly or implicitly when calling Begin().  If called
	//     implicitly with Begin(), we need to disable/unbind when calling End()
	//     rather than waiting for the user to explicitly call Unbind().  This
	//     variable tracks which binding method was used.
	bool m_boundOnBegin;

	// Determines if (when running in xformFeedback mode) we discard rasterized fragments
	bool m_discardFragments;
};


} // End namespace iglu


#endif