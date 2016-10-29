/******************************************************************/
/* igluUniformBuffer.h                                            */
/* -----------------------                                        */
/*                                                                */
/* A class that creates an OpenGL uniform buffer and provides     */
/*    some methods for accessing data inside the uniform buffer.  */
/*                                                                */
/* Chris Wyman (06/13/2011)                                       */
/******************************************************************/

#ifndef __IGLU__UNIFORM_BUFFER__
#define __IGLU__UNIFORM_BUFFER__

#include <GL/glew.h>

#include "../igluBuffer.h"
#include "../igluShaderProgram.h"

namespace iglu { 

class IGLUUniformBuffer : public IGLUBuffer {

public:
	// Create a uniform buffer accessed in shaders with the specified name.  To figure out
	//     what variables are stored inside the buffer and do setup, a valid GLSL program
	//     that uses this uniform buffer must be specified.  This need not be the program
	//     used during rendering!
	// If a program ID of 0 is passed to the constructor, the constructor assumes the currently
	//     bound/enabled program should be used to construct the uniform buffer
	IGLUUniformBuffer( const char *bufName, GLuint glslProgramID=0 );
	IGLUUniformBuffer( const char *bufName, IGLUShaderProgram::Ptr glslProgram );

	// Destructor.
	virtual ~IGLUUniformBuffer();

	// Routines to bind the uniform buffer.  It needs a GLSL program to bind to, unlike other
	//     buffers, as the binding location depends on the linked location of the buffer in
	//     each shader.  
	// If no GLSL shader is passed to Bind, it is assumed that the currently bound/enabled
	//     program should be used.
	virtual void Bind( void );
	virtual void Bind( GLuint glslProgramID );
	virtual void Bind( IGLUShaderProgram::Ptr glslProgram );

	// To unbind the buffer, the inherited code from IGLUBuffer should work.
	//virtual void Unbind( void );

	// Fill buffer data from an array of system memory.  A NULL input array
	//    will size the buffer, but it will be filled with garbage (undefined)
	virtual void SetBufferData   ( GLsizeiptr bufSize,   void *bufData=0,    int use = IGLU_DYNAMIC | IGLU_DRAW );
	virtual void SetBufferSubData( GLintptr   bufOffset, GLsizeiptr subSize, void *subData );

	// A pointer to a IGLUBuffer could have type IGLUBuffer::Ptr
	typedef IGLUUniformBuffer *Ptr;

protected:
	GLuint m_bufProgID;  // All data stored internally is based on some GLSL program.  (This is the ID)
	GLint  m_blkIdx;     // A block index selected based upon location of the uniform buffer in m_bufProgID;
};


}

#endif

