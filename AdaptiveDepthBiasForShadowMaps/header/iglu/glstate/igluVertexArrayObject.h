/******************************************************************/
/* igluVertexArrayObject.h                                        */
/* -----------------------                                        */
/*                                                                */
/* A class for creating and using OpenGL Vertex Array Objects.    */
/*                                                                */
/* Object declaration looks like this:                            */
/*    IGLUVertexArray::Ptr vertArr = 0;                           */
/*                                                                */
/* Basic object initialization looks like this:                   */
/*    float geomData[] = { 0, 0, 0,  1, 0, 0,  0, 1, 0 };         */
/*    vertArr = new IGLUVertexArray();                            */
/*    vertArr->SetVertexArray( sizeof( geomData ), geomData );    */
/*    vertArr->EnableAttribute( IGLU_ATTRIB_VERTEX,               */
/*                            3, GL_FLOAT, 0, BUFFER_OFFSET(0) ); */
/* (Here EnableAttribute() calls glVertexAttribPointer() with the */
/* attribute type [3 floats per vertex], the data stride [0], and */
/* the pointer to the first data [BUFFER_OFFSET(0)] to be passed  */
/* into the shader attribute for vertices (specified in GLSL with */
/* a layout(location = IGLU_VERTEX); IGLU_VERTEX is defined as 0) */
/*                                                                */
/* Basic rendering with a vertex array object looks like:         */
/*    vertArr->DrawArrays( GL_TRIANGLES, 0, 3 );                  */
/*                                                                */
/* Chris Wyman (06/13/2011)                                       */
/******************************************************************/

#ifndef __IGLU__VERTEX_ARRAY_OBJECT__
#define __IGLU__VERTEX_ARRAY_OBJECT__

#include <GL/glew.h>
#include "iglu/igluBuffer.h"
#include "iglu/igluShaderVariable.h"
#include "iglu/glstate/igluTransformFeedback.h"

namespace iglu { 

class IGLUVertexArray {
public:
	IGLUVertexArray();
	virtual ~IGLUVertexArray();

	// Routines to bind / unbind the array
	virtual void Bind( void );
	virtual void Unbind( void );

	// Fill in vertex or element buffer data from an array from system memory.  
	//    A NULL input array will size the buffer(s), but they will be filled 
	//    with garbage (undefined). 
	// An element array need not be created unless using element-based drawing
	//    (i.e., array-based draw calls do not need an element array)
	virtual void SetVertexArray       ( GLsizeiptr bufSize,   void *bufData=0,    int use = IGLU_DYNAMIC | IGLU_DRAW );
	virtual void SetVertexArraySubset ( GLintptr   bufOffset, GLsizeiptr subSize, void *subData );

	virtual void SetElementArray      ( GLenum elemType, GLsizeiptr bufSize,   void *bufData=0,    int use = IGLU_DYNAMIC | IGLU_DRAW );
	virtual void SetElementArraySubset( GLintptr   bufOffset, GLsizeiptr subSize, void *subData );

	// Routines to map / unmap the arrays to CPU-side pointers
	virtual void *MapVertexArray  ( AccessMode mode = IGLU_READ );
	virtual void *MapElementArray ( AccessMode mode = IGLU_READ );
	virtual void UnmapVertexArray ( void );
	virtual void UnmapElementArray( void );

	// Methods for setting up vertex array prior to usage
	void EnableAttribute( const IGLUShaderVariable &attrib,         // A shader attribute to bind this to (e.g., shader["vertex"])
		                  GLint size, GLenum type,                  // 'size' values of the specified type are sent down per vertex
						  GLsizei stride=0,                         // Vertices are stride bytes apart.  Default (=0) means densely packed data 
						  const GLvoid *pointer=BUFFER_OFFSET(0),   // Pointer to start of array data.  Default is the start of the buffer (0 bytes)
						  int instancesDivisor=-1 );                // If instancing (and > 0), this is passed to glVertexAttribDivisor

	void EnableAttribute( int glAttribIdx,                          // Typically an built-in IGLU define like IGLU_ATTRIB_VERTEX (matching a GLSL layout specifier)
		                  GLint size, GLenum type,                  // 'size' values of the specified type are sent down per vertex
						  GLsizei stride=0,                         // Vertices are stride bytes apart.  Default (=0) means densely packed data 
						  const GLvoid *pointer=BUFFER_OFFSET(0),   // Pointer to start of array data.  Default is the start of the buffer (0 bytes)
						  int instancesDivisor=-1 );                // If instancing (and > 0), this is passed to glVertexAttribDivisor


	// If you're using primitive restart (e.g., to specify a stream of triangle strips) give the restart index here.
	void EnablePrimitiveRestart( unsigned int restartIndex );

	// Drawing calls.  These all work as described by the gl* equivalent on the underlying data... 
	//    Note that some substantially similar OpenGL draw calls have been combined into one.
	//    (e.g., the spec says glDrawArraysInstanced is the same as glDrawArraysInstancedBaseInstance 
	//           with baseInstance set to 0.  Both use the method DrawArraysInstanced() below.)
	void DrawArrays( GLenum mode, GLint first, GLsizei count );
	void MultiDrawArrays( GLenum mode, GLint *first, GLsizei *count, GLsizei primCount );
	void DrawArraysInstanced( GLenum mode, GLint first, GLsizei count, GLsizei primCount, GLuint baseInstance=0 );
	void DrawElements( GLenum mode, GLsizei count, GLuint bufOffsetBytes=0 );

	// Transform feedback draw routine.  This should only be used if the internal buffers were
	//    populated using OpenGL's transform feedback mechanism.  They draw the number of vertices
	//    that were output in the last transform feedback stage.  ( e.g., if transform feedback  
	//    ID #15 output 45 vertices then DrawTransformFeedback( GL_TRIANGLES, 15 ) is equivalent  
	//    to DrawArrays( GL_TRIANGLES, 0, 45). )
	// Note:  This function behaves as any of the glDrawTransformFeedback*() functions, depending
	//    on if the optional parameters are modified.  To use the equivalent of glTransformFeedback,
	//    leave the default parameters as is.  
	void DrawTransformFeedback( GLenum mode, const IGLUTransformFeedback::Ptr &feedback, GLuint feedbackStream=0, GLsizei instances=1 );

	// Accessors to see if the buffer is mapped to a host pointer
	inline bool IsMapped( void ) const					{ return m_vertArray->IsMapped() || (m_elemArray && m_elemArray->IsMapped()); }

	// Accessors to get the OpenGL state handle
	inline GLuint GetArrayID( void ) const              { return m_arrayID; }
	inline GLuint GetVertexBufferID( void ) const       { return m_vertArray->GetBufferID(); }
	inline GLuint GetElementBufferID( void ) const      { return m_elemArray ? m_elemArray->GetBufferID() : 0; } 

	// Grab the underlying buffers.  Note the element array may be NULL!
	inline const IGLUBuffer::Ptr &GetVertexBuffer( void ) const    { return m_vertArray; }
	inline const IGLUBuffer::Ptr &GetElementBuffer( void ) const   { return m_elemArray; }

	// A pointer to a IGLUVertexArray should have type IGLUVertexArray::Ptr
	typedef IGLUVertexArray *Ptr;

private:

	// The OpenGL ID of the vertex array object
	GLuint m_arrayID;

	// Check if this is array is currently bound
	bool m_isBound;

	// Rendering with primitive restart?
	bool m_withPrimitiveRestart;
	GLuint m_restartIndex;

	// The IGLUBuffers storing the (interleaved) vertex array and the (optional) element array
	IGLUBuffer::Ptr m_vertArray;
	IGLUBuffer::Ptr m_elemArray;

	// When using the data in the element array, we often need to know what type it is...  Store that
	GLenum m_elemType;

	// There's some tricky state issues depending on the order of the Set*Array() and 
	//    EnableAttribute() calls, given that IGLUBuffer doesn't behave quite the way
	//    we want for this class.  So as soon as we start calling EnableAttribute() we
	//    may need to explicitly leave the arrays bound at the end of the Set*Array()....
	//    Currently this is not done, but the variable here would allow it.
	bool m_enableAttribCalled;

	// We only need to bind out element array to the VAO once.  However, we do this lazily,
	//    since the user may add an element array at any time, so keep track of this.
	bool m_elementArrayBound;

	// You cannot call SetVertexArraySubset() before callsing SetVertexArray().  Store state
	//    to detect this;  ideally this will throw a nice error.  Now it doesn't.
	bool m_vertArrayInit;

	// Internal functions used at beginning/end of draw routines for enabling
	//    primitive restart (if needed)
	bool Internal_InitPrimRestart( void );
	void Internal_StopPrimRestart( bool wasEnabled );

	//A query for 32-bit int GL types
	bool IsIntType( GLenum type );

};


} // End namespace iglu


#endif