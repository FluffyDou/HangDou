/******************************************************************/
/* igluBuffer.h                                                   */
/* -----------------------                                        */
/*                                                                */
/* A class that creates and stores and OpenGL buffer object       */
/*                                                                */
/* Chris Wyman (06/13/2011)                                       */
/******************************************************************/

#ifndef __IGLU__BUFFER__
#define __IGLU__BUFFER__

#include <GL/glew.h>

namespace iglu { 

// Type enum declarations
enum BufferType { IGLU_ARRAY, IGLU_ELEMENT_ARRAY, IGLU_PIXEL_PACK, IGLU_PIXEL_UNPACK, IGLU_TEXTURE, IGLU_UNIFORM };
enum AccessMode { IGLU_READ = 0x01, IGLU_WRITE, IGLU_READ_WRITE };
enum { IGLU_DRAW = 0x02, IGLU_COPY = 0x04, IGLU_STREAM = 0x10, IGLU_STATIC = 0x20, IGLU_DYNAMIC = 0x40 };


class IGLUBuffer {

public:
	// Default constructor.  Creates an OpenGL buffer
	IGLUBuffer( BufferType type = IGLU_ARRAY );

	// Default destructor.
	virtual ~IGLUBuffer();

	// Routines to bind / unbind the buffer
	virtual void Bind( void );
	virtual void Unbind( void );

	// Fill buffer data from an array of system memory.  A NULL input array
	//    will size the buffer, but it will be filled with garbage (undefined)
	virtual void SetBufferData   ( GLsizeiptr bufSize,   void *bufData=0,    int use = IGLU_DYNAMIC | IGLU_DRAW );
	virtual void SetBufferSubData( GLintptr   bufOffset, GLsizeiptr subSize, void *subData );

	// Routines to map / unmap the buffer
	virtual void *Map( AccessMode mode = IGLU_READ );
	virtual void Unmap( void );

	// Accessors to see if the buffer is mapped to a host pointer
	inline bool IsMapped( void ) const					{ return m_mapped; }

	// Accessors to get the OpenGL buffer handle
	inline GLuint GetBufferID( void ) const             { return m_bufID; }

	// Accessor to buffer info
	inline GLsizeiptr GetBufferSize( void ) const       { return m_bufSize; }

	// A pointer to a IGLUBuffer could have type IGLUBuffer::Ptr
	typedef IGLUBuffer *Ptr;

protected:
	// Is the class currently initialized?  Mapped to a host pointer?
	bool m_mapped;

	// Is this buffer currently bound to OpenGL?
	bool m_isBound;

	// The OpenGL ID for this buffer.
	GLuint m_bufID;

	// The default type of buffer (though this is generally flexible)
	GLenum m_type;

	// A pointer to the mapped buffer (or NULL if unmapped)
	void *m_bufPtr;

	// Information about this buffer
	GLsizeiptr m_bufSize;

protected:
	// Conversions from internal class enums to OpenGL enums
	GLenum ConvertToGLAccessMode( AccessMode mode );
	GLenum ConvertToGLType      ( BufferType type );
	GLenum ConvertToUsageMode   ( int use );
};

// Since we're creating buffer objects, we'll likely need buffer offsets.  Use this macro.
#ifndef BUFFER_OFFSET
#define BUFFER_OFFSET(i) ((char *)NULL + (i))
#endif

}

#endif

