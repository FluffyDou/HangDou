/******************************************************************/
/* igluTextureBuffer.h                                            */
/* -----------------------                                        */
/*                                                                */
/* The file defines an OpenGL texture buffer class that allows    */
/*    the creation of a texture buffer object that can be used    */
/*    with shader programs and supports reading of arbitrary      */
/*    IGLUBuffers inside GLSL programs using a samplerBuffer.     */
/*                                                                */
/* Chris Wyman (01/06/2012)                                       */
/******************************************************************/

#ifndef IGLU_TEXTUREBUFFER_H
#define IGLU_TEXTUREBUFFER_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluTexture.h"
#include "igluImage.h"

namespace iglu {

class IGLUBuffer;


class IGLUTextureBuffer : public IGLUTexture
{
public:
    IGLUTextureBuffer( bool initializeImmediately=true );
    virtual ~IGLUTextureBuffer();

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void );

	// Bind an IGLUBuffer to this texture buffer.  
	//     Note: IGLU does some run-time type checking when assigning textures
	//     to samplers in a GLSL shader.  One can break this type checking if
	//     assigning/changing the texture's bound buffer after the texture has
	//     already been assigned to a GLSL sampler.  By default, this method
	//     will warn when the user changes types in this way.  This behavior
	//     is usually desirable unless you really know what you're doing!
	void BindBuffer( GLenum bufTexelType, const IGLUBuffer *buffer, bool warnOnTypeChange=true );

	// Unbind current IGLUBuffer from the texture.  Usually this is UNNECESSARY!
	void UnbindBuffer( void );

	// Returns the OpenGL type of this particular texture (needs to be defined!)
	virtual GLenum GetTextureType( void ) const         { return GL_TEXTURE_BUFFER; }

	// Query detailed information about this texture type
	virtual bool Is1DTexture( void ) const              { return false; }
	virtual bool Is2DTexture( void ) const              { return false; }
	virtual bool Is3DTexture( void ) const              { return false; }
	virtual bool IsBufferTexture( void ) const          { return true; }
	virtual bool IsRectTexture( void ) const            { return false; }
	virtual bool IsCubeTexture( void ) const            { return false; }
	virtual bool IsMultisampleTexture( void ) const     { return false; }
	virtual bool IsArrayTexture( void ) const           { return false; }

	// The internal type of the currently bound buffer
	virtual GLenum GetTextureFormat( void ) const       { return m_internalFormat; }

	// A pointer to a IGLUTexture2D could have type IGLUTexture2D::Ptr
	typedef IGLUTextureBuffer *Ptr;

protected:
	// Our currently bound buffer
	const IGLUBuffer *m_texBuffer;

	// The internal texel format (e.g., GL_RGBA8 or GL_RG32UI)
	GLenum m_internalFormat;
};


// End namespace iglu
}


#endif
