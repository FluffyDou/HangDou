/******************************************************************/
/* igluRenderTexture2dMultisample.h                               */
/* --------------------------------                               */
/*                                                                */
/* The file defines a class derived from IGLURenderTexture that   */
/*     is a renderable 2D multisample texture, which can also be  */
/*     used as an input for a GLSL multisample sampler.           */
/*                                                                */
/* Chris Wyman (01/06/2012)                                       */
/******************************************************************/

#ifndef IGLU_RENDERTEXTURE2DMULTISAMPLE_H
#define IGLU_RENDERTEXTURE2DMULTISAMPLE_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluRenderTexture.h"

namespace iglu {


class IGLURenderTexture2DMultisample : public IGLURenderTexture
{
public:
    IGLURenderTexture2DMultisample( int width, int height, int numSamples, GLenum format = GL_RGBA8, 
		                            bool fixedSamples=true, unsigned int flags = IGLU_TEXTURE_DEFAULT, 
									bool initializeImmediately=true );
    virtual ~IGLURenderTexture2DMultisample();

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void );

	// Returns the OpenGL type of this particular texture (needs to be defined!)
	virtual GLenum GetTextureType( void ) const         { return GL_TEXTURE_2D_MULTISAMPLE; }

	// Query detailed information about this texture type
	virtual bool Is1DTexture( void ) const              { return false; }
	virtual bool Is2DTexture( void ) const              { return true; }
	virtual bool Is3DTexture( void ) const              { return false; }
	virtual bool IsBufferTexture( void ) const          { return false; }
	virtual bool IsRectTexture( void ) const            { return false; }
	virtual bool IsCubeTexture( void ) const            { return false; }
	virtual bool IsMultisampleTexture( void ) const     { return true; }
	virtual bool IsArrayTexture( void ) const           { return false; }

	// There are no mipmaps to generate with a multisample buffer.  Get over it.
	virtual void GenerateMipmaps( bool leaveBound )     {}

	// A pointer to a IGLURenderTexture2DMultisample could have type IGLURenderTexture2DMultisample::Ptr
	typedef IGLURenderTexture2DMultisample *Ptr;

protected:
	// Multisampling settings
	GLsizei m_numSamples;
	GLboolean m_fixedSampLocations;

	// OpenGL texture settings
	GLint m_minFilter, m_magFilter;
	GLint m_sWrap, m_tWrap;
};


// End namespace iglu
}


#endif
