/******************************************************************/
/* igluPixelBufferTexture2D.h                                     */
/* -----------------------                                        */
/*                                                                */
/*       EXPERIMENTAL                                             */
/*                                                                */
/* This texture is a 2D texture using a pixel buffer object for   */
/*     internal data storage.  There may be lots of caveats to    */
/*     using this class...  I hope to list them here, eventually. */
/*                                                                */
/* Chris Wyman (12/11/2012)                                       */
/******************************************************************/

#ifndef IGLU_PIXEL_BUFFER_TEXTURE2D_H
#define IGLU_PIXEL_BUFFER_TEXTURE2D_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluTexture.h"
#include "igluImage.h"

namespace iglu {


class IGLUPixelBufferTexture2D : public IGLUTexture
{
public:
	IGLUPixelBufferTexture2D( GLenum format, int width, int height, 
		                      unsigned int flags = IGLU_TEXTURE_DEFAULT, 
							  bool initializeImmediately=true );
    virtual ~IGLUPixelBufferTexture2D();

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void );

	// Returns the OpenGL type of this particular texture (needs to be defined!)
	virtual GLenum GetTextureType( void ) const         { return GL_TEXTURE_2D; }

	// Query detailed information about this texture type
	virtual bool Is1DTexture( void ) const              { return false; }
	virtual bool Is2DTexture( void ) const              { return true; }
	virtual bool Is3DTexture( void ) const              { return false; }
	virtual bool IsBufferTexture( void ) const          { return false; }
	virtual bool IsRectTexture( void ) const            { return false; }
	virtual bool IsCubeTexture( void ) const            { return false; }
	virtual bool IsMultisampleTexture( void ) const     { return false; }
	virtual bool IsArrayTexture( void ) const           { return false; }

	// Does this texture need updates? [ Not sure if updating this is the right way...]  
	virtual void Update( void )                       {}
	virtual void Update( float frameTime )            {}
	virtual bool NeedsUpdates( void ) const	          { return false; }

	// Set/update the texture parameters for this texture
	virtual void SetTextureParameters( unsigned int flags );

	// This reads from a texture image.  Use the IGLUImage-specified pixel format
	virtual GLenum GetTextureFormat( void ) const       { return m_pixelFormat; }

	// A pointer to a IGLUTexture2D could have type IGLUTexture2D::Ptr
	typedef IGLUPixelBufferTexture2D *Ptr;

protected:
	// OpenGL texture settings
	GLint m_minFilter, m_magFilter;
	GLint m_sWrap, m_tWrap;
	GLint m_width, m_height;

	// The pixel format.  This should be a GL format like GL_RGBA8 or GL_R32F
	GLenum m_pixelFormat;

	bool m_mipmapsNeeded;

	void UpdateTextureParameters( void );
};


// End namespace iglu
}


#endif
