/******************************************************************/
/* igluTextureLightprobeCubemap.h                                 */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class that loads a 2D image that     */
/*     contains a lightprobe stored as a cubemap in an unfolded   */
/*     box (i.e., cross) image format.  This loads the file and   */
/*     constructs an OpenGL cubemap texture that can be used.     */
/*                                                                */
/* Chris Wyman (01/30/2012)                                       */
/******************************************************************/

#ifndef IGLU_TEXTURELIGHTPROBECUBEMAP_H
#define IGLU_TEXTURELIGHTPROBECUBEMAP_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluTexture.h"
#include "igluImage.h"

namespace iglu {


class IGLUTextureLightprobeCubemap : public IGLUTexture
{
public:
    IGLUTextureLightprobeCubemap( char *filename, unsigned int flags = IGLU_TEXTURE_DEFAULT, bool initializeImmediately=true );
    virtual ~IGLUTextureLightprobeCubemap();

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void );

	// Returns the OpenGL type of this particular texture (needs to be defined!)
	virtual GLenum GetTextureType( void ) const         { return GL_TEXTURE_CUBE_MAP; }

	// Query detailed information about this texture type
	virtual bool Is1DTexture( void ) const              { return false; }
	virtual bool Is2DTexture( void ) const              { return false; }
	virtual bool Is3DTexture( void ) const              { return false; }
	virtual bool IsBufferTexture( void ) const          { return false; }
	virtual bool IsRectTexture( void ) const            { return false; }
	virtual bool IsCubeTexture( void ) const            { return true; }
	virtual bool IsMultisampleTexture( void ) const     { return false; }
	virtual bool IsArrayTexture( void ) const           { return false; }

	// Query to see if the texture is generating mipmaps (to allow mipmap access)
	virtual bool IsTextureMipmappable( void ) const     { return m_mipmapsNeeded; }

	// Set/update the texture parameters for this texture
	virtual void SetTextureParameters( unsigned int flags );

	// This reads from a texture image.  Types are pretty basic. 
	//   --> If you need more complex types, you might want to overload IGLUTexture yourself!
	virtual GLenum GetTextureFormat( void ) const       { return ( m_pixelFormat == GL_RGBA ? GL_RGBA8 : GL_RGB8 ); }

	// A pointer to a IGLUTexture2D could have type IGLUTexture2D::Ptr
	typedef IGLUTextureLightprobeCubemap *Ptr;

protected:
	// Our image-loader class
	IGLUImage *m_texImg;

	// OpenGL texture settings
	GLint m_minFilter, m_magFilter;

	// The pixel format.  For these images either GL_RGB or GL_RGBA
	GLenum m_pixelFormat;

	bool m_mipmapsNeeded;
};


// End namespace iglu
}


#endif
