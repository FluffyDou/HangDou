/******************************************************************/
/* igluRenderTexture.h                                            */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class that stores a 2D image texture */
/*     in a format for easy display for OpenGL.                   */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef IGLU_RENDERTEXTURE_H
#define IGLU_RENDERTEXTURE_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluTexture.h"
#include "igluImage.h"

namespace iglu {


class IGLURenderTexture : public IGLUTexture
{
public:
	IGLURenderTexture( GLenum format = GL_RGBA8 );
	virtual ~IGLURenderTexture()                        {}

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void )                     { m_initialized = true; }
	
	// Returns the OpenGL format of this particular rendering texture! (e.g., GL_RGBA32F)
	virtual GLenum GetTextureFormat( void ) const       { return m_format; }

	// Returns if this renderable texture needs to be mipmapped
	virtual bool IsTextureMipmappable( void ) const     { return m_mipmapsNeeded; } 

	// Generate the mipmaps, if available
	virtual void GenerateMipmaps( bool leaveBound )     = 0;

	// A pointer to a IGLURenderTexture could have type IGLURenderTexture::Ptr
	typedef IGLURenderTexture *Ptr;

protected:
	// What is the format of each pixel/texel/voxel?  GL_RGBA8?  GL_R32UI?  
	GLenum m_format;

	// What is the type of data for each texels?  GL_FLOAT?  GL_UNSIGNED_BYTE?  GL_INT?
	GLenum m_dataType;

	// Should we use mipmaps when we render?
	bool m_mipmapsNeeded;
};


// End namespace iglu
}


#endif
