/******************************************************************/
/* igluRenderTexture2DArray.h                                     */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class that stores a 2D array texture */
/*     in a renderable format for easy OpenGL usage.              */
/*                                                                */
/* Chris Wyman (01/30/2012)                                       */
/******************************************************************/

#ifndef IGLU_RENDERTEXTURE2DARRAY_H
#define IGLU_RENDERTEXTURE2DARRAY_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluRenderTexture.h"

namespace iglu {


class IGLURenderTexture2DArray : public IGLURenderTexture
{
public:
    IGLURenderTexture2DArray( int width, int height, int layers, GLenum format = GL_RGBA8, unsigned int flags = IGLU_TEXTURE_DEFAULT, bool initializeImmediately=true );
    virtual ~IGLURenderTexture2DArray();

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void );

	// Returns the OpenGL type of this particular texture (needs to be defined!)
	virtual GLenum GetTextureType( void ) const         { return GL_TEXTURE_2D_ARRAY; }

	// Query detailed information about this texture type
	virtual bool Is1DTexture( void ) const              { return false; }
	virtual bool Is2DTexture( void ) const              { return true; }
	virtual bool Is3DTexture( void ) const              { return false; }
	virtual bool IsBufferTexture( void ) const          { return false; }
	virtual bool IsRectTexture( void ) const            { return false; }
	virtual bool IsCubeTexture( void ) const            { return false; }
	virtual bool IsMultisampleTexture( void ) const     { return false; }
	virtual bool IsArrayTexture( void ) const           { return true; }

	// Generate the mipmaps, if available
	virtual void GenerateMipmaps( bool leaveBound=false ); 

	// Set/update the texture parameters for this texture
	virtual void SetTextureParameters( unsigned int flags );

	// A pointer to a IGLURenderTexture2D could have type IGLURenderTexture2D::Ptr
	typedef IGLURenderTexture2DArray *Ptr;

protected:
	// OpenGL texture settings
	GLint m_minFilter, m_magFilter;
	GLint m_sWrap, m_tWrap;
};


// End namespace iglu
}


#endif
