/******************************************************************/
/* igluVideoTexture2D.h                                           */
/* -----------------------                                        */
/*                                                                */
/* The file defines a texture class that encapsulates a video     */
/*    texture loaded using the IGLUVideo class.  The IGLUVideo    */
/*    class (and hence the IGLUVideoTexture2D class) require      */
/*    linking with the FFmpeg libraries.  See igluVideo.h for     */
/*    more details.                                               */
/*                                                                */
/* Chris Wyman (10/06/2011)                                       */
/******************************************************************/

#ifndef IGLU_VIDEOTEXTURE2D_H
#define IGLU_VIDEOTEXTURE2D_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluTexture.h"
#include "igluVideo.h"

namespace iglu {


class IGLUVideoTexture2D : public IGLUTexture
{
public:
    IGLUVideoTexture2D( char *filename, unsigned int flags = IGLU_TEXTURE_DEFAULT, bool initializeImmediately=true );
    virtual ~IGLUVideoTexture2D();

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

	// This reads from a video.  Types are pretty basic. 
	//   --> If you need more complex types, you might want to overload IGLUTexture yourself!
	virtual GLenum GetTextureFormat( void ) const       { return ( m_pixelFormat ); } // == GL_RGBA ? GL_RGBA8 : GL_RGB8 ); }

	// Does this texture need updates?  
	virtual void Update( void );                       
	virtual void Update( float frameTime );          
	virtual bool NeedsUpdates( void ) const	          { return true; }

	// A pointer to a IGLUVideoTexture2D could have type IGLUVideoTexture2D::Ptr
	typedef IGLUVideoTexture2D *Ptr;

protected:
	// Our image-loader class
	IGLUVideo *m_videoTex;

	// The pixel format.  For videos either GL_RGB or (possibly) GL_RGBA
	GLenum m_pixelFormat;

	// OpenGL texture settings
	GLint m_minFilter, m_magFilter;
	GLint m_sWrap, m_tWrap;
	bool  m_repeatVideo, m_mipmapsNeeded;
};


// End namespace iglu
}


#endif
