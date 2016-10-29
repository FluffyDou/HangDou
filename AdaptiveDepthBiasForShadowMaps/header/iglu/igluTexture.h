/******************************************************************/
/* igluTexture.h                                                  */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class that stores a texture.         */
/*     This is very similar to the image class, but also defines  */
/*     access patterns for the texture with interpolation.        */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef IGLU_TEXTURE_H
#define IGLU_TEXTURE_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include <cassert>

namespace iglu {

class IGLUShaderVariable;

class IGLUTexture
{
	friend class IGLUShaderVariable;
public:
    IGLUTexture();
    virtual ~IGLUTexture();

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void )                   { m_initialized = true; }

	// Returns a boolean that notes if Initialize() has been called (successfully)
	inline bool IsValid() const                       { return m_initialized; }

	// Does this texture need updates?  
	virtual void Update( void )                       {}
	virtual void Update( float frameTime )            {}
	virtual bool NeedsUpdates( void ) const	          { return false; }

	// Bind or unbind this texture.  To enable these methods to be const, which
	//    has some benefits internally in IGLU, the correct texture unit must be
	//    passed in upon binding *and* unbinding. 
	void Bind( GLenum textureUnit=GL_TEXTURE0 ) const;
	void Unbind( GLenum textureUnit=GL_TEXTURE0 ) const;

	// Bind or unbind this texture to an image unit.
	void BindToImageUnit(int imageUnit = 0) const;
	void UnbindFromImageUnit(int imageUnit = 0) const;
	
	//If we bind this texture to an image unit, how will we access it?
	void SetImageAccess(unsigned int flags); 
	
	//Returns the GLenum for the access type (e.g., GL_READ_ONLY for IGLU_IMAGE_READ)
	GLenum GetImageAccess( void ) const {return m_imageAccess;}

	// Returns the GL handle for this texture
	inline GLuint GetTextureID() const                { assert(m_initialized); return m_texID; }

	// Returns the OpenGL type (e.g., GL_TEXTURE_2D) of this particular texture (needs to be defined!)
	virtual GLenum GetTextureType( void ) const         = 0;

	// Returns the OpenGL format of this particular rendering texture! (e.g., GL_RGBA32F)
	virtual GLenum GetTextureFormat( void ) const       = 0;

	// Query the texture to find out some information about its type.  These vary with the texture
	//    type (e.g., GetTextureType() can tell you these), but these provide specific information
	//    in a slightly more user-friendly manner
	virtual bool Is1DTexture( void ) const              = 0;
	virtual bool Is2DTexture( void ) const              = 0;
	virtual bool Is3DTexture( void ) const              = 0;
	virtual bool IsBufferTexture( void ) const          = 0;
	virtual bool IsRectTexture( void ) const            = 0;
	virtual bool IsCubeTexture( void ) const            = 0;
	virtual bool IsMultisampleTexture( void ) const     = 0;
	virtual bool IsArrayTexture( void ) const           = 0;

	// Query to see if the texture is generating mipmaps (to allow mipmap access)
	virtual bool IsTextureMipmappable( void ) const   { return false; }

	// Get size of the texture
	inline int GetWidth()  const                      { return m_width;  }               
	inline int GetHeight() const                      { return m_height; }              
	inline int GetDepth()  const                      { return m_depth;  }            

	// Return the filename (if any) of the texture
	inline const char *GetFilename( void ) const      { return m_filename; }

	// Set and/or get an internal, arbitrary process-specific texture name. 
	//      This is entirely optional for the use of this class, but I find
	//      naming my resources can often be very, very useful.
	inline void SetTextureName( char *newName )       { if (m_texName) free(m_texName); m_texName = strdup( newName ); }
	inline const char *GetTextureName( void ) const   { return m_texName; }

	// Set methods to change various OpenGL texture modes.  
	//    This can change the wrap mode and/or filtering mode
	virtual void SetTextureParameters( unsigned int flags )              {}
	virtual void SetTextureCompareMode( unsigned int flags, bool isTexBound=true ) const;

	// A pointer to a IGLUTexture could have type IGLUTexture::Ptr
	typedef IGLUTexture *Ptr;

protected:
	char *m_filename, *m_texName;
    int m_width, m_height, m_depth;
	GLuint m_texID;
	bool m_initialized, m_compressTex;

	//Stores how we want to access this texture when it is bound to image units
	GLenum m_imageAccess;

	// Utility functions to convert from the IGLU flags to OpenGL params
	GLint  GetMinification( unsigned int flags ) const;
	GLint  GetMagnification( unsigned int flags ) const;
	GLint  GetSWrap( unsigned int flags ) const;
	GLint  GetTWrap( unsigned int flags ) const;

	// Utility functions to answer questions based on input flags
	bool   UsingMipmaps( unsigned int flags ) const;

	// A utility function to convert our image types (typically GL_RGB and 
	//     GL_RGBA) into a corresponding compressed format.
	GLenum GetCompressedFormat( GLenum imgType ) const;
	bool   CanCompressFormat( GLenum imgType ) const;

	// Get datatype (given GL_RGBA8 or GL_R32UI get the data type, e.g., GL_FLOAT)
	GLenum GetFormatDataType( GLenum imgFormat ) const;

	// Check if the image format specified is renderable (e.g., for a renderbuffer)
	bool IsFormatRenderable( GLenum imgFormat ) const;

	// Take a sized internal format and return a simple internal format without size
	//    (One of: GL_RED, GL_RG, GL_RGB, GL_RGBA, GL_DEPTH_COMPONENT, GL_DEPTH_STENCIL)
	GLenum GetBaseInternalFormat( GLenum imgFormat ) const;

	// A utility to determine if this format is a depth buffer
	bool IsDepthFormat( GLenum imgFormat ) const;
};


// Some common flags we'll be using for out texture classes.
enum {
	IGLU_FLAGS_NONE            = 0x000000,
	IGLU_MAG_NEAREST           = 0x000001,
	IGLU_MAG_LINEAR            = 0x000002,
	IGLU_REPEAT_S              = 0x000004,
	IGLU_MIRROR_REPEAT_S       = 0x000008,
	IGLU_CLAMP_S               = 0x000010,
	IGLU_CLAMP_TO_EDGE_S       = 0x000020,
	IGLU_CLAMP_TO_BORDER_S     = 0x000040,
	IGLU_REPEAT_T              = 0x000080,
	IGLU_MIRROR_REPEAT_T	   = 0x000100,
	IGLU_CLAMP_T               = 0x000200,
	IGLU_CLAMP_TO_EDGE_T       = 0x000400,
	IGLU_CLAMP_TO_BORDER_T     = 0x000800,
	IGLU_MIN_NEAREST           = 0x001000,
	IGLU_MIN_LINEAR            = 0x002000,
	IGLU_MIN_NEAR_MIP_NEAR     = 0x004000,
	IGLU_MIN_NEAR_MIP_LINEAR   = 0x008000,
	IGLU_MIN_LINEAR_MIP_NEAR   = 0x010000,
	IGLU_MIN_LINEAR_MIP_LINEAR = 0x020000,
	IGLU_REPEAT_VIDEO          = 0x040000,
	IGLU_COMPRESS_TEXTURE      = 0x080000,
	IGLU_TEX_COMPARE_REF       = 0x100000,
	IGLU_TEX_COMPARE_NONE      = 0x200000,
	IGLU_IMAGE_READ 		   = 0x400000,
	IGLU_IMAGE_WRITE           = 0x800000,
	IGLU_IMAGE_READ_WRITE	   = IGLU_IMAGE_READ | IGLU_IMAGE_WRITE,
	IGLU_TEXTURE_DEFAULT       = IGLU_COMPRESS_TEXTURE | IGLU_MIN_LINEAR | IGLU_MAG_LINEAR | IGLU_CLAMP_TO_EDGE_S | IGLU_CLAMP_TO_EDGE_T
};

// End namespace iglu
}


#endif
