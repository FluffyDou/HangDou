/******************************************************************/
/* igluRandomTexture3D.h                                          */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class for passing a 3D array of      */
/*    random unsigned intergers down to a GLSL program.  This     */
/*    is based of the internal IGLU random number generator to    */
/*    create random seeds.                                        */
/*                                                                */
/* Note: When using the random seeds on the GPU, ensure you use   */
/*    a random number generator that uses *different* adder and   */
/*    multipliers than the built-in IGLU generator, otherwise you */
/*    may get strange correlation artifacts between adjacent      */
/*    pixels!                                                     */
/*                                                                */
/* Chris Wyman (12/07/2011)                                       */
/******************************************************************/

#ifndef IGLU_RANDOMTEXTURE3D_H
#define IGLU_RANDOMTEXTURE3D_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluTexture.h"

namespace iglu {

class IGLURandom;

class IGLURandomTexture3D : public IGLUTexture
{
public:
    IGLURandomTexture3D( uint width, uint height, uint depth,                       // Texture width, height, depth
		                 GLenum       textureType = GL_INT,                         // Either GL_INT or GL_FLOAT texture
		                 unsigned int randomSeed = 0,                               // Seed for random generator
		                 unsigned int flags = IGLU_MIN_NEAREST | IGLU_MAG_NEAREST,  // Use nearest interp by default
						 bool initializeImmediately=true );                         // If no GL context yet, should be false
    virtual ~IGLURandomTexture3D();

	// Initialize() must be called after an OpenGL context has been created.
	virtual void Initialize( void );

	// Returns the OpenGL type of this particular texture (needs to be defined!)
	virtual GLenum GetTextureType( void ) const         { return GL_TEXTURE_3D; }

	// Query detailed information about this texture type
	virtual bool Is1DTexture( void ) const              { return false; }
	virtual bool Is2DTexture( void ) const              { return false; }
	virtual bool Is3DTexture( void ) const              { return true; }
	virtual bool IsBufferTexture( void ) const          { return false; }
	virtual bool IsRectTexture( void ) const            { return false; }
	virtual bool IsCubeTexture( void ) const            { return false; }
	virtual bool IsMultisampleTexture( void ) const     { return false; }
	virtual bool IsArrayTexture( void ) const           { return false; }

	// Each random pixel has either 1 random float or 1 random int, depending on constructor parameters
	virtual GLenum GetTextureFormat( void ) const       { return ( m_randType == GL_FLOAT ? GL_R32F : GL_R32I ); }

	// A pointer to a IGLUTexture2D could have type IGLUTexture2D::Ptr
	typedef IGLURandomTexture3D *Ptr;

	//A method for resizing the texture
	void Resize( int width, int height, int depth );

protected:
	// OpenGL texture settings
	GLint m_minFilter, m_magFilter;
	GLint m_sWrap, m_tWrap, m_uWrap;

	// Random seed
	GLuint m_rndSeed;
	IGLURandom *m_rng;

	// What type of random texture is this?  GL_INT or GL_FLOAT?
	GLenum m_randType;

	// Not sure if mipmapping works (not tested)...  But let's leave it in.  It might be desirable.
	bool m_mipmapsNeeded;
};


// End namespace iglu
}


#endif
