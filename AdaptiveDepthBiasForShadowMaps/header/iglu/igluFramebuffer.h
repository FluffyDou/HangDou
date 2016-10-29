/******************************************************************/
/* igluFramebuffer.h                                              */
/* -----------------------                                        */
/*                                                                */
/* The file defines a class that encapsulates framebuffer objects */
/*    for rendering into texture or multiple textures.            */
/*                                                                */
/* Chris Wyman (10/11/2011)                                       */
/******************************************************************/

#ifndef IGLU__FRAMEBUFFER_H
#define IGLU__FRAMEBUFFER_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "igluRenderTexture.h"

namespace iglu {

class IGLUFramebuffer 
{
public:
	IGLUFramebuffer();
    virtual ~IGLUFramebuffer();

	// Bind/unbind this framebuffer.  When bound, this is the buffer being drawn to. 
	virtual void Bind( bool storePrevBinding=false );
	virtual void Unbind( void );

	// Clears the framebuffer.  By default, clears all valid buffers
	//    inside the FBO.  However, passing in a bitmask will only clear
	//    the appropriate buffers (if they exist in this FBO)
	virtual void Clear( uint bitmask=0xFFFFFFFFu );

	// Can use when buffer has MRT.  Uses the first <numBufs> bound buffers.  This is a 
	//   pretty simplistic control over which buffers are drawn to.  More complex versions
	//   may be desirable.  
	// NOTE:  Bind() will automatically call glDrawBufs() to use all bound textures when
	//   m_setDrawBufsOnBind is true.  SetDrawbuffers() is a way to use *fewer* than the
	//   number of buffers bound.
	void SetDrawbuffers( int numBufs );

	// Checks if drawing to this framebuffer works.  Is false if some step during creation fails
	inline bool IsDrawable( void ) const	        { return m_isRenderable; }

	// Query properties of the FBo
	inline bool HasDepth( void ) const              { return m_depthBufs.Size() > 0u; }
	inline bool HasStencil( void ) const            { return m_stencilBufs.Size() > 0u; }
	inline uint NumColorAttachments( void ) const   { return m_colorBufs.Size(); }

	// Returns an internal renderable texture.
	IGLURenderTexture &operator[]   ( int buffer );
	IGLURenderTexture *GetAttachment( int buffer );

	// Get the width & height of the FBO.  Note, this may not be the w & h of all (or any)
	//    of the individual attached render targets, since the renderable area is the
	//    union of all the attachement sizes.
	inline int GetRenderableWidth( void ) const     { return m_width; }
	inline int GetRenderableHeight( void ) const    { return m_height; }

	// Attach a renderable texture to this buffer
	bool AttachTexture( int buffer, IGLURenderTexture *tex );

	//Resize attachments of this framebuffer
	void Resize( int width, int height );

	// A pointer to a IGLUFramebuffer could have type IGLUFramebuffer::Ptr
	typedef IGLUPtr<IGLUFramebuffer,IGLURenderTexture&,IGLURenderTexture&> Ptr;

	////////////////////////////////////////////////////////////////////////////////////
	// Helper methods to create useful callbacks for IGLUFramebuffers.
	//

	void Callback_UniformResize( IGLUVariable *intSize );  // Do not call directly.  Resizes to intSize x intSize
	IGLUCallback *GetResizeCallback( void ) { return new IGLUMemberCallback<IGLUFramebuffer>( *this, &IGLUFramebuffer::Callback_UniformResize ); }

	////////////////////////////////////////////////////////////////////////////////////
	// Some advanced usage methods.  Don't use these unless you actually know what you're doing!
	//

	// This should only need to be called if you *manually* create a MRT IGLUFramebuffer.
	//    It tells IGLU to try to call glDrawBuffers() appropriately when calling Bind().
	// The automatic call to glDrawBuffers() is somewhat simplistic, so if you want to draw to
	//    the textures on strange attachment combinations (e.g., buffers #0, 3, and 5), you will
	//    need to do this manually by calling glDrawBuffers() yourself.
	inline void EnableAutomaticSetDrawbuffers( void )     { m_setDrawBufsOnBind = true; }



	////////////////////////////////////////////////////////////////////////////////////
	// Static helper methods for constructing more complex IGLUFramebuffers 
	//

	// Create a single color-plane with a specified internal format, and (optionally) a depth plane 
	//    and/or stencil plane.  If both depth & stencil are enabled, it creates a depth-stencil texture.
	static IGLUFramebuffer::Ptr Create( GLenum internalFormat, int width, int height, bool hasDepth=false, bool hasStencil=false );

	// Create a FBO with multiple color planes, and optionally a depth plane 
	static IGLUFramebuffer::Ptr CreateMRT( GLenum internalFormat, int width, int height, int numColorTargets, bool hasDepth=false );

	// Create a FBO used for shadow mapping.  This guarantees that the buffer has an appropriate depth format.
	//    It may also create a color buffer (it currently does so), since many OSes require a color buffer to
	//    be "framebuffer complete".  However, it will setup the FBO class such that upon Bind(), a color mask
	//    will be set to disable rendering to the color buffer.  The color mask will be re-enabled when Unbind()ing.
	//    Additionally, the interpolation mode on the depth buffers are set to GL_NEAREST.
	// The type parameter specifies the type of shadow map.  Default is IGLU_FBO_BASIC_SHADOW_MAP.
	//    The addlNum parameter is ignored for IGLU_FBO_BASIC_SHADOW_MAP, specifies the number of layers in
	//    a IGLU_FBO_ARRAY_SHADOW_MAP, and the number of samples in a IGLU_FBO_MULTISAMPLE_SHADOW_MAP
	static IGLUFramebuffer::Ptr CreateShadowMap( int width, int height, uint type=0xFFFFFFFEu, int addlNum=1 );

	// Create a single color-plane with multisampling and aspecified internal format, and (optionally) a 
	//    depth and/or stencil plane.  If both depth & stencil are enabled, it creates a depth-stencil texture.
	static IGLUFramebuffer::Ptr CreateMultisample( GLenum internalFormat, int width, int height, 
		                                           int numSamples, bool hasDepth=false, bool hasStencil=false );

	// Create a single color render target consisting of a texture array with a specified internal format, and 
	//    (optionally) a depth plane and/or stencil texture array.  If both depth & stencil are enabled, 
	//    it creates a depth-stencil texture.
	static IGLUFramebuffer::Ptr CreateArray( GLenum internalFormat, int width, int height, int layers,
	                                         bool hasDepth=false, bool hasStencil=false );

protected:
	GLuint m_clearBits;
	bool   m_isBound;
	bool   m_isRenderable;
	bool   m_isShadowMapFBO;    // Special flag which sets/unsets a colormask inside Bind()/Unbind()
	bool   m_setDrawBufsOnBind; // Special flat which sets/unsets the draw buffers inside Bind()/Unbind()
	GLuint m_fboID;
	int    m_maxColorBufs;
	GLuint m_fboTextureType;
	GLint  m_prevBinding;
	int    m_width, m_height;
	float  m_prevViewport[4];

	IGLUArray1D<IGLURenderTexture *> m_colorBufs;
	IGLUArray1D<IGLURenderTexture *> m_depthBufs;   // Currently, should be exactly 0 or 1
	IGLUArray1D<IGLURenderTexture *> m_stencilBufs; // Currently, should be exactly 0 or 1

	// Some internal utility/helper functions
	bool CheckTextureType( GLuint texType );
	void ErrorTextureTypeMismatch( GLuint inputType );
};

enum {
	IGLU_COLOR         =  0,
	IGLU_COLOR0        =  0,
	IGLU_COLOR1        =  1,
	IGLU_COLOR2        =  2,
	IGLU_COLOR3        =  3,
	IGLU_COLOR4        =  4,
	IGLU_COLOR5        =  5,
	IGLU_COLOR6        =  6,
	IGLU_COLOR7        =  7,
	IGLU_COLOR_MAX     =  8,
	IGLU_DEPTH         = 50,
	IGLU_STENCIL       = 60,
	IGLU_DEPTH_STENCIL = 70
};

enum {
	// This can be passed as an 'internalFormat' to Create(), CreateArray(), and 
	//     CreateMultisample() to create a z-buffer (shadow map) only FBO.  It 
	//     also sets the FBO up to automatically set the texture appropriately 
	//     when Bind()ing it so it can be used in a 'texture*Shadow' sampler in GLSL.
	IGLU_FBO_SHADOW_MAP = 0xFFFFFFFFu,

	// These can be used as the optional 3rd parameter to CreateShadowMap().  They're
	//     mostly meant for internal IGLU use, but you can use them directly.
	IGLU_FBO_BASIC_SHADOW_MAP       = 0xFFFFFFFEu,
	IGLU_FBO_ARRAY_SHADOW_MAP       = 0xFFFFFFFDu,
	IGLU_FBO_MULTISAMPLE_SHADOW_MAP = 0xFFFFFFFCu,
};

// End namespace iglu
}


#endif