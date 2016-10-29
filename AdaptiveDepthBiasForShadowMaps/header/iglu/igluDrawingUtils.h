/***************************************************************************/
/* igluDrawingUtils.h                                                      */
/* -------------------                                                     */
/*                                                                         */
/* This "class" is simply a collection of drawing utilities, each is a     */
/*     static method.  In some sense, these really belong as part of other */
/*     classes (e.g., the full-screen texture draw could belong to         */
/*     IGLUTexture).  But these are not really integral to the behavior    */
/*     of the associated classes and carry a lot of overhead (e.g.,        */
/*     shaders) that would clutter up other classes.                       */
/*                                                                         */
/* Chris Wyman (10/13/2011)                                                */
/***************************************************************************/

#ifndef IGLU_DRAWINGUTILS_H
#define IGLU_DRAWINGUTILS_H


#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>

namespace iglu {

class IGLUTexture;
class vec4;

class IGLUDraw
{
private:
	IGLUDraw() {};
	~IGLUDraw() {};

public:
	// Draw full screen textures
	static int Fullscreen( const IGLUTexture &texture, uint flags=0 );		
	static int Fullscreen(       IGLUTexture *texture, uint flags=0 );

	// Draws a full screen texture using a shader specified by the developer.  
	//    'texVar' specifies the name of the texture uniform in the shader.
	// The quad that is drawn has vertex positions (in NDC, i.e., [-1..1]),
	//     a default normal (0,0,-1), and texture coords (ranging from [0..1])
	//     at the standard IGLU attribute locations (#0, 1, and 2, respectively)
	static int Fullscreen( IGLUShaderProgram::Ptr &shader, const IGLUTexture &texture, const char *texVar="inputTex" );
	static int Fullscreen( IGLUShaderProgram::Ptr &shader,       IGLUTexture *texture, const char *texVar="inputTex" );

	// Draw a texture on the screen at some location
	static int Texture( const IGLUTexture &texture, float left, float bottom, float right, float top, uint flags=0 );
	static int Texture(       IGLUTexture *texture, float left, float bottom, float right, float top, uint flags=0 );

	// Draw text onto the screen (in a particular viewport).  The viewport is queried for positioning.
	static int DrawText( int igluFont, int x_left, int y_bottom, const char *text, float scaleFactor=1.0, int targetViewport=0 );
	static int DrawColorText( int igluFont, int x_left, int y_bottom, const vec4 &color, const char *text, float scaleFactor=1.0, int targetViewport=0 );

};

enum {
	IGLU_DRAW_NONE        = 0x0000,
	IGLU_DRAW_FLIP_X      = 0x0001,
	IGLU_DRAW_FLIP_Y      = 0x0002,
	IGLU_FONT_FIXED       = 0x0004,
	IGLU_FONT_VARIABLE    = 0x0008,

	IGLU_DRAW_TEX_BLACK_AS_TRANSPARENT = 0x0010,
	IGLU_DRAW_TEX_SIMPLE_BRIGHTEN      = 0x0020,
};


// End iglu namespace
}




#endif

