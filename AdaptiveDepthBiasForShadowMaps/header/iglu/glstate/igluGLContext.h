/******************************************************************/
/* igluGLContext.h                                                */
/* -----------------------                                        */
/*                                                                */
/* A class encapsulating basic information and settings for the   */
/*    current OpenGL context.                                     */
/*                                                                */
/*  --> Currently, everything only works for a single GL context  */
/*                                                                */
/* Chris Wyman (06/13/2011)                                       */
/******************************************************************/


#ifndef IGLU_GL_CONTEXT_H
#define IGLU_GL_CONTEXT_H

namespace iglu { 

class IGLUContext {
private:
	IGLUContext()  {}
	~IGLUContext() {}

	static int  m_glAttribList[128];
	static char m_glVerStr[256];
public:
	// Some enum'd bitflags
	enum {
		FLAGS_NONE            = 0x0000,
		FLAGS_DEBUG           = 0x0001,
		FLAGS_FORWARD_COMPAT  = 0x0002   // Generally, avoid this.
	};

	enum {
		PROFILE_CORE          = 0x0001,
		PROFILE_COMPATIBILITY = 0x0002
	};


	// Prior to creating an OpenGL context, you may set these parameters
	//   which will control which types of OpenGL context is created
	static int  glMajorVer, glMinorVer;
	static uint glProfileBits, glFlags;

	// A helper to set the values above
	static void SetGLVersion( int majorVer, int minorVer, uint profileBits=PROFILE_COMPATIBILITY );

	// Gives a string that describes the current OpenGL context (if initialized)
	static char *GetGLVersion( void );

	// Given the current settings of glMajorVer, glMinorVer, glProfileBits, and glFlags,
	//   create (and return a pointer to) an attriblist that can be passed to either
	//   wglCreateContextAttribs() or glXCreateContextAttribs()
	static int  *GetContextAttribList( void );

	// Returns true/false if the current context has certain OpenGL properties
	static bool  IsCoreContext( void );
	static bool  IsCompatibilityContext( void );
	static bool  IsDebugContext( void );

};

}  // end namespace iglu {



#endif