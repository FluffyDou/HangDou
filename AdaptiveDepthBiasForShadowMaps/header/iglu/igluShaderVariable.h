/**************************************************************************
** igluShaderVariable.h                                                  **
** -----------------                                                     **
**                                                                       **
** This header points to helper classes that should not need to be       **
**   manipulated directly by the user...  These are used as part of the  **
**   functionality of the IGLUShaderProgram class.  In particular the [] **
**   operator allows you to select shader variables to set via the =     **
**   operator.  There is an intermediate IGLUShaderVariable class that   **
**   actually catches the assignment operator.  Similarly, there's a     **
**   IGLUShaderTexture class that is used to remember textures between   **
**   shader uses so they can be correctly enabled and disabled.          **
**                                                                       **
** Chris Wyman (09/27/2011)                                              **
**************************************************************************/

#ifndef __IGLU_SHADER_VARIABLE_H
#define __IGLU_SHADER_VARIABLE_H

// Other needed headers
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>

//#include "igluShaderProgram.h"

namespace iglu {

class IGLUMatrix4x4;
class IGLUTexture;
class IGLUShaderProgram;
class IGLUFloat;
class IGLUInt;
class IGLUBool;

class IGLUShaderVariable
{
public:
	// Default constructor
	IGLUShaderVariable() : m_varIdx(-1), m_parent(0) {}

	// Constructor for the shader encapsulation
	IGLUShaderVariable( GLuint programID, const char *varName, IGLUShaderProgram *parent );
	IGLUShaderVariable( const IGLUShaderVariable &copy );

	// Queries to this shader variable
	inline bool IsUniform( void ) const	            { return (m_varIdx >= 0) && !m_isAttribute; }
	inline bool IsAttribute( void ) const           { return (m_varIdx >= 0) && m_isAttribute; }

	// Query the uniform/attribute location
	inline GLuint GetVariableIndex( void ) const    { return m_varIdx; }

	// Assignment operators, to set the uniform value in the shader 

	// From a single scalar
	void operator= (         float val );
	void operator= (        double val );
	void operator= (           int val );
	void operator= (  unsigned int val );
	void operator= (    IGLUFloat *val );  // So they need not be dereferenced
	void operator= (      IGLUInt *val );  // So they need not be dereferenced

	// From a 2-dimensional vector
	void operator= (        float2 val );  // Also takes fvec2 or vec2
	void operator= (       double2 val );  // Also takes dvec2
	void operator= (          int2 val );  // Also takes ivec2
	void operator= (         uint2 val );  // Also takes uivec2

	// From a 3-dimensional vector
	void operator= (        float3 val );  // Also takes fvec3 or vec3
	void operator= (       double3 val );  // Also takes dvec3
	void operator= (          int3 val );  // Also takes ivec3
	void operator= (         uint3 val );  // Also takes uivec3

	// From a 4-dimensional vector
	void operator= (        float4 val );  // Also takes fvec4 or vec4
	void operator= (       double4 val );  // Also takes dvec4
	void operator= (          int4 val );  // Also takes ivec4
	void operator= (         uint4 val );  // Also takes uivec4

	// From a 4x4 matrix
	void operator= ( const IGLUMatrix4x4 &val );
	void operator= ( const IGLUMatrix4x4 *val ); // So they need not be dereferenced

	// Texture assignment
	void operator= (   const IGLUTexture &val );  
	void operator= (   const IGLUTexture *val );  

private:
	GLsizei m_varNameLen;
	GLint   m_varSize;
	GLenum  m_varType;
	char    m_varName[32];
	GLint   m_varIdx;
	char    m_errorVal[64];
	bool    m_isAttribute;
	IGLUShaderProgram *m_parent;

	// Function to print out various user errors detected at runtime
	bool TypeMismatch( const char *type );                 // Always returns true
	bool UnknownVariable( const char *val );               // Always returns true
	bool AssignmentToAttribute( const char *type );        // Always returns true
	void TypeToChar( GLenum varType, char *retTypeName );

	// Type checking between textures in GLSL and input IGLUTexture is tricky 
	//    and is used in at least two different assignment operators.  Separate
	//    this code out so it's easier to incrementally update, fix, and doesn't
	//    muck of the legibility of the assignment operators.
	bool HasTextureTypeMismatch( const IGLUTexture *rhsTex );
	bool HasImageTypeMismatch( const IGLUTexture *rhsTex );

	// If we have an error, we want to print out the names of the types we are assigning.
	//    The GLSL type names are complex to compute, split this off to a method so we
	//    only do it when a type mismatch is detected.
	void GetRHSTextureTypeName( const IGLUTexture *rhsTex, int dataType, bool isDepth, char *typeName );

	// Utilities to query properties about GLSL texture uniforms
	static bool Is1DSampler( GLenum uniformType );
	static bool Is2DSampler( GLenum uniformType );
	static bool Is3DSampler( GLenum uniformType );
	static bool IsArraySampler( GLenum uniformType );
	static bool IsShadowSampler( GLenum uniformType );
	static bool IsCubeSampler( GLenum uniformType );
	static bool IsIntSampler( GLenum uniformType );
	static bool IsUIntSampler( GLenum uniformType );
	static bool IsRectSampler( GLenum uniformType );
	static bool IsBufferSampler( GLenum uniformType );
	static bool IsMultiSampleSampler( GLenum uniformType );

	static bool IsImage(GLenum uniformType);
	static bool Is1DImage( GLenum uniformType );
	static bool Is2DImage( GLenum uniformType );
	static bool Is3DImage( GLenum uniformType );
	static bool IsArrayImage( GLenum uniformType );
	static bool IsShadowImage( GLenum uniformType );
	static bool IsCubeImage( GLenum uniformType );
	static bool IsIntImage( GLenum uniformType );
	static bool IsUIntImage( GLenum uniformType );
	static bool IsRectImage( GLenum uniformType );
	static bool IsBufferImage( GLenum uniformType );
	static bool IsMultiSampleImage( GLenum uniformType );
};


// This is a helper class that is used by IGLUShaderProgram to remember information
//    about textures previously assigned to this shader so they can be enabled and
//    disabled appropriately.
class IGLUShaderTexture {
public:
	IGLUShaderTexture() : m_tex(0), m_glslVarLoc(-1) {}
	~IGLUShaderTexture() {} 

	char   m_texName[64];
	GLint  m_glslVarLoc;
	GLint  m_texUnit;
	const IGLUTexture *m_tex;
};


// End namespace iglu
}

#endif